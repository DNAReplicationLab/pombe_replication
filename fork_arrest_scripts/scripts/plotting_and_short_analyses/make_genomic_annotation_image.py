import argparse
import drawsvg
import json
import os
import pandas as pd
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from itertools import repeat


def check_and_adjust_import_intervals(genomic_data: dict) -> tuple[str, int, int]:
    """ Check and adjust the various import intervals if needed.
    The user may ask for plotting in some interval e.g.: chr1:1000-2000, but may ask for
    imports from bed files in other regions e.g. chr2:500-1000, which is not consistent.
    Or they may ask for plotting in some interval e.g.: chr1:1000-2000, but may ask for imports
    from a larger region e.g. chr1:500-3000, which is also not consistent. This function checks
    for consistency, complaining if the problems cannot be fixed or adjusting the import intervals
    if they can be fixed.
    We specifically look for:
        - a dictionary with type set to "reference_line" and associated start, stop (mandatory) and
          contig (optional) keys. These tell us the genomic region we are plotting in.
        - any other dictionary with an import_from_bed key. This key must be a list of two elements:
          the bed file name and a bed region in the format contig:start-end or just contig.
        After gathering this information, we check for consistency between these intervals
        as described above, adjusting the import intervals if needed.

    Args:
        genomic_data: a dictionary with a key "genomic_elements" that contains a list of dictionaries,
            each of which represents a genomic element like an origin, gene etc. We adjust this in place
            as described above.

    Returns:
        tuple of contig, start, stop

    """
    if "genomic_elements" not in genomic_data or not isinstance(genomic_data["genomic_elements"], list):
        raise ValueError("Invalid genomic_data format: missing 'genomic_elements' list.")

    reference_line = list(filter(lambda x: isinstance(x, dict) and x.get("type") == "reference_line",
                                 genomic_data["genomic_elements"]))

    if len(reference_line) != 1:
        raise ValueError("At least one and at most one reference_line element is needed!")
    else:
        check_contig = reference_line[0].get("contig", "")
        check_start = int(reference_line[0].get("start"))
        check_stop = int(reference_line[0].get("stop"))
        if check_start >= check_stop:
            raise ValueError("Reference line start must be less than stop!")

    # Parse genomic elements
    for element in (y for y in genomic_data["genomic_elements"] if isinstance(y, dict)
                    and y.get("type") != "reference_line" and "import_from_bed" in y):
        bed_info = element["import_from_bed"]
        if not (isinstance(bed_info, list) and (len(bed_info) >= 2)):
            raise ValueError("Malformed import_from_bed element!")
        else:
            bed_region = bed_info[1]
            if ":" in bed_region and "-" in bed_region:
                bed_contig, range_part = bed_region.split(":")
                bed_start, bed_stop = map(int, range_part.split("-"))
                if bed_start >= check_stop or bed_stop <= check_start:
                    bed_start = bed_stop = 0
                else:
                    bed_start = max(bed_start, check_start)
                    bed_stop = min(bed_stop, check_stop)
            else:
                bed_contig = bed_region
                bed_start = check_start
                bed_stop = check_stop

            if bed_contig == "":
                raise ValueError("Contig must be specified in import_from_bed!")
            elif check_contig != bed_contig:
                # if no contig was specified in reference_line, then set it to the contig in import_from_bed
                # otherwise, complain about the mismatch
                if check_contig == "":
                    check_contig = bed_contig
                else:
                    raise ValueError("Contig mismatch between reference line and import_from_bed!")
            else:
                element["import_from_bed"][1] = f"{bed_contig}:{bed_start}-{bed_stop}"

    return check_contig, check_start, check_stop


def import_from_bed_if_needed(genomic_element: dict) -> dict:
    """ If the genomic element contains an import_from_bed key, then import genomic elements from bed file.
    The import_from_bed key must be a list of 2-3 elements: the bed file name, a bed region in the format
    contig:start-end or just contig, and an optional third element, an integer interpreted as minimum score
    to filter the bed file. Only elements that overlap (even partially) with the bed region will be imported.
    The bed file must have at least three columns and must be in the standard format. To use score-based filtering,
    we need at least five columns in the bed file.

    For instance,
    { "type": "origin", "import_from_bed": ["/a/b/c/d.bed", "chrI:1000-2000"], "modifier": "inactive"}
    might get converted to something like
    [{ "type": "origin", "start": 1100, "stop": 1200, "modifier": "inactive", label: "origin1"},
    { "type": "origin", "start": 1400, "stop": 1300, "modifier": "inactive", label: "origin2"}]
    (start < end or start > end depend on the orientation of the bed region, and labels correspond
    to the name column of the bed file if it exists, or is set to None).

    Args:
        genomic_element: Genomic elements like origins etc.

    Returns:
        the same genomic element but with the import_from_bed key processed if it exists
    """

    if "import_from_bed" in genomic_element:

        if len(genomic_element["import_from_bed"]) < 2:
            raise ValueError("Malformed import_from_bed element!")
        elif len(genomic_element["import_from_bed"]) == 2:
            bed_file, bed_region = genomic_element["import_from_bed"]
            min_score = None
        elif len(genomic_element["import_from_bed"]) == 3:
            bed_file, bed_region, min_score = genomic_element["import_from_bed"]
        else:
            raise ValueError("Malformed import_from_bed element!")

        if ":" in bed_region and "-" in bed_region:
            contig, bed_region = bed_region.split(":")
            start, end = bed_region.split("-")
            start, end = int(start), int(end)
        else:
            contig = bed_region
            start = 0
            end = float("inf")

        # ensure start < end
        if start > end:
            raise ValueError("Start must be less than end")
        elif start == end:
            return []

        item = dict(genomic_element)
        item["start"] = []
        item["stop"] = []
        item["label"] = []

        with open(bed_file, 'r') as bed_file_pointer:

            for line in filter(lambda x: not x.startswith(("#", "track", "browser")), bed_file_pointer):

                line = line.strip().split()

                current_contig = line[0]
                current_start = int(line[1])
                current_end = int(line[2])
                current_name = line[3] if len(line) > 3 else None
                current_score = float(line[4]) if len(line) > 4 else None
                current_strand = line[5] if len(line) > 5 else "+"

                if min_score is not None and current_score is not None and current_score < min_score:
                    continue

                if current_strand not in ["+", "-"]:
                    raise ValueError("Strand must be + or -")

                if (current_start < current_end and current_contig == contig and
                        current_end > start and current_start < end):

                    item_start = current_start
                    item_end = current_end

                    # flip start and stop if the strand is negative
                    if current_strand == "-":
                        item_start, item_end = item_end, item_start

                    item["start"].append(item_start)
                    item["stop"].append(item_end)
                    item["label"].append(current_name)

        return item
    else:
        return genomic_element


def import_from_bedgraphs(bedgraphs_data: list[dict], ref_contig: str, ref_start: int, ref_end: int) -> list[dict]:
    """
    Import data from bedgraph files and process it into a format suitable for plotting.

    Args:
        bedgraphs_data: List of dictionaries containing bedgraph file information.
                        Each dictionary should have the key "bedgraph" containing the path to the bedgraph file.
        ref_contig: Reference contig to filter the bedgraph data.
        ref_start: Start position of the reference range to filter the bedgraph data.
        ref_end: End position of the reference range to filter the bedgraph data.

    Returns:
        Same list of dictionaries as input, but the "bedgraph" key has been replaced with:
            - "x": List of midpoints of the start and end positions from the bedgraph data.
            - "y": List of values from the bedgraph data.
            - "start": List of start positions from the bedgraph data.
            - "end": List of end positions from the bedgraph data.

    Raises:
        FileNotFoundError: If any of the specified bedgraph files do not exist.
    """
    processed_data = []

    # iterate through each entry in bedgraphs_data
    for entry in bedgraphs_data:

        bedgraph_file = entry.pop("bedgraph")

        # complain if the bedgraph file does not exist
        if not os.path.exists(bedgraph_file):
            raise FileNotFoundError(f"Bedgraph file {bedgraph_file} not found!")

        # Read the bedgraph file into a DataFrame
        df = pd.read_csv(bedgraph_file, sep=' ', header=None, names=["contig", "start", "end", "value"], comment="#")

        # Filter the DataFrame for the specified contig and range
        df_filtered = df[(df["contig"] == ref_contig) & (df["start"] >= ref_start) & (df["end"] <= ref_end)].copy()

        # Calculate the midpoint of start and end
        df_filtered.loc[:, "x"] = (df_filtered["start"] + df_filtered["end"]) / 2

        # Extract the x, y, start, and end values
        x = df_filtered["x"].tolist()
        y = df_filtered["value"].tolist()
        start = df_filtered["start"].tolist()
        end = df_filtered["end"].tolist()

        # Create the processed entry
        processed_entry = {
            "x": x,
            "y": y,
            "start": start,
            "end": end
        }

        # Retain any other keys from the original entry
        for key, value in entry.items():
            processed_entry[key] = value

        # Append the processed entry to the processed data
        processed_data.append(processed_entry)

    return processed_data


def convert_genomic_lists_if_needed(genomic_element: dict) -> list[dict]:
    """ If the genomic element is a list element, then convert it into list of genomic elements.

    For instance,
    { "type": "origin", "position": [7.625, 9.5], "modifier": "inactive"}
    will get converted to
    [{ "type": "origin", "position": 7.625, "modifier": "inactive"},
    { "type": "origin", "position": 9.5, "modifier": "inactive"}]

    and
    { "type": "origin", "position": 10, "modifier": "inactive"}
    gets converted to
    [{ "type": "origin", "position": 10, "modifier": "inactive"}]

    Args:
        genomic_element: Genomic elements like origins etc.

    Returns:
        list of genomic elements
    """
    output_list = []

    if "position" in genomic_element and isinstance(genomic_element["position"], list):
        for k in zip(genomic_element["position"],
                     genomic_element["label"] if "label" in genomic_element else repeat(None)):
            item = dict(genomic_element)
            item["position"] = k[0]
            item["label"] = k[1]
            output_list.append(item)
    elif "start" in genomic_element and isinstance(genomic_element["start"], list):
        for k in zip(genomic_element["start"], genomic_element["stop"],
                     genomic_element["label"] if "label" in genomic_element else repeat(None)):
            item = dict(genomic_element)
            item["start"] = k[0]
            item["stop"] = k[1]
            item["label"] = k[2]
            output_list.append(item)
    else:
        output_list.append(genomic_element)

    return output_list


def process_genomic_elements_into_shapes(genomic_element: dict, scale: float,
                                         no_stroke: bool, use_dna_features_viewer: bool = False) -> list[dict]:
    """ Convert genetic elements into shapes to be used in annotations

    Args:
        genomic_element: elements like origin, fork stall site etc.
        scale: scale factor associated with shape
        no_stroke: if True, then no stroke or boundaries in svg image. Default is False.
        use_dna_features_viewer: if True, then use the DNA Features Viewer library to draw the image. Default is False.

    Returns:
        list with elements corresponding to graphical objects that can be used to make pictures
    """

    # some elements below use the "start", "stop" notation, whereas others use "position", "orientation" notation.
    # we convert "start" "stop" notation to "position" "orientation" notation if the former is present,
    # and check for consistency if both styles are present.
    if "start" in genomic_element and "stop" in genomic_element:

        if genomic_element["start"] == genomic_element["stop"]:
            raise ValueError("Genomic element has zero length!")

        position = (genomic_element["start"] + genomic_element["stop"]) / 2
        orientation = "right" if genomic_element["start"] < genomic_element["stop"] else "left"

        # ensure genomic position and orientation if specified match the calculated values above.
        # if not specified, then set them to the calculated values
        genomic_element["position"] = genomic_element.get("position", position)
        genomic_element["orientation"] = genomic_element.get("orientation", orientation)

        if genomic_element["position"] != position:
            raise ValueError("Genomic element has inconsistent position!")

        if genomic_element["orientation"] != orientation:
            raise ValueError("Genomic element has inconsistent orientation!")

    if genomic_element["type"] == "origin":

        # set characteristics depending on origin active or not
        # we assume active by default
        color_options = {"active": "#ffff00", "inactive": "#fcfcf2"}
        color = color_options[genomic_element.get("modifier", "active")]

        if use_dna_features_viewer:
            return [GraphicFeature(
                start=genomic_element["position"] - scale,
                end=genomic_element["position"] + scale,
                color=color,
                thickness=0.4,
                linewidth=1.0 if not no_stroke else 0,
                linecolor="#000000",
                type="origin"
            )]
        else:
            return [
                {
                    "type": "circle",
                    "center": genomic_element["position"],
                    "radius": scale,
                    "color": color,
                    "stroke": "black",
                    "stroke_width": scale / 6 if not no_stroke else 0
                }
            ]

    elif genomic_element["type"] == "gene":

        # set characteristics depending on whether gene is marked as "emphasize" or "de-emphasize"
        color_options = {"emphasize": "#007a3d", "de-emphasize": "#9bffcd"}
        color = color_options[genomic_element.get("modifier", "emphasize")]

        # fix start, stop if orientation is not consistent with it
        if ("orientation" in genomic_element and
                ((genomic_element["orientation"] in ["+", "right", "fwd"] and
                  genomic_element["start"] > genomic_element["stop"]) or
                 (genomic_element["orientation"] in ["-", "left", "rev"] and
                  genomic_element["start"] < genomic_element["stop"]))):
            genomic_element["start"], genomic_element["stop"] = genomic_element["stop"], genomic_element["start"]

        flag = -1 if genomic_element["start"] < genomic_element["stop"] else 1

        if use_dna_features_viewer:
            if genomic_element.get("representation", "on_line") == "displaced":
                print("Warning: Displaced representation not supported in DNA Features Viewer mode!")
            return [GraphicFeature(
                start=min(genomic_element["start"], genomic_element["stop"]),
                end=max(genomic_element["start"], genomic_element["stop"]),
                strand=-flag if genomic_element.get("remove_strand", False) is False else 0,
                color=color,
                linewidth=1.0 if not no_stroke else 0,
                linecolor="#000000",
                label=genomic_element.get("label", None) if genomic_element.get("remove_label", False) is False else None
            )]
        else:
            # set x values depending on representation chosen
            if genomic_element.get("representation", "on_line") == "on_line":
                arrow_x_vals = [genomic_element["start"], genomic_element["stop"] + flag * scale,
                                genomic_element["stop"]]
                if not (arrow_x_vals[0] < arrow_x_vals[1] < arrow_x_vals[2] or
                        arrow_x_vals[0] > arrow_x_vals[1] > arrow_x_vals[2]):
                    arrow_x_vals[1] = arrow_x_vals[0]
            else:
                arrow_x_vals = [genomic_element["start"], genomic_element["stop"], genomic_element["stop"]]

            # set y levels depending on representation chosen
            # we assume on_line representation by default
            arrow_y_options = {"on_line": [scale/2, 0, -scale/2],
                               "displaced": [0, -scale/2, -scale] if flag == -1 else [scale, scale/2, 0]}
            arrow_y_vals = arrow_y_options[genomic_element.get("representation", "on_line")]

            # some elements are not supported in this mode, so print warnings
            if genomic_element.get("label", None) is not None:
                print("Warning: Label not supported in this mode!")
            if genomic_element.get("remove_label", False):
                print("Warning: remove_label not supported in this mode; no labels are drawn anyway!")
            if genomic_element.get("remove_strand", False):
                print("Warning: remove_strand not supported in this mode!")

            return [
                {
                    "type": "polygon",
                    "vertex_positions": [arrow_x_vals[k] for k in [0, 1, 2, 1, 0]],
                    "vertex_offsets": [arrow_y_vals[k] for k in [2, 2, 1, 0, 0]],
                    "color": color,
                    "stroke": "black",
                    "stroke_width": scale / 6 if not no_stroke else 0
                }
            ]

    elif genomic_element["type"] == "reference_line":

        if use_dna_features_viewer:
            raise NotImplementedError("Reference line not supported in DNA Features Viewer mode!")
        else:
            return [
                {
                    "type": "rectangle",
                    "start": genomic_element["start"],
                    "end": genomic_element["stop"],
                    "thickness": scale / 3,
                    "color": "#000000",
                    "stroke": "black",
                    "stroke_width": 0
                }
            ]

    elif genomic_element["type"] == "fork_stall":

        # set characteristics depending on whether gene is marked as "emphasize" or "de-emphasize"
        color_options = {"emphasize": "#444444", "de-emphasize": "#dddddd"}
        color = color_options[genomic_element.get("modifier", "emphasize")]

        # get the orientation of stall site and draw
        flag_options = {"left": -1, "right": 1}
        flag = flag_options[genomic_element["orientation"]]

        if use_dna_features_viewer:
            return [GraphicFeature(
                start=genomic_element["position"] - scale,
                end=genomic_element["position"] + scale,
                strand=-flag,
                color=color,
                thickness=0.4,
                linewidth=1.0 if not no_stroke else 0,
                linecolor="#000000",
                type="fork_stall"
            )]
        else:
            shape_x_vals = [genomic_element["position"], genomic_element["position"] + flag * scale,
                            genomic_element["position"] - flag * scale / 2]
            shape_y_vals = [2 * flag * scale, 0, -2 * flag * scale]

            return [
                {
                    "type": "polygon",
                    "vertex_positions": [shape_x_vals[k] for k in [0, 1, 2, 2, 1]],
                    "vertex_offsets": [shape_y_vals[k] for k in [1, 2, 2, 0, 0]],
                    "color": color,
                    "stroke": "black",
                    "stroke_width": scale / 6 if not no_stroke else 0
                }
            ]

    elif genomic_element["type"] == "model_pause":

        if use_dna_features_viewer:

            # set characteristics depending on whether pauses are marked lead or lag.
            color_options = {"default": "#ff0000", "lead": "#000000", "lag": "#ffffff"}
            color = color_options[genomic_element.get("modifier", "default")]

            return [GraphicFeature(
                start=genomic_element["position"] - scale,
                end=genomic_element["position"] + scale,
                strand=+1 if genomic_element.get("orientation", "right") == "right" else -1,
                color=color,
                thickness=0.4,
                linewidth=1.0 if not no_stroke else 0,
                linecolor="#000000",
                type="model_pause"
            )]
        else:

            if genomic_element.get("modifier", "default") not in ["default"]:
                raise NotImplementedError("Model pause modifier not supported in this mode!")

            # 2.83 below is roughly sqrt(8). we have chosen this so that the lengths of the sides of the triangle
            # representing the model pause are 3 scl, 3 scl and 2 scl, with the shorter side lying parallel
            # to the reference line.
            shape_x_vals = [genomic_element["position"], genomic_element["position"] - scale,
                            genomic_element["position"] + scale]
            if genomic_element.get("orientation", "right") == "left":
                shape_y_vals = [2.83 * scale, 0]
            else:
                shape_y_vals = [-2.83 * scale, 0]

            return [
                {
                    "type": "polygon",
                    "vertex_positions": shape_x_vals,
                    "vertex_offsets": [shape_y_vals[k] for k in [1, 0, 0]],
                    "color": "red",
                    "stroke": "black",
                    "stroke_width": scale/6 if not no_stroke else 0
                }
            ]
    else:
        raise ValueError("Element unsupported!")


def draw_with_dna_features_viewer(genomic_data: dict, filename_svg: str, no_stroke: bool) -> None:
    """ Use the DNA Features Viewer library to draw the image.

    Args:
        filename_svg: output svg filename
        genomic_data: contains annotations used in genomic plot.
        no_stroke: if True, then no stroke or boundaries in svg image. Default is False.

    Returns:
        None
    """

    features = []
    genomic_start = genomic_stop = 0

    for k in genomic_data["genomic_elements"]:
        for item in convert_genomic_lists_if_needed(import_from_bed_if_needed(k)):
            if item['type'] == 'reference_line':
                genomic_start = item['start']
                genomic_stop = item['stop']
            else:
                features.extend(process_genomic_elements_into_shapes(item, genomic_data["scale_factor"],
                                                                     no_stroke, use_dna_features_viewer=True))

    if genomic_start == genomic_stop:
        raise ValueError("Genomic start and stop are the same!")
    else:
        graphic_record = GraphicRecord(
            sequence_length=max([k.end for k in features]),
            sequence=None,
            features=features
        )

    crop_start = genomic_start
    crop_stop = min(max([k.end for k in features]), genomic_stop)
    figure_width = genomic_data.get("img_width", 10)

    if genomic_data.get("bedgraph_data", []) == []:

        ax, _ = graphic_record.crop((crop_start, crop_stop)).plot(figure_width=figure_width)
        ax.figure.tight_layout()
        ax.figure.savefig(filename_svg)

    else:

        # if bedgraph data is present, then plot bedgraph data as well
        n_bedgraphs = len(genomic_data["bedgraph_data"])
        fig, bedgraph_axes = plt.subplots(
            n_bedgraphs + 1, 1, figsize=(figure_width, n_bedgraphs + 3),
            sharex=True, gridspec_kw={"height_ratios": [3] + [1] * n_bedgraphs}
        )
        graphic_record.crop((crop_start, crop_stop)).plot(figure_width=figure_width, ax=bedgraph_axes[0], with_ruler=False)

        # check that all bedgraph data have levels
        if not all("level" in k for k in genomic_data["bedgraph_data"]):
            raise ValueError("Bedgraph data must have levels!")
        
        # check that levels are a permutation of 1, 2, ..., n_bedgraphs
        if set(k["level"] for k in genomic_data["bedgraph_data"]) != set(range(1, n_bedgraphs + 1)):
            raise ValueError("Levels must be a permutation of 1, 2, ..., n_bedgraphs!")

        for bedgraph_data in genomic_data["bedgraph_data"]:
            level = bedgraph_data.get("level", -1)
            if level <= 0 or level > n_bedgraphs:
                raise ValueError("Invalid level!")
            else:
                x_values = []
                y_values = []

                if ("start" in bedgraph_data and "end" in bedgraph_data):
                    for counter in range(len(bedgraph_data["start"])):

                        if counter > 0 and not(x_values[-1] == bedgraph_data["start"][counter]):
                            x_values.append(x_values[-1])
                            y_values.append(0)
                            x_values.append(bedgraph_data["start"][counter])
                            y_values.append(0)

                        x_values.append(bedgraph_data["start"][counter])
                        x_values.append(bedgraph_data["end"][counter])
                        y_values.append(bedgraph_data["y"][counter])
                        y_values.append(bedgraph_data["y"][counter])
                            
                elif ("x" in bedgraph_data and "y" in bedgraph_data):
                    x_values = bedgraph_data["x"]
                    y_values = bedgraph_data["y"]
                else:
                    raise ValueError("Bedgraph data must have x and y values!")
                
                # check that x is in non-decreasing order
                if not all(x_values[i] <= x_values[i + 1] for i in range(len(x_values) - 1)):
                    raise ValueError("x values must be in non-decreasing order!")
                
                bedgraph_axes[level].fill_between(x_values, y_values, alpha=0.3)
                bedgraph_axes[level].set_ylabel(bedgraph_data.get("label", f"signal_{level}"))
                bedgraph_axes[level].set_ylim(bottom=0)
                bedgraph_axes[level].set_xlabel("Position (bp)")

        fig.tight_layout()
        fig.savefig(filename_svg)


def draw_genomic_annotation(width: float, height: float, pos1: tuple[float, float],
                            pos2: tuple[float, float], filename_svg: str,
                            genomic_data: dict, no_stroke: bool) -> None:
    """ Prepare a genomic-annotation-shapes image to accompany a genomic-data image

    Args:
        width: width of svg image in arbitrary units, equal to width of genomic-data image
        height: height of svg image in arbitrary units. As svgs are scalable, one of the lengths we specify here
            (width or height or some other length) is redundant. But we don't care.
        pos1: tuple of (actual genomic coordinate, image x coordinate it corresponds to)
        pos2: tuple of (actual genomic coordinate, image x coordinate it corresponds to)
        filename_svg: output svg filename
        genomic_data: contains annotations used in genomic plot.
        no_stroke: if True, then no stroke or boundaries in svg image. Default is False.

    Returns:
        None
    """

    # define transform from genomic to image coordinates
    scale_factor = (pos2[1] - pos1[1]) / (pos2[0] - pos1[0])

    def coord_transform(x):
        return pos1[1] + (x - pos1[0]) * scale_factor

    # initialize canvas
    d = drawsvg.Drawing(width, height, displayInline=False)

    # read plot elements file, iterate through each row and plot different shapes
    genomic_data["shapes"] = []

    for k in genomic_data["genomic_elements"]:
        for item in convert_genomic_lists_if_needed(import_from_bed_if_needed(k)):
            genomic_data["shapes"].extend(process_genomic_elements_into_shapes(
                item, genomic_data["scale_factor"], no_stroke))

    for row in genomic_data["shapes"]:

        if row["type"] == "rectangle":
            thickness = row["thickness"] * scale_factor
            d.append(drawsvg.Rectangle(coord_transform(row["start"]), height / 2 - thickness / 2,
                                       coord_transform(row["end"]) - coord_transform(row["start"]),
                                       thickness, fill=row["color"],
                                       stroke=row["stroke"],
                                       stroke_width=row["stroke_width"] * scale_factor))
        elif row["type"] == "circle":
            d.append(drawsvg.Circle(coord_transform(row["center"]), height / 2,
                                    row["radius"] * scale_factor, fill=row["color"],
                                    stroke=row["stroke"],
                                    stroke_width=row["stroke_width"] * scale_factor))
        elif row["type"] == "polygon":
            vertex_positions = [coord_transform(k) for k in row["vertex_positions"]]
            vertex_offsets = [height / 2 + k * scale_factor for k in row["vertex_offsets"]]
            d.append(drawsvg.Lines(*[val for pair in zip(vertex_positions, vertex_offsets) for val in pair],
                                   close=True,
                                   fill=row["color"],
                                   stroke=row["stroke"],
                                   stroke_width=row["stroke_width"] * scale_factor
                                   ))

    # save as svg
    d.save_svg(filename_svg)


if __name__ == "__main__":
    desc = """    Makes genomic annotation image.

        Sample usage:
            # standalone
            python programName.py --annotate annotate.json --op op.svg
        """

    # get options
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--annotate', type=str, required=True,
                        help='json file with genomic annotations')
    parser.add_argument('--op', type=str, required=True,
                        help='output svg file name')
    parser.add_argument('--no-stroke', action='store_true', required=False, default=False,
                        help='no stroke or boundaries in svg image')
    parser.add_argument('--use-dna-features-viewer', action='store_true', required=False,
                        default=False, help='use the DNA Features Viewer library to draw the image.')
    parser.add_argument('--bedgraphs', type=str, required=False,
                        help='json file with bedgraphs data')

    # parse arguments
    args = parser.parse_args()

    # load data from annotation file
    if not os.path.exists(args.annotate):
        raise FileNotFoundError(f"File {args.annotate} not found!")
    else:
        with open(args.annotate, 'r') as f:
            data = json.load(f)

    # check and adjust the various import intervals if needed
    ref_contig, ref_start, ref_end = check_and_adjust_import_intervals(data)

    # load bedgraphs data if provided
    if args.bedgraphs:
        if not (os.path.exists(args.bedgraphs) and args.use_dna_features_viewer and len(ref_contig) > 0):
            raise ValueError("Either the bedgraphs file does not exist, or the DNA Features Viewer mode "
                             "is not enabled, or the reference line is not specified in the annotation file!")
        else:
            with open(args.bedgraphs, 'r') as f:
                bedgraphs_data = json.load(f)
            data["bedgraph_data"] = import_from_bedgraphs(bedgraphs_data, ref_contig, ref_start, ref_end)

    # get the full path of the directory containing the annotation file
    annotate_dir = os.path.dirname(args.annotate)

    # if any import_from_bed keys are present, and they are specified relative to the annotation file,
    # then convert them to absolute paths
    for genomic_element in data["genomic_elements"]:
        if "import_from_bed" in genomic_element:
            if not os.path.isabs(genomic_element["import_from_bed"][0]):
                genomic_element["import_from_bed"][0] = os.path.abspath(os.path.join(
                    annotate_dir, genomic_element["import_from_bed"][0]))

    # draw annotated image
    if not args.use_dna_features_viewer:

        img_width = data["img_width"]
        left, left_grid_line = data["mapping_coordinates_1"]
        right, right_grid_line = data["mapping_coordinates_2"]

        draw_genomic_annotation(img_width, img_width / 6,
                                (left, left_grid_line), (right, right_grid_line),
                                args.op, data, args.no_stroke)
    else:
        if any([k in data for k in ["mapping_coordinates_1", "mapping_coordinates_2"]]):
            print("Warning: mapping_coordinates_1, mapping_coordinates_2 are ignored in "
                  "the DNA Features Viewer mode.")
        draw_with_dna_features_viewer(data, args.op, args.no_stroke)
