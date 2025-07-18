import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import os
import sys


# Script written with Google Gemini

def convert_png_red_to_viridis(input_png_path: str, output_png_path: str):
    """
    Reads a PNG image, converts red shades to the viridis colormap based on
    red intensity, and saves the result as a new PNG file.

    Args:
        input_png_path (str): The path to the input PNG file.
        output_png_path (str): The path where the output PNG will be saved.
    """
    # --- Viridis Colormap Preparation ---
    # Get the 'viridis' colormap from matplotlib.
    # We create a lookup table of 256 colors from the colormap.
    viridis_map = plt.get_cmap('viridis', 256)
    viridis_lookup = (viridis_map(np.linspace(0, 1, 256)) * 255).astype(np.uint8)

    try:
        # Open the input PNG file
        with Image.open(input_png_path) as img:
            print(f"Processing image '{input_png_path}'...")

            # Convert the image to RGBA mode to ensure we have an alpha channel
            rgba_image = img.convert('RGBA')

            # Convert the PIL Image into a NumPy array for fast processing.
            data = np.array(rgba_image)
            output_data = data.copy()

            # Create separate arrays for Red, Green, Blue, and Alpha channels
            r, g, b, a = data[:, :, 0], data[:, :, 1], data[:, :, 2], data[:, :, 3]

            # Define the condition for a pixel to be considered a "red shade".
            # We also check that the pixel is not fully transparent.
            is_red_mask = (r > g) & (r > b) & (r > 50) & (a > 0)

            # --- Map Red Intensity to Viridis (Corrected Logic) ---
            # The original logic failed on gradients where the R channel is constant (e.g., 255).
            # To fix this, we calculate a "redness" metric that measures how much
            # the red channel dominates the green and blue channels.
            # Metric: Redness = R - max(G, B)
            # We use floating point numbers for the calculation to avoid data type issues.
            r_float = r.astype(np.float32)
            g_float = g.astype(np.float32)
            b_float = b.astype(np.float32)

            # Convert the blue obstacle into red
            is_blue_mask = (b > g) & (b > r) & (b > 250) & (a > 0)
            output_data[is_blue_mask] = [255, 0, 0, 255]

            # Calculate the redness metric for all pixels in the image.
            redness = r_float - np.maximum(g_float, b_float)

            # Select the "redness" values only for the pixels that match our mask.
            redness_values_in_mask = redness[is_red_mask]

            if redness_values_in_mask.size > 0:
                # Normalize the collected redness values to span the full colormap (0-1).
                min_val = np.min(redness_values_in_mask)
                max_val = np.max(redness_values_in_mask)

                # Avoid division by zero if all the red shades have the same intensity.
                if max_val > min_val:
                    normalized_values = (redness_values_in_mask - min_val) / (max_val - min_val)
                else:
                    # If all values are the same, map them to the middle of the colormap.
                    normalized_values = np.full(redness_values_in_mask.shape, 0.5)

                # Scale the normalized values to get indices (0-255) for our lookup table.
                # We set this to 225 as we do not want the highest possible colours on the
                # viridis scale as we don't see 100% incorporation in our data.
                indices = (normalized_values * 255).astype(np.uint8)
                indices = np.clip(indices, 0, 225)  # Ensure indices are within bounds

                # Get the new colors from the viridis lookup table.
                new_colors = viridis_lookup[indices]

                # Preserve the original alpha values of the modified pixels.
                new_colors[:, 3] = a[is_red_mask]

                # Apply the new viridis colors to the output data array.
                output_data[is_red_mask] = new_colors

            # Convert the processed NumPy array back into a PIL Image.
            processed_image = Image.fromarray(output_data, 'RGBA')

            # --- Save the New PNG ---
            print(f"Saving new image to '{output_png_path}'...")
            processed_image.save(output_png_path, 'PNG')
            print("Done!")

    except FileNotFoundError:
        print(f"Error: The file '{input_png_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def create_dummy_png(file_path, size=(200, 100)):
    """Creates a simple test PNG with a red gradient."""
    print(f"Creating a dummy PNG for testing at '{file_path}'...")
    # Create a gradient from light pink to pure red
    img = Image.new('RGBA', size, (255, 255, 255, 0))  # Transparent background
    pixels = img.load()

    for x in range(size[0]):
        # The G and B channels decrease as x increases, making the red purer.
        green_blue_val = 220 - int((x / size[0]) * 220)
        for y in range(size[1]):
            pixels[x, y] = (255, green_blue_val, green_blue_val, 255)

    img.save(file_path, 'PNG')
    print("Dummy PNG created.")


if __name__ == '__main__':
    # --- Configuration ---
    input_filename = sys.argv[1] if len(sys.argv) > 1 else ""
    output_filename = sys.argv[2] if len(sys.argv) > 2 else ""

    # if either input or output filename is not provided, error out
    if not input_filename or not output_filename:
        print("Usage: python recolour_fork_schematic_to_viridis.py <input_png> <output_png>")
        sys.exit(1)

    # --- Run the Conversion ---
    convert_png_red_to_viridis(input_filename, output_filename)
