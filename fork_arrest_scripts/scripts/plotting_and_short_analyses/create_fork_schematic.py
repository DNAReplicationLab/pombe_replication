import sys
import itertools
import argparse

if __name__ == '__main__':

    desc = """Usage: python create_fork_schematic.py [-obstacle] [-reverse] <tex_file_prefix>
    We will create several files with the filename <tex_file_prefix>.number.tex
    The number will go 0000, 0001, etc.
    These tex files when converted to png will form the frames of the animation of a moving fork
    schematic where nascent DNA is being synthesized with a time-varying BrdU incorporation density.
    The optional flag -obstacle will create a fork schematic with an obstacle at the left fork.
    The optional flag -reverse will create a fork schematic but with gradients reversed.
    To specify this flag, remove the square brackets and type -obstacle before the tex_file_prefix.
    If you do not want this, do not specify the flag i.e. omit the square brackets and everything inside.
    """
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-obstacle", action="store_true", help="Include an obstacle in the schematic.")
    parser.add_argument("-reverse",  action="store_true", help="Reverse the gradients in the schematic.")
    parser.add_argument("tex_file_prefix", type=str, help="The prefix of the TeX file to process.")

    args = parser.parse_args()

    is_obstacle = args.obstacle
    is_reverse = args.reverse
    tex_file_prefix = args.tex_file_prefix

    # create the tex files with the given prefix using the following template
    tex_file_string = r'''\documentclass[tikz,border=3.14mm]{standalone}
\usetikzlibrary{arrows.meta}
\usetikzlibrary{fadings}
\tikzfading[name=fade right, left color=transparent!$opacity_global1, right color=transparent!$opacity_global2]
\tikzfading[name=fade left, right color=transparent!$opacity_global1, left color=transparent!$opacity_global2]
\tikzfading[name=fade left 1, right color=transparent!$opacity_global1, left color=transparent!$opacity_level1]
\tikzfading[name=fade left 2, right color=transparent!$opacity_level2, left color=transparent!$opacity_global2]
\begin{document}

\begin{tikzpicture}

    % Set boundaries of picture
    \clip (-6.5,-5) rectangle (6.5,3);

    % Draw nascent strand of replication fork
    \draw[<-,red,line width=3pt,shift={(0,-0.5)}] (0,2) -- ($start_right,2) -- ($end_right,0.5);
    \draw[->,red,line width=3pt,rotate around x=180,shift={(0,-0.5)}] (0,2) -- ($start_right,2) -- ($end_right,0.5);
    \draw[->,red,line width=3pt,rotate around y=180,shift={(0,-0.5)}] (0,2) -- ($start_left,2) -- ($end_left,0.5);
    \draw[<-,red,line width=3pt,rotate around y=180,rotate around x=180,shift={(0,-0.5)}] (0,2) -- ($start_left,2) -- ($end_left,0.5);
    
    % Draw gradient for nascent strand
    \fill[white,path fading=fade right] (0,2) rectangle (6.5,-2);
    %remove_if_no_obstacle \fill[white,path fading=fade left,rotate around y=180] (0,2) rectangle (6.5,-2);
    %remove_if_obstacle \fill[white,path fading=fade left 1,rotate around y=180] (0,2) rectangle ($obstacle_pos,-2);
    %remove_if_obstacle \fill[white,path fading=fade left 2,rotate around y=180] ($obstacle_pos,2) rectangle (6.5,-2);
    % The above set up of two gradients does not work very well mathematically speaking
    % (because the second gradient in the presence of an obstacle is stretched out more than the first gradient)
    % but it looks ok to the eye and we are fine with this as this is just a schematic.
    
    %$obstacle_line\draw[line width=3pt,blue,rotate around y=180] ($obstacle_pos,2) -- ($obstacle_pos,-2);

    % Draw parent strand of replication fork
    \draw[line width=3pt] (0,2) -- ($start_right,2) -- ($end_right,0.5) -- (6.5,0.5);
    \draw[line width=3pt,rotate around x=180] (0,2) -- ($start_right,2) -- ($end_right,0.5) -- (6.5,0.5);
    \draw[line width=3pt,rotate around y=180] (0,2) -- ($start_left,2) -- ($end_left,0.5) -- (6.5,0.5);
    \draw[line width=3pt,rotate around y=180,rotate around x=180] (0,2) -- ($start_left,2) -- ($end_left,0.5) -- (6.5,0.5);

    % Draw clock
    \draw[ultra thick] (-1,-3.45) circle [radius=0.7cm];
    \fill[gray, ultra thick] (-1,-3.45) circle [radius=0.04cm];
    \draw[gray, ultra thick,rotate around={-$angle:(-1,-3.45)}] (-1,-3.45) -- (-1,-2.95); %changing

    % Draw square
    \fill[red,opacity=$opacity] (0.9,-3) rectangle (2.05,-4.15);

    % Annotate replication fork, clock, and square
    \node at (0,2.55) {\huge Parental};
    \node[red] at (0,-0.85) {\huge Nascent};
    \node[black] at (-1,-4.65) {\huge Time};
    \node[black] at (1.45,-4.65) {\huge BrdU};

\end{tikzpicture}
\end{document}
'''

    # set overall parameters for the movie
    step_size = 0.03
    duration = 250
    gradient_duration = 100

    # set obstacle parameters
    obstacle_frame = 20
    obstacle_duration = 50 if is_obstacle else 0
    obstacle_pos = (1 + 1.5 + obstacle_frame * step_size)

    # set opacity levels for gradients
    initial_opacity = 73
    final_opacity = 10
    if is_reverse:
        initial_opacity, final_opacity = final_opacity, initial_opacity # swap the gradients if needed
    max_opacity = 100
    min_opacity = 10
    opacity_step = (initial_opacity - final_opacity) / gradient_duration
    opacity_level1 = initial_opacity - obstacle_frame * opacity_step
    opacity_level2 = initial_opacity - (obstacle_frame + obstacle_duration) * opacity_step

    # set iterators for left and right fork positions etc.
    index_right_iter = itertools.chain.from_iterable([range(0, duration), itertools.repeat(duration, 20)])
    if not is_obstacle:
        index_left_iter, index_right_iter = itertools.tee(index_right_iter)
    else:
        index_left_iter = itertools.chain.from_iterable([range(0, obstacle_frame),
                                                         itertools.repeat(obstacle_frame, obstacle_duration),
                                                         range(obstacle_frame, duration - obstacle_duration),
                                                         itertools.repeat(duration - obstacle_duration, 20)])

    for cnt, index_left, index_right in zip(itertools.count(), index_left_iter, index_right_iter):

        # set filename and open file
        tex_file_name = tex_file_prefix + "." + str(cnt).zfill(4) + ".tex"

        with open(tex_file_name, "w") as f:

            # write the tex file with the given parameters
            start_pos_left = 1 + index_left * step_size
            end_pos_left = start_pos_left + 1.5

            start_pos_right = 1 + index_right * step_size
            end_pos_right = start_pos_right + 1.5

            angle = index_right / duration * 300
            opacity = max((initial_opacity - index_right * opacity_step), final_opacity)/100

            f.write(tex_file_string.replace("$start_left", str(start_pos_left)).
                    replace("$end_left", str(end_pos_left)).
                    replace("$start_right", str(start_pos_right)).
                    replace("$end_right", str(end_pos_right)).
                    replace("$angle", str(angle)).
                    replace("$opacity_level1", str(opacity_level1)).
                    replace("$opacity_level2", str(opacity_level2)).
                    replace("$opacity_global1", str(min_opacity) if is_reverse else str(max_opacity)).
                    replace("$opacity_global2", str(max_opacity) if is_reverse else str(min_opacity)).
                    replace("$opacity", str(opacity)).
                    replace("$obstacle_pos", str(obstacle_pos)).
                    replace("%$obstacle_line",
                            "" if (is_obstacle and index_right < obstacle_frame + obstacle_duration) else "%").
                    replace("%remove_if_obstacle ", "" if is_obstacle else "%do_not_draw_this%").
                    replace("%remove_if_no_obstacle ", "" if not is_obstacle else "%do_not_draw_this%"))
