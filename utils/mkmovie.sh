#!/bin/bash

# copy this script in the img directory and run there

ffmpeg -framerate 15 -i fig%03d.png  -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
