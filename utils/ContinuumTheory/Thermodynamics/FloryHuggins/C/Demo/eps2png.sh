#!/bin/bash

fname=$1
fname=${fname%.*}

gs -dSAFER -dEPSCrop -r600 -sDEVICE=pngalpha -o $fname.png $fname.eps
