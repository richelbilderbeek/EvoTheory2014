#!/bin/bash
maxima_input_file="Day23.txt"
tex_output_file="Day23.tex"

if [ -e $tex_output_file ]
then
  rm $tex_output_file
fi

maxima -b $maxima_input_file
pdflatex $tex_output_file
pdflatex $tex_output_file

