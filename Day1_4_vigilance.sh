#!/bin/bash
maxima_input_file="Day1_4_vigilance.txt"
tex_output_file="Day1_4_vigilance_output.tex"

if [ -e $tex_output_file ]
then
  rm $tex_output_file
fi

maxima -b $maxima_input_file
pdflatex $tex_output_file
#Do this twice, so pdflatex can fill in the references
pdflatex $tex_output_file