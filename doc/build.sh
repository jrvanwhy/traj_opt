#!/bin/bash

# This script builds the documentation for the project. It also updates README in the folder itself (the
# root folder of this repository).

# Abort upon encountering an error
#set -e

# Make sure we're in the right directory
cd `dirname $0`

# Search for pdflatex
PDFLATEX=`which pdflatex`
if [ "x$PDFLATEX" == "x" ]
then
	# We did not find pdflatex. Spit out an error to the console and exit.
	echo 1>&2 "Could not find pdflatex! Documentation build cannot continue without pdflatex."
	echo On Ubuntu-based distributions, installing the texlive-latex-recommended package will bring in
	echo all the necessary dependencies to build the documentation for this project.
	exit
fi
echo "Found pdflatex at: $PDFLATEX"

# Create the PDF of the documentation
echo "Compiling LaTeX documentation into a PDF"
$PDFLATEX documentation.tex
