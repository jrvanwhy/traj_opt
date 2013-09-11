#!/bin/bash

# This script builds the documentation for the project. It also updates README in the folder itself (the
# root folder of this repository).

# Abort upon encountering an error
set -e

# Make sure we're in the right directory
cd `dirname $0`

# Search for pdflatex
PDFLATEX=`which pdflatex`
if [ "x$PDFLATEX" == "x" ]
then
	# We did not find pdflatex. Spit out an error to the console and exit.
	echo 1>&2 "Could not find pdflatex! Aborting."
fi
echo "Found pdflatex at: $PDFLATEX"

# Create the PDF of the documentation
echo "Compiling LaTeX documentation into a PDF"
$PDFLATEX documentation.tex

# Generate README.md from the LaTeX sources using pandoc
# If pandoc fails to read the latex, then we should switch it to read from
# the PDF instead
# Search for pandoc
PANDOC=`which pandoc`
if [ "x$PANDOC" == "x" ]
then
	# We did not find pandoc. Inform the user of this fact.
	echo "Could not find pandoc, unable to update README.md"
fi
echo "Found pandoc at: $PANDOC"

# Create the readme
echo "Compiling LaTeX documentation into README.md"
$PANDOC -f latex documentation.tex -t markdown_github -o ../README.md
