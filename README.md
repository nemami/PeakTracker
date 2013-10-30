PeakTracker
================================================================================

Peak value comparisons across UCSC .[narrow|broad]peak genome track files<br>
Author: Nima C. Emami, UCSF BMI (Bioinformatics)<br>
License: github.com/nemami/PeakTracker/blob/master/LICENSE

================================================================================

Dependencies
============

Python 2.7      (http://www.python.org/getit/releases/2.7/)<br>
Python modules: argparse, numpy

File Format Specifications
==========================

narrowpeak:     http://genome.ucsc.edu/FAQ/FAQformat.html#format12<br>
broadpeak:      http://genome.ucsc.edu/FAQ/FAQformat.html#format13

Usage
=====

Usage case: You have multiple .narrowpeak and/or .broadpeak files from the
            UCSC genome browser. Without adding the tracks to the browser
            GUI and scanning manually through the genome, you want to know
            which loci have common peaks for your tracks of interest, above
            some signalValue stringency threshold. This is precisely why I
            wrote this script -- it takes multiple .narrowpeak / .broadpeak
            files, maps them to set "windows" in their common reference genome,
            identifies how many hits have mapped to a particular window, while
            satisfying other user stringency criteria, and reports windows
            with associated mean / standard deviation data for signal Values,
            to an output file.

Disclaimer
==========

This package operates on narrow/broadpeak files, which can be relatively large.
If you plan on comparing many of these files using this program, this script 
might crash your laptop -- instead, I recommend running it on a desktop/server 
with more RAM.

Flags
=====

Please refer to:

            python PeakTracker.py
                    OR
            python PeakTracker.py -h

Feedback
========

Feel free to fork & customize! Also, please contact me at the address listed 
on my github with any feedback.