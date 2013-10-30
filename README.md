PeakTracker
===========

Peak value comparisons across UCSC .[narrow|broad]peak genome track files

Author: Nima C. Emami, UCSF BMI (Bioinformatics)

File Format Specification: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
						   http://genome.ucsc.edu/FAQ/FAQformat.html#format13

Dependency: Python argparse module

Python version: written for Python 2.7(.5)

License: github.com/nemami/PeakTracker/blob/master/LICENSE

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

Flags: For detailed information about flags and default values, please invoke:
		
			python PeakTracker.py
					OR
			python PeakTracker.py -h

	   on the command line.

Feedback: Feel free to fork & customize! Also, feel free to contact me at the
		  address listed on my github with any feedback.