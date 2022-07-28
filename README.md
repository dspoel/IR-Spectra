# IR-Spectra
Code and examples to compute IR spectra from normal mode analysis

IR spectra can be computed based on normal mode analysis from either quantum chemistry or force fields. In this repository we gather information and scripts to compute spectra in the context of the GROMACS (http://www.gromacs.org) simulation software.

Force field calculations of frequency-dependent properties have been published here: https://pubs.acs.org/doi/10.1021/acs.jpca.8b09867 and eventually the resulting spectra will be published on http://virtualchemistry.org.

We have published two papers related to this work:
+ Theoretical Infrared Spectra: Quantitative Similarity Measures and Force Fields
Henning Henschel, Alfred T. Andersson, Willem Jespers, Mohammad Mehdi Ghahremanpour, and David van der Spoel, https://doi.org/10.1021/acs.jctc.0c00126
+ An Intuitively Understandable Quality Measure for Theoretical Vibrational Spectra
Henning Henschel and David van der Spoel, https://doi.org/10.1021/acs.jpclett.0c01655

Note for using the scripts:
===========================

GROMACS should be compiled in double precision using:

cmake -DGMX_DOUBLE=ON 

Important literature
====================

https://pubs.acs.org/doi/10.1021/acs.chemrev.9b00007 Review focusing on high-accuracy computational spectroscopy.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4877123/ About the NIST Quantitative Infrared Database.

https://link.springer.com/article/10.1007/s00249-014-1005-6 Calculation of the infrared spectra of proteins.

https://link.springer.com/article/10.1007/s00249-013-0927-8 Fast calculation of the infrared spectra of large biomolecules.

Other related projects
======================

The Virtual Multifrequency Spectrometer project http://dreamslab.sns.it/vms/index.php
