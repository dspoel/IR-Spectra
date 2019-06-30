# IR-Spectra
Code and examples to compute IR spectra from normal mode analysis

IR spectra can be computed based on normal mode analysis from either quantum chemistry or force fields. In this repository we gather information and scripts to compute spectra in the context of the GROMACS (http://www.gromacs.org) simulation software.

Force field calculations of frequency-dependent properties have been published here: https://pubs.acs.org/doi/10.1021/acs.jpca.8b09867 and eventually the resulting spectra will be published on http://virtualchemistry.org.

Note for using the scripts:
===========================

GROMACS should be compiled in double precision using:

cmake -DGMX_DOUBLE=ON 

Important literature
====================

https://pubs.acs.org/doi/10.1021/acs.chemrev.9b00007 Review focusing on high-accuracy computational spectroscopy.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4877123/

https://link.springer.com/article/10.1007/s00249-014-1005-6

https://link.springer.com/article/10.1007/s00249-013-0927-8

