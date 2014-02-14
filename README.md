PlosCompBio_2013
================

This repository contains code for generating Figures 2-5 and Table 1 from the following paper.

Perkins TA, Scott TW, Le Menach A, Smith DL (2013) **Heterogeneity, Mixing, and the Spatial Scales of Mosquito-Borne Pathogen Transmission**. *PLoS Computational Biology* 9(12): e1003327. doi:[10.1371/journal.pcbi.1003327](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003327)

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/).

================

### Figure 2

Within the folder ./fig2/, execute `source('figKernels.R')` in **R** to produce ./fig2/output/fig2.pdf.

### Figures 3 & 4

Within the folder ./fig34/, first execute `source('runPhylogram.R')` in **R** and then `source('figMix.R')` and `source('figHet.R')` to produce ./fig34/output/mixing.pdf and ./fig34/output/heterogeneity.pdf, respectively.

### Figure 5

Within the folder ./fig5/, first execute `source('runInfectionsOverTime.R')` in **R** and then `source('figInfectionsOverTime.R')` to produce ./fig5/output/fig5.pdf.

### Table 1

Within the folder ./tab1/, execute `source('makeTable.R')` in **R** to produce a matrix called `tab` that resides within ./tab1/output/table.RData.
