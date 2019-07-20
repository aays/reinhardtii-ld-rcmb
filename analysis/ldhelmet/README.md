
## `ldhelmet`

Scripts to calculate the recombination landscape in _C. reinhardtii_ with
LDhelmet 1.9 (Chan et al. 2012). 

### `run_ldhelmet.sh`

Given a block penalty, will run LDhelmet over all fastas in hardcoded dir.

### `find_hotspots.py`

Adapted from Singhal et al. (2015, Science). Summarises LDhelmet
output in windows of specified size, while also calculating recombination
rates in regions flanking each window (with modifiable flank size)

### `landscape_analysis.Rmd`

R Markdown file containing analysis of overall recombination landscape.

### `hotspot_analysis.Rmd`

R Markdown file containing analysis of hotspots in _C. reinhardtii_.

### `frequency_analysis.Rmd`

R Markdown file containing estimation of the frequency of sex in _C. reinhardtii_. 
