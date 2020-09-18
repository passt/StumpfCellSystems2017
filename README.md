
## Stem Cell Differentiation as a Non-Markov Stochastic Process

Source code related to Stumpf et al. Cell Systems (2017) https://doi.org/10.1016/j.cels.2017.08.009


### Data availability

Single-cell data: http://dx.doi.org/10.17632/g2md5gbhz7.1

Microarray data: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5861/


### Network inference & visualisation

1) Apply k-means clustering with k=3 (see 'NeuronalDifferentiation.Rmd')
2) Run network inference on the complete data (ESC/EPI/NPC) and split data into (ESC/EPI) and (EPI/NPC) timepoints using cell-state assignment from k-means (data from both biological replicates).
3) Prune bottom 95% of edges (maintaining the top 5% PID scores).
4) Run community detection on pruned (ESC/EPI/NPC) network (method: PMID: 20615936).
5) Visualise network w. edge color mapped from ESC/EPI (white-blue; low-high) or EPI/NPC (white-red; low-high) PID scores (or average if present in both; white-grey; low-high). Results are shown in Figure S3b.
6) Prune network, removing modules with no significant changes in (intra-module) PID scores (early vs late). Results are shown in Figure 2E.
7) Conduct network analysis on (6). Results are shown in Figure 2F.
