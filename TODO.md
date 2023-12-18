- [ ] Add --eval=reference.gtf
- [ ] From Ross regarding post-processing:

> [9:49 am] Ross Crowhurst
Here is an easy one: BLATSp vs swissprot & Arabidpsis and check query is with set thresholds of reference - if so accept; If not move to BLASTp vs Uniref90 or Refeq (or some other predetermined model species) - same deal accept if within threshold limits. Else BLASTn of cds vs NCBI nt (really scrapping the bottom of the barrel here). If not a hit to anything then chances are its garbage and should be removed. Some ppl might try to claim its a unique protein to the genotype but in 20 years I have never seen one of those be supported - mostly this category is garbage. The screen agains NCBI nt also assists to classify "bits" as well retroposonss etc. Idea being you want to remove garbage predictions - as this does take time you can see why some papers just filter out by size.

- [ ] From Cecilia:

> https://github.com/zhaotao1987/SynNet-Pipeline

- [ ] From Ross:

> https://www.biorxiv.org/content/10.1101/096529v2.full.pdf

> Don't use `-exclude_partial`