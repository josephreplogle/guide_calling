# guide_calling
Tools for sgRNA calling in direct capture Perturb-seq data

This repository contains scripts and notebooks used in "Scalable single-cell CRISPR screens by direct guide RNA capture and targeted library enrichment", Replogle JM et al., Nature Biotechnology 2020. 

At present, these are simply presented as used to process a single lane of 10x data in the manuscript, with notebooks serving as examples.

The inputs to mixed_model_guide_calling.ipynb are:

(1) a table where each row represents a read aligned to a unique guide identity, cellranger-corrected cell barcode, and cellranger-corrected UMI

(2) the valid cell barcodes as determined by cellranger based on the scRNA-seq data

An example can be downloaded: https://www.dropbox.com/sh/vvswv4t8tdophbl/AABPbg8yZC_gnvxDwOF9c2W6a?dl=0. This is the data from UPR experiment with sgRNA constant region 1 with capture sequence 1 in the stem loop.

These datasets are merged and then collapsed to a table of guide UMI counts for each cell.

Then, for each guide, we fit a Poisson-Gaussian mixed model to determine thresholds for sgRNA positive and negative cells.

The output table "cell_identities.csv" contains for a row for every positive sgRNA - cell call. Cells containing multiple guides will be represented in multiple rows. This table is formatted to serve as input to Tom Norman's https://github.com/thomasmaxwellnorman/perturbseq_demo.




