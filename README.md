Rcode used for analysis and comparison of genetic variability using scaled eigenvectors, for the preprint manuscript "North Pontic crossroads: Mobility in Ukraine from the Bronze Age to the early modern period". This manuscript can be viewed [here](https://www.biorxiv.org/content/10.1101/2024.05.24.595769v1), and any questions regarding any other contents of this study should be directed to the first author, Dr [Lehti Saag](mailto:lehti.saag@ut.ee).
Details regarding the analysis undertaken in this Rscript are provided in the Methods section of the [preprint](https://www.biorxiv.org/content/10.1101/2024.05.24.595769v1.full.pdf) (presently page 29 in the pdf).
Included are the eigenvector and corresponding eigenvalue data for the 100PCs generated with smartPCA in the course of this study ("data.evec") and the metadata table with details for the samples ("ukrad_metadata.csv"), as raw data files used by the Rscript. The Rscript will also read into memory (and then immediately remove) the large metadata table for the [AADR dataset](https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data) to obtain dates and coordinates for reference data used.
