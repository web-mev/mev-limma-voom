# mev-limma-voom

This repository contains a WebMeV-compatible tool for performing differential expression analysis using Bioconductor's Limma+voom package. It uses the `voom` function to estimate the mean-variance relationship prior to use of limma, which was originally developed for microarray platforms. More details at https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29

Note that the default usage assumes you are performing a simple contrast, without consideration for more complex designs. 

The outputs are:
- A tab-delimited file of the differential expression results merged with the normalized counts
- A tab-delimited file of just the normalized counts.

The concatenation of the differential expression results and counts is for convenience with the WebMeV frontend interface since it avoids pulling data from two different files.

---

### To run external of WebMeV:

Either:
- build the Docker image using the contents of the `docker/` folder (e.g. `docker build -t myuser/limma:v1 .`) 
- pull the docker image from the GitHub container repository (see https://github.com/web-mev/mev-limma-voom/pkgs/container/mev-limma-voom)

To run, change into the directory containing your raw/integer count matrix you wish to use. Then:
```
Rscript /usr/local/bin/limma.R \
    /work/<raw/integer counts> \
    <base/control condition samples as CSV-string> \
    <experimental condition samples as CSV-string> \
    <base/control condition name> \
    <experimental condition name> \
    <minimum number of reads cutoff>
```
Note that we mounted your current directory to `/work` inside the Docker container. Hence, your file is relative to `/work`.

The call to the script assumes the following:
- The input file of expression counts is tab-delimited format and contains only integer entries
- The samples in either the control or experimental groups are given as comma-delimited strings and correspond to the column names contained in the raw count matrix. As an example, if the base/control samples are A, B, and C, specify: `"A,B,C"`. The wrapping quotations are not necessary unless the sample names contain whitespace.
- The condition names are regular strings used to help with naming the output file.
- Counts lower than the minimum number of reads (typically 1) will be filtered out.
- The output files will be written to the same directory where the input/raw counts file is located.

The choice of samples in each group should be a proper subset of the samples represented in the samples in the count matrix; you are *not* required to subset the count matrix to include only those samples involved in the contrast of interest.