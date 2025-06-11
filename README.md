# Genomic evolution of Salmonella Dublin in cattle and humans in the United States

## Publication 

## Repository Structure and Analysis Pipeline 

### Directory contents

*   `bash` contains all code used to perform the analyis and generate output tables used in RStudio
  
      * subdir `condaenv_yaml` contains yaml files for all conda environments needed for the analysis

*    `R` contains all code used to perform downstream analysis and data visualizations used in the manuscript

      * `abricate` : abricate results 
      * `amrfinder` : amrfinder results 
      * `assembqc` : assembly qc tables 
      * `functionalanno` : 
      * `meta` : metadata for final dataset
      * `pangenome` : 
      * `phylo` : 
      * `readqc` : read qc tables
          * `combinedkreports.xlsx` is too large for Github and can be found at XXX
      * `script` : all R scripts

### Analysis Pipeline 

#### **HPC Component**

Bash scripts should be run in roughly this order: 

1. accessionqc.sh - download reads and generate qc reports
2. trimqc.sh - qc reads and generate qc reports
3. classifyreads.sh - taxonomic classification to check for contamination
4. filterreads.sh - filter for *Enterobacteriaceae* reads
5. assembly.sh - genome assembly
6. assembqc.sh - assembly qc

Once assemblies have been constructed all of the following can be run somewhat concurrently:

* sistr.sh - *in silico* serotyping
* seqsero2.sh - *in silico* serotyping
* mlst.sh - cgMLST assignment
* prokka.sh - genome annotation
* abricate.sh - genome annotation (virulence factors and plasmids only)

The following can be run after Prokka completes
* combineprokkasummary.sh - summarize Prokka results
* amrfinder.sh - AMR annotation
* roary.sh - pangenome annotation and core genome alignment

The following can be run after Roary completes
1. modeltest.sh - determine best nucleotide substitution model for phylogenetic tree building
2. iqtree.sh - build phylogenetic tree

#### **RStudio Component**

Recommended order:

1. meta.R
2. readqc.R
3. assembqc.R
4. abricate.R
5. amrfinder.R
   





