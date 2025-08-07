# Genomic evolution of *Salmonella* Dublin in cattle and humans in the United States

## Publication 

## Repository Structure and Analysis Pipeline 

### Directory contents

*   `bash` contains all code used to perform the analysis and generate output tables used in RStudio
  
      * subdir `condaenv_yaml` contains yaml files for all conda environments needed for the analysis

*    `R` contains all code used to perform downstream analysis and data visualizations used in the manuscript

      * `abricate` : abricate results 
      * `amrfinder` : amrfinder results 
      * `assembqc` : assembly qc tables 
      * `functionalanno` : annotation files, relevant gene lists, and reference files for gene annotation
      * `meta` : metadata for final dataset
      * `pangenome` :  raw data files and subset matrices from permutational analysis
      * `phylo` : SNP distance matrix, and tree file
      * `readqc` : read qc tables
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

1. meta.R - plot Figure 1
2. readqc.R - aggregate read qc details in supplementary tables
     * `combinedkreports.xlsx` can be downloaded [here](https://zenodo.org/records/15652582)
4. assembqc.R - aggregate assembly qc details in supplementary tables
5. abricate.R - create Tables 1 and 4 (including associated statistical tests)
6. amrfinder.R - create Tables 2 and 3; plot Figure 2 (including associated statistical tests)
7. pangenome.R - create Table 5; plot Figure 3A-C (including associated statistical tests)
     * full dataset `gene_presence_absence` tables can be downloaded [here](https://zenodo.org/records/15652582)
9. anno.R - plot Figure 3D
10. snp.R - plot Figure 4B-D (Figure 4A generated with [iTOL](https://itol.embl.de/))
     * As the way the chord diagrams are plotted imply transmission directionality that cannot be reliably inferred, plots were edited in [Flourish](https://flourish.studio/) to resolve this and generate the figures seen in the manuscript.    





