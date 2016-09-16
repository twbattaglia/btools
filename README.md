# btools
A set of R functions that help faciliate a lot of tedious processing  
[![Travis-CI Build Status](https://travis-ci.org/twbattaglia/btools.svg?branch=master)](https://travis-ci.org/twbattaglia/btools)

## Install
```bash
install.packages('devtools')
devtools::install_github('twbattaglia/btools')
library(btools)
```

## Import QIIME to phyloseq
The input OTU table in (.biom) format, must be JSON formatted before import. This issue maybe resolved in future updates.

```R
# Install necessary packages
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
biocLite("ggplot2")
biocLite("vegan")
biocLite("randomforests")

# Load necessary packages
library(phyloseq)
library(ggplot2)

# Import OTU table and tree
# Either using function or individual steps
phylo <- create_phylo(biom_fp = "otu_table_json.biom",
                      mappingfile_fp = "mapping_file.txt",
                      tree_fp = "rep_set.tre")


# Or
otutable <- import_biom(BIOMfilename = 'otu_table_json.biom', 
                        treefilename = 'rep_set.tre', 
                        parseFunction = parse_taxonomy_greengenes)

# Import mapping file
mapping <- import_qiime_sample_data(mapfilename = 'mapping_file.txt')

# Merge map and otu table into once phyloseq object
phylo <- merge_phyloseq(otutable, mapping)
```


## List of Functions
#### Alpha Diversity
```R
compare_alpha_diversity(phylo, x = "Time", 
                        group = "Treatment", 
                        diversity = "Observed",
                        test_type = "nonparametric", 
                        write = T, 
                        filename = "adiv_results") 
```

#### Beta Diversity
```R
compare_beta_diversity(phylo, 
                       x = "Time",
                       group = "Treatment",
                       bdiv = "unweighted",
                       test = "adonis", 
                       write = T, 
                       fdr = TRUE,
                       filename = "bdiv_results")
```

#### PICRUSt
```R
contributions <- analyze_contributions(contributions_fp = "metagenomic_contributions.tab", 
                                       mappingfile_fp = "mapping_file.txt")
```

#### Pairwise distances
```R
jaccard <- diversity_comparison(phylo, distance = "jaccard")
jsd <- diversity_comparison(phylo, distance = "jsd")
unweighted <- diversity_comparison(phylo, distance = "unifrac")
weighted <- diversity_comparison(phylo, distance = "wunifrac")

```


