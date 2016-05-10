# btools
A set of R functions that help faciliate a lot of tedious processing

## Install
```bash
install.packages('devtools')
devtools::install.github('twbttaglia/btools')
library(btools)
```

# Import QIIME to phyloseq
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
otutable <- import_biom(BIOMfilename = 'otu_table_json.biom', 
                        treefilename = 'rep_set.tre', 
                        parseFunction = parse_taxonomy_greengenes)

# Import mapping file
mapping <- import_qiime_sample_data(mapfilename = 'mapping_file.txt')

# Merge map and otu table into once phyloseq object
phylo <- merge_phyloseq(otutable, mapping)
```

## Usage
### Beta diversity
```R
compare_beta_diversity(phylo, x = "Time", group = "Treatment", bdiv = "unweighted", test = "adonis")
```


## List of Functions

#### Alpha Diversity
compare_alpha_diversity()  
plot_nice_boxplot()

#### Beta Diversity
compare_beta_diversity()

#### PICRUSt
analyze_contributions()
plot_contributions()
