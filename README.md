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


# Or the long wayyy
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

# Plot contributions
contributions %>% 
  group_by(Gene, Treatment) %>%
  mutate(Contribution_perc = ContributionPercentOfAllSamples * 100) %>%
  filter(Contribution_perc >= 0) %>%
  select(Gene, family, Contribution_perc) %>%
  mutate(Contribution = Contribution_perc/sum(Contribution_perc) * 100) %>%
  ggplot(arrange(pp_pathway_norm, family), aes(x = Treatment, y = Contribution, fill = family)) + 
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(size = 6)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  theme_light(base_size = 18) + 
  scale_fill_brewer(palette = "Set1")
```

#### Pairwise distances
```R
jaccard <- diversity_comparison(phylo, distance = "jaccard")
jsd <- diversity_comparison(phylo, distance = "jsd")
unweighted <- diversity_comparison(phylo, distance = "unifrac")
weighted <- diversity_comparison(phylo, distance = "wunifrac")
```

#### Import NanoString data
Thanks to NanoStringNorm
```R
genes <- import_rcc("cel_files/")
```

#### BF ratio
```R
phyloseq <- bf_ratio(phyloseq)

# View log2 BF ratio's
phyloseq$log2_bf_ratio
```

#### Plot 3D PCA 
Thanks to DESeq2
```R
plotPCA3D(deseq2, intgroup = "Treatment")
```
