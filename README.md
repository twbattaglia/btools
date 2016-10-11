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

```R
# Install necessary packages
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
biocLite("ggplot2")
biocLite("vegan")

# Load necessary packages
library(phyloseq)
library(ggplot2)

# Import OTU table + tree + map
phylo <- create_phylo(biom_fp = "otu_table.biom",
                      mappingfile_fp = "mapping_file.txt",
                      tree_fp = "rep_set.tre")
```


## List of Functions
#### Alpha diversity significance
```R
compare_alpha_diversity(phylo, x = "Time", 
                        group = "Treatment", 
                        diversity = "Observed",
                        test_type = "nonparametric", 
                        write = T, 
                        filename = "adiv_results") 
```

#### Beta diversity significance
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

#### PICRUSt metagenomic contributions table + grpah
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

#### Pairwise distances table
```R
jaccard <- diversity_comparison(phylo, distance = "jaccard")
jsd <- diversity_comparison(phylo, distance = "jsd")
unweighted <- diversity_comparison(phylo, distance = "unifrac")
weighted <- diversity_comparison(phylo, distance = "wunifrac")
```

#### Import NanoString data with corrected sample names
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

#### Plot 3D PCA with plotly
Thanks to DESeq2
```R
plotPCA3D(deseq2, intgroup = "Treatment")
```
