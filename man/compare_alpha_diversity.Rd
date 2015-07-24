\name{compare_alpha_diversity}
\alias{alpha_div}
\title{Perform stastistical testing on alpha diversity measures in phyloseq}
\usage{
compare_alpha_diversity(phyloseq,
x = "Day",
group = "Treatment",
diversity = "Observed",
test_type = "nonparametric",
num_perm = 999,
multiple_corrections = F,
write = F,
filename = "results",
... )
}
\description{
Outputs a table of p-values, means and standard deviations found similarliy in QIIME's 'compare_alpha_diversity.py' script.
}
\examples{
compare_alpha_diversity()
}
