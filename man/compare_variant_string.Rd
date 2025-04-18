% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{compare_variant_string}
\alias{compare_variant_string}
\title{Compares variant strings to look for a match}
\usage{
compare_variant_string(target_string, comparison_strings)
}
\arguments{
\item{target_string}{a single variant string that we want to compare. Cannot
contain any heterozygous calls.}

\item{comparison_strings}{a vector of variant strings against which the
target is compared.}
}
\description{
Compares a target variant string against a vector of comparison strings. A
match is found if every amino acid at every codon position in every gene of
the target is also found within the comparison. Note that ambiguous matches
may occur if there are multiple heterozygous loci in the comparison. In this
case, the target may or may not be within this sample. A match is recorded
but a second output also flags this as an ambiguous match.
}
