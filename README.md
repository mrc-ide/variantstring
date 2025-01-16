
# variantstring

An R package for working with genetic information encoded in *variant string format*.

More specifically, we can define two types of strings that we interested in:

- **variant strings** contain position information (gene and codon) along with corresponding amino acids. Optionally they may also include read counts for each corresponding amino acid.
- **position strings** are a stripped down version of a variant string, containing just the gene and codon position information.

In brief, this package contains functions that...

- Check for correctly formatted variant and position strings.
- Extract a position string from a variant string.
- Compare two variant strings to look for a match (useful in numerator of prevalence calculation). Reports if this is an exact match or an ambiguous match.
- Compare a position string against a variant string to look for a match (useful in denominator of prevalence calculation).
- Convert between string format and a long-form data.frame format.

There are also a few more utility functions not listed here - see the package help for a complete list of functions.

## Release history

The current version is 1.5.1, released 16 Jan 2025.
