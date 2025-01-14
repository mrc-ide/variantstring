
# variantstring

An R package for working with genetic information encoded in *variant string format*.

More specifically, we can define two types of strings that we interested in:

- **variant strings** contain position information (gene and codon) along with corresponding amino acids.
- **position strings** are a stripped down version of a variant string, containing just the position information but no amino acids.

Variant strings are useful when obtaining the numerator in a prevalence calculation. The corresponding position string is useful when obtaining the denominator.

This package contains functions that...

- Check for a valid variant and position strings.
- Extract a position string from a variant string.
- Compare two variant strings to look for a match (useful in numerator calculation). Reports whether the match is ambiguous or not, as in the case of mixed infections.
- Compare a position string against a variant string to look for a match (useful in denominator calculation).
- Convert a variant and position strings into an ordered versions.
- Extract all single-locus variants from a variant string.

## Release history

The current version is 1.0.0, released 14 Jan 2025.
