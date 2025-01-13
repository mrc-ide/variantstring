
# variantstring

An R package for working with genetic information encoded in *variant string format*.

More specifically, we can define two variations of strings:

- **variant strings** contain position information (gene and codon) along with corresponding amino acids.
- **position strings** contain just the position information, and can be useful in denominator calculation.

This package contains functions to...

- Check for a valid variant string
- Check for a valid position string
- Extract a position string from a variant string
- Compare two variant strings to look for a match (with or without ambiguity)
- Compare a position string against a variant string to look for a match
- Convert a variant string into an ordered version
- Convert a position string into an ordered version
- Extract all single-locus variants from a variant string
