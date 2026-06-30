# variantstring 1.8.7

## Bug fixes

* `get_component_variants()` no longer errors on valid strings that combine a
  phased heterozygous locus with a homozygous locus (e.g. `pfcrt:72_73:C|S_V`).
  The homozygous allele is now correctly carried across all phased components.
* `check_position_string()` now validates every distinct input. Previously, when
  the input contained duplicates, an invalid string positioned after the
  duplicates could be skipped and wrongly pass validation.
* `check_variant_string()` and `check_position_string()` now report the reason
  for failure in the error condition itself (via `conditionMessage()`), instead
  of printing it to the console and raising an empty error. This makes the
  reason available to calling code that catches the error.

## Behaviour changes

* `drop_read_counts()` and `position_from_variant_string()` now preserve the
  format of the input rather than re-normalising it. In particular, position
  ranges (e.g. `72-74`) are kept compact rather than expanded to `72_73_74`, and
  the concise amino-acid notation (e.g. `CV`) is preserved rather than rewritten
  as `C_V`. Genes are still returned in alphabetical order.

## Internal

* Substantial refactor of the internals for speed and to remove duplicated code,
  with no change to the documented behaviour of the affected functions. The
  largest gains are in `drop_read_counts()` and `position_from_variant_string()`
  (now pure string operations), with smaller gains across most other functions.
* Removed the unused `tidyr` dependency.
