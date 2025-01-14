

#------------------------------------------------
#' @title Check for a valid variant string
#'
#' @description
#' Checks that an input string (or a vector of strings) matches the required
#' format for a variant string.
#'
#' @param x a character string or vector of character strings.
#'
#' @importFrom stringr str_split
#' @importFrom stats na.omit
#'
#' @export

check_variant_string <- function(x) {

  # CHECKS
  # - is character string
  # - no empty genes
  # - contains exactly three characteristics per gene
  # - gene name (first characteristic):
  #     - contains only allowed characters
  #     - no duplicated gene names
  # - codon position (second characteristic):
  #     - contains only allowed characters
  #     - positions are non-zero
  #     - positions are not duplicated
  #     - positions are sorted increasing
  # - amino acid (third characteristic):
  #     - contains only allowed characters
  #     - cannot start or end with / or | symbols
  #     - cannot contain adjacent / or | symbols
  #     - number of amino acid loci matches number of codons
  #     - a single amino acid group cannot contain both | and / (no mixed phasing)
  #     - all amino acid loci with phased heterozygotes contain the same number of amino acids. This is true both within and between genes.
  #     - hets cannot have the same amino acid multiple times

  # is character vector
  stopifnot(all(is.character(x)))

  # create objects for storing which rows fail. This allows for more
  # informative error messages in which several issues can be listed at once
  # rather than breaking at the first issue
  n <- length(x)
  valid <- rep(TRUE, n)
  reason <- rep(NA, n)

  # get list of valid amino acid characters. Do this once here to avoid
  # repetition for element of x
  IUPAC_df <- allowed_amino_acids()
  valid_amino_characters <- paste0("^[", paste(IUPAC_df$IUPAC_amino_acid_code, collapse = ""), "_/|]+$")

  # input may be a single string or a vector of strings. Loop through elements
  for (i in 1:n) {

    # split into genes
    s1 <- stringr::str_split(x[i], ";", n = Inf, simplify = FALSE)[[1]]
    n_genes <- length(s1)

    if (any(s1 == "")) {
      valid[i] <- FALSE
      reason[i] <-"contains one or more empty genes"
      next()
    }

    # keep track of gene names
    gene_names <- rep(NA, n_genes)

    # keep track of the number of distinct amino acid calls in phased
    # heterozygous loci. This must be the same over all genes
    gene_n_phased <- rep(NA, n_genes)

    # loop over genes
    for (j in seq_along(s1)) {

      # split this gene into characteristics
      z <- stringr::str_split(s1[j], ":", n = Inf, simplify = FALSE)[[1]]
      gene_names[j] <- z[1]

      if (length(z) != 3) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene %s does not contain three characteristics separated by a colon", j)
        next()
      }

      if (!grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", z[1])) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene name %s contains invalid characters", z[1])
        next()
      }

      if (!grepl("^[0-9_]+$", z[2])) {
        valid[i] <- FALSE
        reason[i] <- "codon positions contain invalid characters"
        next()
      }

      # split second characteristic by codon and make numeric
      y <- stringr::str_split(z[2], "_", n = Inf, simplify = FALSE)[[1]] |>
        as.numeric()

      if (any(y == 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon position cannot be zero"
        next()
      }

      if (any(duplicated(y))) {
        valid[i] <- FALSE
        reason[i] <- "codon contains duplicated positions"
        next()
      }

      if (!all(diff(y) > 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon positions must be sorted increasing"
        next()
      }

      if (!grepl(valid_amino_characters, z[3])) {
        valid[i] <- FALSE
        reason[i] <- "amino acid sequence contains invalid characters. See ?allowed_amino_acids()"
        next()
      }

      # remove underscores
      no_underscore <- gsub("_", "", z[3])

      if (grepl("^[/|]|[/|]$", no_underscore)) {
        valid[i] <- FALSE
        reason[i] <- "amino acid sequence cannot start or end with / or | symbols"
        next()
      }

      if (grepl("[/|]{2,}", no_underscore)) {
        valid[i] <- FALSE
        reason[i] <- "amino acid sequence cannot contain adjacent / or | symbols"
        next()
      }

      # split third characteristic using a regular expression to group
      # sequences of letters connected by '/' or '|'
      a <- unlist(strsplit(no_underscore, "(?<=[A-Z])(?=[^/|])", perl = TRUE))

      if (length(a) != length(y)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("number of amino acid loci (%s) must equal the number of codon positions (%s)",
                             length(a), length(y))
        next()
      }

      if (any(grepl("/", a) & grepl("\\|", a))) {
        valid[i] <- FALSE
        reason[i] <- paste("a single amino acid locus cannot contain both the / and",
                           "the | symbols. Heterozygous loci must be either entirely phased or",
                           "entirely unphased", collapse = "")
        next()
      }

      # checks on phased heterozygous loci (if present)
      phased_het <- a[grepl("\\|", a)]
      if (length(phased_het) != 0) {

        n_phased <- mapply(length, strsplit(phased_het, "\\|", ))
        if (length(unique(n_phased)) != 1) {
          valid[i] <- FALSE
          reason[i] <- paste("if there are phased heterozygous loci within a gene then the number of",
                             "distinct amino acids must be the same over all of these loci. Variable",
                             "degrees of phasing are not allowed", collapse = "")
          next()
        }

        # store number of distinct amino acids in phased calls for this
        # gene. This will be tested for consistency across genes
        gene_n_phased[j] <- n_phased[1]
      }

      # split into distinct amino acids
      q <- strsplit(a, "[/|]")

      if (any(mapply(function(x) any(duplicated(x)), q))) {
        valid[i] <- FALSE
        reason[i] <- "the same amino acid cannot be present more than once in a heterozygous call"
        next()
      }

    } # end loop over genes

    if (length(unique(na.omit(gene_n_phased))) > 1) {
      valid[i] <- FALSE
      reason[i] <- paste("if there are phased heterozygous loci then the same number of distinct amino",
                         "acids must be present at all of these loci AND over all genes", collapse = "")
      next()
    }

    if (any(duplicated(gene_names))) {
      valid[i] <- FALSE
      reason[i] <- "duplicated gene names"
      next()
    }

  } # end loop over elements of x

  # if any invalid entries then print all issues before exit
  if (any(!valid)) {
    message("The following issues were found:")
    for (i in seq_along(valid)) {
      if (!valid[i]) {
        message(sprintf("  - entry %s: %s", i, reason[i]))
      }
    }

    stop()
  }

  invisible(TRUE)
}

#------------------------------------------------
# Expand variant strings into their separate elements
#' @noRd
expand_variant_string <- function(x) {

  # checks
  check_variant_string(x)

  # remove underscores from amino acids

  # expand each element into a list
  ret <- mapply(function(s1) {
    mapply(function(s2) {
      list(name = s2[1],
           pos = as.numeric(strsplit(s2[2], "_")[[1]]),
           aa = strsplit(gsub("_", "", s2[3]), "(?<=[A-Z])(?=[^/|])", perl = TRUE)[[1]])
    }, strsplit(s1, ":"), SIMPLIFY = FALSE)
  }, strsplit(x, ";"), SIMPLIFY = FALSE)

  return(ret)
}

#------------------------------------------------
#' @title Check for a valid position string
#'
#' @description
#' Checks that an input string (or a vector of strings) matches the required
#' format for a position string. This is equivalent to a full variant string but
#' with the amino acid information removed, so just giving the gene name(s) and
#' position(s).
#'
#' @param x a character string or vector of character strings.
#'
#' @export

check_position_string <- function(x) {

  # CHECKS
  # - is character string
  # - no empty genes
  # - contains exactly two characteristics per gene
  # - gene name (first characteristic):
  #     - contains only allowed characters
  #     - no duplicated gene names
  # - codon position (second characteristic):
  #     - contains only allowed characters
  #     - positions are non-zero
  #     - positions are not duplicated
  #     - positions are sorted increasing

  # is character vector
  stopifnot(all(is.character(x)))

  # create objects for storing which rows fail. This allows for more
  # informative error messages in which several issues can be listed at once
  # rather than breaking at the first issue
  n <- length(x)
  valid <- rep(TRUE, n)
  reason <- rep(NA, n)

  # input may be a single string or a vector of strings. Loop through elements
  for (i in 1:n) {

    # split into genes
    s1 <- stringr::str_split(x[i], ";", n = Inf, simplify = FALSE)[[1]]
    n_genes <- length(s1)

    if (any(s1 == "")) {
      valid[i] <- FALSE
      reason[i] <-"contains one or more empty genes"
      next()
    }

    # keep track of gene names
    gene_names <- rep(NA, n_genes)

    # loop over genes
    for (j in seq_along(s1)) {

      # split this gene into characteristics
      z <- stringr::str_split(s1[j], ":", n = Inf, simplify = FALSE)[[1]]
      gene_names[j] <- z[1]

      if (length(z) != 2) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene %s does not contain two characteristics separated by a colon", j)
        next()
      }

      if (!grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", z[1])) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene name %s contains invalid characters", z[1])
        next()
      }

      if (!grepl("^[0-9_]+$", z[2])) {
        valid[i] <- FALSE
        reason[i] <- "codon positions contain invalid characters"
        next()
      }

      # split second characteristic by codon and make numeric
      y <- stringr::str_split(z[2], "_", n = Inf, simplify = FALSE)[[1]] |>
        as.numeric()

      if (any(y == 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon position cannot be zero"
        next()
      }

      if (any(duplicated(y))) {
        valid[i] <- FALSE
        reason[i] <- "codon contains duplicated positions"
        next()
      }

      if (!all(diff(y) > 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon positions must be sorted increasing"
        next()
      }

    } # end loop over genes

    if (any(duplicated(gene_names))) {
      valid[i] <- FALSE
      reason[i] <- "duplicated gene names"
      next()
    }

  } # end loop over elements of x

  # if any invalid entries then print all issues before exit
  if (any(!valid)) {
    message("The following issues were found:")
    for (i in seq_along(valid)) {
      if (!valid[i]) {
        message(sprintf("  - entry %s: %s", i, reason[i]))
      }
    }

    stop()
  }

  invisible(TRUE)
}

#------------------------------------------------
# Expand position strings into their separate elements
#' @noRd
expand_position_string <- function(x) {

  # checks
  check_position_string(x)

  # expand each element into a list
  ret <- mapply(function(s1) {
    mapply(function(s2) {
      list(name = s2[1],
           pos = as.numeric(strsplit(s2[2], "_")[[1]]))
    }, strsplit(s1, ":"), SIMPLIFY = FALSE)
  }, strsplit(x, ";"), SIMPLIFY = FALSE)

  return(ret)
}

#------------------------------------------------
#' @title Check that no heterozygous calls within a variant string
#'
#' @description
#' Check that no heterozygous calls within a variant string.
#'
#' @param x a variant string or vector of variant strings.
#'
#' @export

check_nohet_variant_string <- function(x) {

  # checks
  stopifnot(all(is.character(x)))
  check_variant_string(x)

  if (any(grepl("[/|]", x))) {
    stop("input string cannot contain any heterozygous sites")
  }

  invisible(TRUE)
}

#------------------------------------------------
#' @title Extract a position string from a variant string
#'
#' @description
#' Extract a position string from a variant string by stripping the amino acids.
#'
#' @param x a character string or vector of character strings.
#'
#' @export

position_from_variant_string <- function(x) {

  # checks
  stopifnot(all(is.character(x)))
  check_variant_string(x)

  # strip amino acids
  ret <- mapply(function(s1) {
    mapply(function(s2) {
      paste(s2[1:2], collapse = ":")
    }, strsplit(s1, ":")) |>
      paste(collapse = ";")
  }, strsplit(x, ";"))

  return(ret)
}

#------------------------------------------------
#' @title Reorders a variant string in alphabetical order of genes
#'
#' @description
#' Reorders a variant string in alphabetical order of genes. This can be useful
#' when checking for duplicated strings, as it allows us to account for the same
#' genes being present in a different order.
#'
#' @param x a character string or vector of character strings.
#'
#' @export

order_variant_string <- function(x) {

  # checks
  stopifnot(all(is.character(x)))
  check_variant_string(x)

  # sort genes in alphabetical order
  ret <- mapply(function(s1) {
    gene_names <- mapply(function(s2) s2[1], strsplit(s1, ":"))
    s1[order(gene_names)] |>
      paste(collapse = ";")
  }, strsplit(x, ";"))


  return(ret)
}

#------------------------------------------------
#' @title Reorders a position string in alphabetical order of genes
#'
#' @description
#' Reorders a position string in alphabetical order of genes. This can be useful
#' when checking for duplicated strings, as it allows us to account for the same
#' genes being present in a different order.
#'
#' @param x a character string or vector of character strings.
#'
#' @export

order_position_string <- function(x) {

  # checks
  stopifnot(all(is.character(x)))
  check_position_string(x)

  # sort genes in alphabetical order
  ret <- mapply(function(s1) {
    gene_names <- mapply(function(s2) s2[1], strsplit(s1, ":"))
    s1[order(gene_names)] |>
      paste(collapse = ";")
  }, strsplit(x, ";"))


  return(ret)
}

#------------------------------------------------
#' @title Compares variant strings to look for a match
#'
#' @description
#' Compares a target variant string against a vector of comparator strings. A
#' match is found if every amino acid at every codon position in every gene of
#' the target is also found within the comparator. Note that ambiguous matches
#' may occur if there are multiple mixed calls in the comparator, in which case
#' the target may or may not be within this sample. In this case, a match is
#' recorded but a second output also flags this as an ambiguous match.
#'
#' @param target_string a single variant string that we want to compare. Cannot
#'   contain any heterozygous calls.
#' @param comparison_strings a vector of variant strings against which the
#'   target is compared.
#'
#' @export

compare_variant_string <- function(target_string, comparison_strings) {

  # checks
  check_variant_string(target_string)
  stopifnot(length(target_string) == 1)
  check_nohet_variant_string(target_string)
  check_variant_string(comparison_strings)

  # expand both input strings into lists
  target_list <- expand_variant_string(target_string)[[1]]
  comparison_list <- expand_variant_string(comparison_strings)
  n_target_genes <- length(target_list)
  n_comparisons <- length(comparison_list)

  target_names <- mapply(function(x) x$name, target_list)
  match_vec <- rep(FALSE, n_comparisons)
  ambiguous_vec <- rep(FALSE, n_comparisons)
  for (i in 1:n_comparisons) {

    comparison_names <- mapply(function(x) x$name, comparison_list[[i]])
    if (all(target_names %in% comparison_names)) {

      m <- match(target_names, comparison_names)
      match_vec_i <- rep(FALSE, n_target_genes)
      ambiguous_vec_i <- rep(0, n_target_genes)
      for (j in 1:n_target_genes) {

        target_pos <- target_list[[j]]$pos
        comparison_pos <- comparison_list[[i]][[m[j]]]$pos
        if (all(target_pos %in% comparison_pos)) {

          m2 <- match(target_pos, comparison_pos)
          target_aa <- target_list[[j]]$aa
          comparison_aa <- comparison_list[[i]][[m[j]]]$aa[m2]

          match_vec_i[j] <- all(mapply(function(a, b) grepl(a, b), target_aa, comparison_aa))
          ambiguous_vec_i[j] <- sum(grepl("[|/]", comparison_aa))
        }
      }
      match_vec[i] <- all(match_vec_i)
      ambiguous_vec[i] <- (sum(ambiguous_vec_i) > 1)
    }
  }

  return(list(match = match_vec,
              ambiguous = ambiguous_vec))
}

#------------------------------------------------
#' @title Compares a position strings against variant strings to look for a match
#'
#' @description
#' Compares a target position string against a vector of comparator strings. A
#' match is found if every codon position in every gene of the target is also
#' found within the comparator (irrespective of the observed amino acids).
#'
#' @param target_string a single position string that we want to compare.
#' @param comparison_strings a vector of variant strings against which the
#'   target is compared.
#'
#' @export

compare_position_string <- function(target_string, comparison_strings) {

  # checks
  check_position_string(target_string)
  stopifnot(length(target_string) == 1)
  check_variant_string(comparison_strings)

  # expand both input strings into lists
  target_list <- expand_position_string(target_string)[[1]]
  comparison_list <- expand_variant_string(comparison_strings)
  n_target_genes <- length(target_list)
  n_comparisons <- length(comparison_list)

  target_names <- mapply(function(x) x$name, target_list)
  match_vec <- rep(FALSE, n_comparisons)
  for (i in 1:n_comparisons) {

    comparison_names <- mapply(function(x) x$name, comparison_list[[i]])
    if (all(target_names %in% comparison_names)) {

      m <- match(target_names, comparison_names)
      match_vec_i <- rep(FALSE, n_target_genes)
      for (j in 1:n_target_genes) {

        target_pos <- target_list[[j]]$pos
        comparison_pos <- comparison_list[[i]][[m[j]]]$pos
        match_vec_i[j] <- all(target_pos %in% comparison_pos)
      }
      match_vec[i] <- all(match_vec_i)
    }
  }

  return(match = match_vec)
}

#------------------------------------------------
#' @title Extract all single-locus variants from a variant string
#'
#' @description
#' Takes a vector of variant strings, potentially with information at multiple
#' codon positions or genes, and returns variant strings corresponding to all
#' unique single-locus variants within the input. For example, crt:72_73:C_N/V
#' can be extracted to crt:72:C and crt:73:N/V.
#'
#' @param x a vector of variant strings.
#'
#' @export

extract_single_locus_variants <- function(x) {

  # checks
  check_variant_string(x)

  # split into single-locus variants
  ret <- mapply(function(y) {
    mapply(function(z) {
      sprintf("%s:%s:%s", z$name, z$pos, z$aa)
    }, y, SIMPLIFY = FALSE) |>
      unlist()
  }, expand_variant_string(x), SIMPLIFY = FALSE) |>
    unlist() |>
    unique() |>
    sort()

  return(ret)
}

#------------------------------------------------
#' @title List allowed amino acids
#'
#' @description
#' Returns a data.frame of allowed amino acid single-letter codes. These come
#' from IUPAC (International Union of Pure and Applied Chemistry),
#' \href{https://www.bioinformatics.org/sms/iupac.html}{see here} for details.
#'
#' @export

allowed_amino_acids <- function() {

  # get the full path to the file in the extdata directory
  file_path <- system.file("extdata", "IUPAC_codes.rds", package = "variantstring")

  # read the data and return
  ret <- readRDS(file_path)
  return(ret)
}
