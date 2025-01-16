
#------------------------------------------------
#' @title Check for a valid variant string
#'
#' @description
#' Checks that an input string (or a vector of strings) matches the required
#' format for a variant string.
#'
#' @param x a character string or vector of character strings.
#'
#' @importFrom stringr str_split str_c
#' @importFrom stats na.omit
#'
#' @export

check_variant_string <- function(x) {

  # CHECKS
  # - is character string
  # - no empty genes
  # - contains exactly three elements per gene
  # - gene name (first element):
  #     - contains only allowed characters
  #     - no duplicated gene names
  # - codon position (second element):
  #     - contains only allowed characters
  #     - positions are non-zero
  #     - positions are not duplicated
  #     - positions are sorted increasing
  # - amino acid (third element):
  #     - contains only allowed characters
  #     - cannot start or end with / or | symbols
  #     - cannot contain adjacent / or | symbols
  #     - number of amino acid loci matches number of codon positions
  #     - a single amino acid group cannot contain both | and / (no mixed phasing)
  #     - hets cannot have the same amino acid multiple times
  #     - all amino acid loci with phased heterozygotes contain the same number of amino acids. This is true both within and between genes.
  # - read count (optional fourth element):
  #     - contains only allowed characters
  #     - cannot start or end with / or | symbols
  #     - cannot contain adjacent / or | symbols
  #     - number of read count loci matches number of codon positions
  #     - a single read count group cannot contain both | and / (no mixed phasing)
  #     - the pattern of phased and unphased read counts must match that in amino acids
  #     - read counts must be greater than 0
  #     - number of read counts at each locus must match number of amino acids at each locus
  #     - if read counts are present they must be present in all genes, and vice versa

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

    # keep track of gene names and the number of distinct amino acids in phased
    # heterozygous loci
    gene_name_vec <- rep(NA, n_genes)
    n_phased_vec <- rep(NA, n_genes)
    reads_present_vec <- rep(NA, n_genes)

    # loop over genes
    for (j in seq_along(s1)) {

      # split this gene into elements
      s2 <- stringr::str_split(s1[j], ":", n = Inf, simplify = FALSE)[[1]]
      n_elements <- length(s2)

      if (!(n_elements %in% c(3,4))) {
        valid[i] <- FALSE
        reason[i] <- str_c("genes must contain either three or four elements separated ",
                           "by a colon (name, position, amino acid, optionally read counts)")
        next()
      }

      gene_name <- s2[1]
      codon_pos_string <- s2[2]
      codon_aa_string <- s2[3]
      codon_reads_string <- ifelse(n_elements == 4, s2[4], NA)

      # ------------------------------------------------------------
      # check gene name

      if (!grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", gene_name)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene name %s contains invalid characters", gene_name)
        next()
      }

      # store gene name. This will be tested for duplication across genes
      gene_name_vec[j] <- gene_name

      # ------------------------------------------------------------
      # check codon positions

      if (!grepl("^[0-9_]+$", codon_pos_string)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("codon position %s contains invalid characters", codon_pos_string)
        next()
      }

      # split codon position into separate loci and make numeric
      codon_pos <- stringr::str_split(codon_pos_string, "_", n = Inf, simplify = FALSE)[[1]] |>
        as.numeric()

      if (any(codon_pos == 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon position cannot be zero"
        next()
      }

      if (any(duplicated(codon_pos))) {
        valid[i] <- FALSE
        reason[i] <- "codon contains duplicated positions"
        next()
      }

      if (!all(diff(codon_pos) > 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon positions must be sorted increasing"
        next()
      }

      # ------------------------------------------------------------
      # amino acids

      if (!grepl(valid_amino_characters, codon_aa_string)) {
        valid[i] <- FALSE
        reason[i] <- "amino acid sequence contains invalid characters. See ?allowed_amino_acids()"
        next()
      }

      # split amino acid string into loci
      codon_aa <- strsplit(codon_aa_string, "(?<=[A-Z])(?=[^/|])", perl = TRUE)[[1]]
      codon_aa <- gsub("_", "", codon_aa)

      if (any(grepl("^[/|]|[/|]$", codon_aa))) {
        valid[i] <- FALSE
        reason[i] <- "amino acid sequence cannot start or end with / or | symbols"
        next()
      }

      if (any(grepl("[/|]{2,}", codon_aa))) {
        valid[i] <- FALSE
        reason[i] <- "amino acid sequence cannot contain adjacent / or | symbols"
        next()
      }

      if (length(codon_pos) != length(codon_aa)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("number of amino acid loci (%s) must equal the number of codon positions (%s)",
                             length(codon_aa), length(codon_pos))
        next()
      }

      # keep track of which loci contain phased or unphased hets
      codon_unphased <- grepl("/", codon_aa)
      codon_phased <- grepl("\\|", codon_aa)

      if (any(codon_unphased & codon_phased)) {
        valid[i] <- FALSE
        reason[i] <- str_c("a single amino acid locus cannot contain both the / and ",
                           "the | symbols. Heterozygous loci must be either entirely phased or ",
                           "entirely unphased")
        next()
      }

      # split into distinct amino acids
      aa_list <- strsplit(codon_aa, "[/|]")

      if (any(mapply(function(x) any(duplicated(x)), aa_list))) {
        valid[i] <- FALSE
        reason[i] <- "the same amino acid cannot be present more than once in a heterozygous call"
        next()
      }

      # count the number of amino acids at each locus
      n_aa <- mapply(length, aa_list)

      # checks on phased heterozygous loci (if present)
      if (any(codon_phased)) {

        if (length(unique(n_aa[codon_phased])) != 1) {
          valid[i] <- FALSE
          reason[i] <- str_c("if there are phased heterozygous loci within a gene then the number of ",
                             "distinct amino acids must be the same over all of these loci")
          next()
        }

        # store number of distinct amino acids in phased calls for this
        # gene. This will be tested for consistency across genes
        n_phased_vec[j] <- n_aa[codon_phased][1]
      }

      # ------------------------------------------------------------
      # read counts (if present)

      reads_present_vec[j] <- !is.na(codon_reads_string)
      if (!is.na(codon_reads_string)) {

        if (!grepl("^[0-9_/|]+$", codon_reads_string)) {
          valid[i] <- FALSE
          reason[i] <- sprintf("read counts %s contains invalid characters", codon_reads_string)
          next()
        }

        # split reads string into loci
        codon_reads <- strsplit(codon_reads_string, "_")[[1]]

        if (any(grepl("^[/|]|[/|]$", codon_reads))) {
          valid[i] <- FALSE
          reason[i] <- "read counts cannot start or end with / or | symbols"
          next()
        }

        if (any(grepl("[/|]{2,}", codon_reads))) {
          valid[i] <- FALSE
          reason[i] <- "read counts cannot contain adjacent / or | symbols"
          next()
        }

        if (length(codon_reads) != length(codon_pos)) {
          valid[i] <- FALSE
          reason[i] <- sprintf("number of read counts (%s) must equal the number of codon positions (%s)",
                               length(codon_reads), length(codon_pos))
          next()
        }

        # keep track of which loci contain phased or unphased hets
        reads_unphased <- grepl("/", codon_reads)
        reads_phased <- grepl("\\|", codon_reads)

        if (any(reads_unphased & reads_phased)) {
          valid[i] <- FALSE
          reason[i] <- str_c("read counts cannot contain both the / and ",
                             "the | symbols. Heterozygous loci must be either entirely phased or ",
                             "entirely unphased")
          next()
        }

        if (any(reads_unphased != codon_unphased) | any(reads_phased != codon_phased)) {
          valid[i] <- FALSE
          reason[i] <- str_c("if there heterozygous loci then these must correspond between amino acids and read counts. ",
                             "Also phased heterozygous loci must match to phased, and unphased to unphased")
          next()
        }

        # split into distinct read counts and make numeric
        reads_list <- mapply(as.numeric, strsplit(codon_reads, "[/|]"), SIMPLIFY = FALSE)

        if (any(unlist(reads_list) <= 0)) {
          valid[i] <- FALSE
          reason[i] <- "read counts must be greater than 0"
          next()
        }

        # count the number of read counts at each locus
        n_reads <- mapply(length, reads_list)

        if (any(n_reads != n_aa)) {
          valid[i] <- FALSE
          reason[i] <- "number of read counts must match the number of alleles at every locus"
          next()
        }

      }

    } # end loop over genes

    if (valid[i] && (length(unique(na.omit(n_phased_vec))) > 1)) {
      valid[i] <- FALSE
      reason[i] <- str_c("if there are phased heterozygous loci then the same number of distinct amino ",
                         "acids must be present at all of these loci AND over all genes")
      next()
    }

    if (valid[i] && any(duplicated(gene_name_vec))) {
      valid[i] <- FALSE
      reason[i] <- "duplicated gene names"
      next()
    }


    if (valid[i] && (length(unique(reads_present_vec)) > 1)) {
      valid[i] <- FALSE
      reason[i] <- "if read counts are present, they must be present in all genes within a sample"
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
  # - contains exactly two elements per gene
  # - gene name (first element):
  #     - contains only allowed characters
  #     - no duplicated gene names
  # - codon position (second element):
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
    gene_name_vec <- rep(NA, n_genes)

    # loop over genes
    for (j in seq_along(s1)) {

      # split this gene into elements
      s2 <- stringr::str_split(s1[j], ":", n = Inf, simplify = FALSE)[[1]]
      n_elements <- length(s2)

      if (n_elements != 2) {
        valid[i] <- FALSE
        reason[i] <- "genes must contain two elements separated by a colon (name, position)"
        next()
      }

      gene_name <- s2[1]
      codon_pos_string <- s2[2]

      # ------------------------------------------------------------
      # check gene name

      if (!grepl("^[a-zA-Z][a-zA-Z0-9_-]*$", gene_name)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene name %s contains invalid characters", gene_name)
        next()
      }

      # store gene name. This will be tested for duplication across genes
      gene_name_vec[j] <- gene_name

      # ------------------------------------------------------------
      # check codon positions

      if (!grepl("^[0-9_]+$", codon_pos_string)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("codon position %s contains invalid characters", codon_pos_string)
        next()
      }

      # split codon position into separate loci and make numeric
      codon_pos <- stringr::str_split(codon_pos_string, "_", n = Inf, simplify = FALSE)[[1]] |>
        as.numeric()

      if (any(codon_pos == 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon position cannot be zero"
        next()
      }

      if (any(duplicated(codon_pos))) {
        valid[i] <- FALSE
        reason[i] <- "codon contains duplicated positions"
        next()
      }

      if (!all(diff(codon_pos) > 0)) {
        valid[i] <- FALSE
        reason[i] <- "codon positions must be sorted increasing"
        next()
      }

    } # end loop over genes

    if (any(duplicated(gene_name_vec))) {
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
#' @title Expand variant strings into long form data.frames
#'
#' @description
#' Takes a vector of variant strings and expands into a list of data.frames
#' containing the same information in long form.
#'
#' @param x a vector of variant strings.
#'
#' @import dplyr
#'
#' @export

variant_to_long <- function(x) {

  # checks
  check_variant_string(x)

  # expand each element into a list
  ret <- mapply(function(s1) {
    mapply(function(s2) {

      # extract gene name
      gene_name <- s2[1]

      # extract codon positions to vector
      pos_vec <- strsplit(s2[2], "_")[[1]] |>
        as.numeric()

      # extract amino acids to list
      aa_list <- mapply(function(s3) {
        strsplit(s3, "[/|]")[[1]]
      }, strsplit(gsub("_", "", s2[3]), "(?<=[A-Z])(?=[^/|])", perl = TRUE)[[1]], SIMPLIFY = FALSE)
      names(aa_list) <- NULL

      # get number of amino acids at each locus
      n_aa <- mapply(length, aa_list)

      # get if each locus is phased
      is_phased <- mapply(function(s3) {
        grepl("\\|", s3)
      }, strsplit(gsub("_", "", s2[3]), "(?<=[A-Z])(?=[^/|])", perl = TRUE)[[1]])
      names(is_phased) <- NULL

      # make data.frame
      ret <- data.frame(gene = gene_name,
                        pos = rep(pos_vec, times = n_aa),
                        n_aa = rep(n_aa, times = n_aa),
                        het = rep(n_aa > 1, times = n_aa),
                        phased = rep(is_phased, times = n_aa),
                        aa = unlist(aa_list),
                        read_count = NA)

      # optionally append read counts
      if (length(s2) == 4) {

        # extract read counts to list
        read_list <- mapply(function(s3) {
          strsplit(s3, "[/|]")[[1]] |>
            as.numeric()
        }, strsplit(s2[4], "_")[[1]], SIMPLIFY = FALSE)
        names(read_list) <- NULL

        ret$read_count <- unlist(read_list)
      }

      return(ret)
    }, strsplit(s1, ":"), SIMPLIFY = FALSE) |>
      dplyr::bind_rows()
  }, strsplit(x, ";"), SIMPLIFY = FALSE)

  return(ret)
}

#------------------------------------------------
#' @title Take long form information and convert to variant string
#'
#' @description
#' Takes a list of data.frames in long form and converts each to variant string
#' format.
#'
#' @param x a list of data.frames.
#'
#' @export

long_to_variant <- function(x) {

  # check input is list
  stopifnot(is.list(x))

  ret <- mapply(function(y1) {

    # check inputs
    stopifnot(is.data.frame(y1))
    stopifnot(all(names(y1) == c("gene", "pos", "n_aa", "het", "phased", "aa", "read_count")))

    # get into string format within each gene
    df_gene <- y1 |>
      group_by(.data$gene, .data$pos) |>
      summarise(aa = ifelse(.data$phased[1],
                            paste(.data$aa, collapse = "|"),
                            paste(.data$aa, collapse = "/")),
                .groups = "drop") |>
      group_by(.data$gene) |>
      summarise(pos = paste(.data$pos, collapse = "_"),
                aa = paste(.data$aa, collapse = "_"),
                .groups = "drop")

    # behave differently depending on if read count information is available
    ret <- NA
    if (all(is.na(y1$read_count))) { # no read count info

      # concatenate into one string
      ret <- df_gene |>
        mutate(variant = sprintf("%s:%s:%s", .data$gene, .data$pos, .data$aa)) |>
        pull(.data$variant) |>
        paste(collapse = ";")

    } else if (!any(is.na(y1$read_count))) { # yes read count info

      # get read counts corresponding to each position
      df_read_count <- y1 |>
        group_by(.data$gene, .data$pos) |>
        summarise(read_count = ifelse(.data$phased[1],
                                      paste(.data$read_count, collapse = "|"),
                                      paste(.data$read_count, collapse = "/")),
                  .groups = "drop") |>
        group_by(.data$gene) |>
        summarise(pos = paste(.data$pos, collapse = "_"),
                  read_count = paste(.data$read_count, collapse = "_"),
                  .groups = "drop")

      # concatenate into one string
      ret <- df_gene |>
        left_join(df_read_count, by = join_by(.data$gene, .data$pos)) |>
        mutate(variant = sprintf("%s:%s:%s:%s", .data$gene, .data$pos, .data$aa, .data$read_count)) |>
        pull(.data$variant) |>
        paste(collapse = ";")

    } else { # invalid read count info
      stop("")
    }

    return(ret)
  }, x, SIMPLIFY = FALSE) |>
    unlist()

  # final checks on return object
  check_variant_string(ret)

  return(ret)
}

#------------------------------------------------
#' @title Expand position strings into long form data.frames
#'
#' @description
#' Takes a vector of position strings and expands into a list of data.frames
#' containing the same information in long form.
#'
#' @param x a vector of position strings.
#'
#' @import dplyr
#'
#' @export

position_to_long <- function(x) {

  # checks
  check_position_string(x)

  # expand each element into a list
  ret <- mapply(function(s1) {
    mapply(function(s2) {
      data.frame(gene = s2[1],
                 pos = as.numeric(strsplit(s2[2], "_")[[1]]))
    }, strsplit(s1, ":"), SIMPLIFY = FALSE) |>
      dplyr::bind_rows()
  }, strsplit(x, ";"), SIMPLIFY = FALSE)

  return(ret)
}

#------------------------------------------------
#' @title Take long form information and convert to position string
#'
#' @description
#' Takes a list of data.frames in long form and converts each to position string
#' format.
#'
#' @param x a list of data.frames.
#'
#' @export

long_to_position <- function(x) {

  # check input is list
  stopifnot(is.list(x))

  ret <- mapply(function(y1) {

    # check inputs
    stopifnot(is.data.frame(y1))
    stopifnot(all(names(y1) == c("gene", "pos")))

    # get into string format within each gene
    y1 |>
      group_by(.data$gene) |>
      summarise(pos = paste(.data$pos, collapse = "_"),
                .groups = "drop") |>
      mutate(variant = sprintf("%s:%s", .data$gene, .data$pos)) |>
      pull(.data$variant) |>
      paste(collapse = ";")

  }, x, SIMPLIFY = FALSE) |>
    unlist()

  # final checks on return object
  check_position_string(ret)

  return(ret)
}

#------------------------------------------------
#' @title Extract a position string from a variant string
#'
#' @description
#' Extract a position string from a variant string by stripping the amino acids.
#'
#' @param x a character string or vector of character strings.
#'
#' @import dplyr
#'
#' @export

position_from_variant_string <- function(x) {

  # checks
  check_variant_string(x)

  mapply(function(y1) {
    y1 |>
      group_by(.data$gene) |>
      reframe(pos = unique(.data$pos))
  }, variant_to_long(x), SIMPLIFY = FALSE) |>
    long_to_position()

}

#------------------------------------------------
#' @title Reorders a variant string
#'
#' @description
#' Reorders a variant string in alphabetical order of genes, and then
#' alphabetical order of amino acids at each heterozygous locus. This can be
#' useful when checking for duplicated strings as the same information may be
#' presented in a different order.
#'
#' @param x a variant string or vector of variant strings.
#'
#' @export

order_variant_string <- function(x) {

  # checks
  check_variant_string(x)

  mapply(function(y) {
    arrange(y, .data$gene, .data$pos, .data$aa)
  }, variant_to_long(x), SIMPLIFY = FALSE) |>
    long_to_variant()

}

#------------------------------------------------
#' @title Reorders a position string
#'
#' @description
#' Reorders a position string in alphabetical order of genes. This can be useful
#' when checking for duplicated strings as the same information may be presented
#' in a different order.
#'
#' @param x a position string or vector of position strings.
#'
#' @export

order_position_string <- function(x) {

  # checks
  check_position_string(x)

  mapply(function(y) {
    arrange(y, .data$gene, .data$pos)
  }, position_to_long(x), SIMPLIFY = FALSE) |>
    long_to_position()

}

#------------------------------------------------
#' @title Count the number of heterozygous loci in each variant string
#'
#' @description
#' Count the number of heterozygous loci in each variant string and return as a
#' vector.
#'
#' @param x a variant string or vector of variant strings.
#'
#' @export

count_het_loci <- function(x) {

  # checks
  check_variant_string(x)

  mapply(function(y1) {
    y1 |>
      group_by(.data$gene, .data$pos) |>
      summarise(het = .data$het[1],
                .groups = "drop") |>
      pull(.data$het) |>
      sum()
  }, variant_to_long(x))
}

#------------------------------------------------
#' @title Compares variant strings to look for a match
#'
#' @description
#' Compares a target variant string against a vector of comparison strings. A
#' match is found if every amino acid at every codon position in every gene of
#' the target is also found within the comparison. Note that ambiguous matches
#' may occur if there are multiple heterozygous loci in the comparison. In this
#' case, the target may or may not be within this sample. A match is recorded
#' but a second output also flags this as an ambiguous match.
#'
#' @param target_string a single variant string that we want to compare. Cannot
#'   contain any heterozygous calls.
#' @param comparison_strings a vector of variant strings against which the
#'   target is compared.
#'
#' @import dplyr
#'
#' @export

compare_variant_string <- function(target_string, comparison_strings) {

  # checks
  check_variant_string(target_string)
  stopifnot(length(target_string) == 1)
  if (count_het_loci(target_string) != 0) {
    stop("target string cannot contain any heterozygous loci")
  }
  check_variant_string(comparison_strings)

  # get target in long form
  df_target <- variant_to_long(target_string)[[1]]

  # loop over all comparison strings
  n <- length(comparison_strings)
  ret <- data.frame(match = rep(NA, n),
                    ambiguous = NA,
                    prop = NA)
  for (i in 1:n) {

    # get this comparison string in long form
    df_comparison <- variant_to_long(comparison_strings[i])[[1]]

    # get if there is a match and if in a het for each locus
    df_match <- mapply(function(j) {
      match_pos <- (df_target$gene[j] == df_comparison$gene) &
        (df_target$pos[j] == df_comparison$pos)
      match_all <- match_pos & (df_target$aa[j] == df_comparison$aa)

      c(match = sum(match_all),
        het = sum(df_comparison$het * match_all),
        read_num = sum(df_comparison$read_count * match_all),
        read_denom = sum(df_comparison$read_count * match_pos))
    }, 1:nrow(df_target), SIMPLIFY = FALSE) |>
      dplyr::bind_rows()

    # get if a match over all loci, and if match is ambiguous
    ret$match[i] <- all(df_match$match == 1)
    ret$ambiguous[i] <- (sum(df_match$het) > 1)

    if (sum(df_match$het) == 0) {
      ret$prop[i] <- ret$match[i]
    } else if (sum(df_match$het) == 1) {
      w <- which(df_match$het == 1)
      ret$prop[i] <- df_match$read_num[w] / df_match$read_denom[w]
    } else {
      ret$prop[i] <- NA
    }
  }

  return(ret)
}

#------------------------------------------------
#' @title Compares a position strings against variant strings to look for a match
#'
#' @description
#' Compares a target position string against a vector of comparison strings. A
#' match is found if every codon position in every gene of the target is also
#' found within the comparison (irrespective of the observed amino acids).
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

  # get target in long form
  df_target <- position_to_long(target_string)[[1]]

  # get comparisons in long form
  list_comparison <- position_from_variant_string(comparison_strings) |>
    position_to_long()

  # loop over all comparison strings
  n <- length(comparison_strings)
  ret <- rep(NA, n)
  for (i in 1:n) {

    ret[i] <- mapply(function(j) {
      any((df_target$gene[j] == list_comparison[[i]]$gene) &
        (df_target$pos[j] == list_comparison[[i]]$pos))
    }, 1:nrow(df_target), SIMPLIFY = FALSE) |>
      unlist() |>
      all()
  }

  return(ret)
}

#------------------------------------------------
#' @title Extract all single-locus variants from a variant string
#'
#' @description
#' Takes a vector of variant strings, potentially with information at multiple
#' codon positions or genes, and returns variant strings corresponding to all
#' unique single-locus variants within the input. For example, crt:72_73:C_N/V
#' can be extracted to crt:72:C, crt:73:N, and crt:73:V.
#'
#' @param x a vector of variant strings.
#'
#' @export

extract_single_locus_variants <- function(x) {

  # checks
  check_variant_string(x)

  mapply(function(y) {
    sprintf("%s:%s:%s", y$gene, y$pos, y$aa)
  }, variant_to_long(x), SIMPLIFY = FALSE)

}

#------------------------------------------------
#' @title Get all genotypes that are consistent with a variant string
#'
#' @description
#' For a variant string with at most one heterozygous locus we can unambiguously
#' define the genotypes that are present in this mixture. This function returns
#' all such component genotypes.
#'
#' @param x a vector of variant strings.
#'
#' @importFrom tidyr pivot_longer
#' @import dplyr
#'
#' @export

get_consistent_variants <- function(x) {

  # check
  check_variant_string(x)

  # get number of het loci. We will focus on unambiguous, so at most 1 het locus
  n_het_loci <- count_het_loci(x)

  # get all variants into long format
  variant_list <- variant_to_long(x)

  n <- length(x)
  ret <- replicate(n, NA, simplify = FALSE)
  for (i in 1:n) {
    if (n_het_loci[i] == 0) {
      ret[[i]] <- long_to_variant(variant_list[i])
    } else if (n_het_loci[i] == 1) {

      # get amino acids into list nested within data.frame
      df_nested <- variant_list[[i]] |>
        group_by(.data$gene, .data$pos) |>
        summarise(aa = list(.data$aa),
                  .groups = "drop")

      # get all possible combinations of amino acids
      combos <- expand.grid(df_nested$aa) |>
        as.matrix() |>
        t()
      colnames(combos) <- sprintf("combo_%s", 1:ncol(combos))

      # get all possible variant strings
      components <- df_nested |>
        dplyr::bind_cols(combos) |>
        select(-.data$aa) |>
        pivot_longer(cols = starts_with("combo_"),
                     names_to = "combo",
                     values_to = "aa") |>
        group_by(.data$combo, .data$gene) |>
        summarise(pos = paste(.data$pos, collapse = "_"),
                  aa = paste(.data$aa, collapse = "_"),
                  .groups = "drop") |>
        mutate(variant = sprintf("%s:%s:%s", .data$gene, .data$pos, .data$aa)) |>
        group_by(.data$combo) |>
        summarise(variant = paste(.data$variant, collapse = ";"),
                  .groups = "drop") |>
        pull(.data$variant)

      ret[[i]] <- components
    }
  }

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
