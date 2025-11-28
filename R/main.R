
# avoids "no visible bindings" warnings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("gene", "pos", "aa", "het", "phased", "read_count",
                           "combo", "variant"))
}

# --- FUNCTION LIST ---
# check_variant_string
# check_position_string
# variant_to_long
# long_to_variant
# position_to_long
# long_to_position
# position_from_variant_string
# subset_position
# order_variant_string
# order_position_string
# count_unphased_hets
# count_phased_hets
# drop_read_counts
# compare_variant_string
# compare_position_string
# extract_single_locus_variants
# get_component_variants
# allowed_amino_acids

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
  # - contains exactly three or four elements per gene
  # - gene name (first element):
  #     - contains only allowed characters
  # - codon position (second element):
  #     - contains only allowed characters
  #     - if there are hyphens denoting a range, checks that a single internal hyphen only
  #     - positions are non-zero
  #     - positions are not duplicated
  #     - positions are sorted increasing
  # - amino acid (third element):
  #     - contains only allowed characters
  #     - cannot start or end with / or | symbols
  #     - cannot contain adjacent / or | symbols
  #     - number of amino acid loci matches number of codon positions
  #     - a single amino acid group cannot contain both | and / (no mixed phasing)
  #     - unphased hets cannot have the same amino acid multiple times (phased hets can)
  #     - all amino acid loci with phased hets contain the same number of amino acids within a gene
  # - read count (optional fourth element):
  #     - contains only allowed characters
  #     - cannot start or end with / or | symbols
  #     - cannot contain adjacent / or | symbols
  #     - number of read count loci matches number of codon positions
  #     - a single read count group cannot contain both | and / (no mixed phasing)
  #     - the pattern of phased and unphased read counts must match that in amino acids
  #     - read counts must be greater than 0
  #     - number of read counts at each locus must match number of amino acids at each locus
  # - comparison across all elements
  #     - no duplicated gene names
  #     - if read counts are present they must be present in all genes, and vice versa
  #     - all amino acid loci with phased hets contain the same number of amino acids between genes

  # is character vector
  stopifnot(all(is.character(x)))

  # only work on unique elements in x
  x_unique <- unique(x)
  n_unique <- length(x_unique)

  # create objects for storing which rows fail. This allows for more
  # informative error messages in which several issues can be listed at once
  # rather than breaking at the first issue
  valid <- rep(TRUE, n_unique)
  reason <- rep(NA, n_unique)

  # get list of valid amino acid characters. Do this once here to avoid
  # repetition for every element of x
  #IUPAC_df <- allowed_amino_acids()
  #valid_amino_characters <- paste0("^[", paste(IUPAC_df$IUPAC_amino_acid_code, collapse = ""), "_/|]+$")
  valid_amino_characters <- "^[ACDEFGHIKLMNPQRSTVWY_/|]+$"

  # input may be a single string or a vector of strings. Loop through elements
  for (i in 1:n_unique) {

    # split into genes
    s1 <- stringr::str_split(x_unique[i], ";", n = Inf, simplify = FALSE)[[1]]
    n_genes <- length(s1)

    if (any(s1 == "")) {
      valid[i] <- FALSE
      reason[i] <-"contains one or more empty genes"
      next()
    }

    # keep track of elements for comparison across all genes
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

      if (!grepl("^[a-zA-Z][a-zA-Z0-9._-]*$", gene_name)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("gene name %s contains invalid characters", gene_name)
        next()
      }

      # store gene name. This will be tested for duplication across genes
      gene_name_vec[j] <- gene_name

      # ------------------------------------------------------------
      # check codon positions
      if (!grepl("^[0-9_]$|^[0-9_][0-9_-]*[0-9_]$", codon_pos_string)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("codon position %s contains invalid characters", codon_pos_string)
        next()
      }

      # split codon position into separate loci
      codon_pos <- stringr::str_split(codon_pos_string, "_", n = Inf, simplify = FALSE)[[1]]

      # make numeric, dealing with ranges where needed
      if (any(grepl("-", codon_pos))) {
        if (any(stringr::str_count(codon_pos, "-") > 1)) {
          valid[i] <- FALSE
          reason[i] <- sprintf("codon position %s does not define a valid range", codon_pos)
          next()
        }
        codon_pos <- mapply(function(x) {
          r <- as.numeric(strsplit(x, "-")[[1]])
          if (length(r) > 1) {
            r <- r[1]:r[2]
          }
          r
        }, codon_pos, SIMPLIFY = FALSE) |>
          unlist()
      }
      codon_pos <- as.numeric(codon_pos)

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
        reason[i] <- "amino acid sequence contains invalid characters. See ?variantstring::allowed_amino_acids()"
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

      if (any(mapply(function(x) any(duplicated(x)), aa_list[which(!codon_phased)]))) {
        valid[i] <- FALSE
        reason[i] <- str_c("the same amino acid cannot be present more than once in an ",
                           "unphased heterozygous call")
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

    # ------------------------------------------------------------
    # comparison over genes

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
        message(sprintf("  - %s: %s", x_unique[i], reason[i]))
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
  #     - if there are hyphens denoting a range, a single internal hyphen only
  #     - positions are non-zero
  #     - positions are not duplicated
  #     - positions are sorted increasing

  # is character vector
  stopifnot(all(is.character(x)))

  # only work on unique elements in x
  x_unique <- unique(x)
  n_unique <- length(x_unique)

  # create objects for storing which rows fail. This allows for more
  # informative error messages in which several issues can be listed at once
  # rather than breaking at the first issue
  valid <- rep(TRUE, n_unique)
  reason <- rep(NA, n_unique)

  # input may be a single string or a vector of strings. Loop through elements
  for (i in 1:n_unique) {

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

      if (!grepl("^[a-zA-Z][a-zA-Z0-9._-]*$", gene_name)) {
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
      codon_pos <- stringr::str_split(codon_pos_string, "_", n = Inf, simplify = FALSE)[[1]]

      # make numeric, dealing with ranges where needed
      if (any(grepl("-", codon_pos))) {
        if (any(stringr::str_count(codon_pos, "-") > 1)) {
          valid[i] <- FALSE
          reason[i] <- sprintf("codon position %s does not define a valid range", codon_pos)
          next()
        }
        codon_pos <- mapply(function(x) {
          r <- as.numeric(strsplit(x, "-")[[1]])
          if (length(r) > 1) {
            r <- r[1]:r[2]
          }
          r
        }, codon_pos, SIMPLIFY = FALSE) |>
          unlist()
      }
      codon_pos <- as.numeric(codon_pos)

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
        message(sprintf("  - %s: %s", x_unique[i], reason[i]))
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
#' @param x a vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @import dplyr
#'
#' @export

variant_to_long <- function(x) {

  # basic checks
  stopifnot(is.character(x))

  # work on unique elements in x and map back at end
  x_unique <- unique(x)

  # expand each element into a list
  l_unique <- mapply(function(s1) {
    mapply(function(s2) {

      # extract gene name
      gene_name <- s2[1]

      # extract codon positions to vector
      pos_vec <- strsplit(s2[2], "_")[[1]]
      if (any(grepl("-", pos_vec))) {
        pos_vec <- mapply(function(x) {
          r <- as.numeric(x)
          if (length(r) > 1) {
            r <- r[1]:r[2]
          }
          r
        }, strsplit(pos_vec, "-"), SIMPLIFY = FALSE) |>
          unlist()
      }
      pos_vec <- as.numeric(pos_vec)

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
  }, strsplit(x_unique, ";"), SIMPLIFY = FALSE)

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Take long form information and convert to variant string
#'
#' @description
#' Takes a list of data frames in long form and converts each to variant string
#' format.
#'
#' @param x a list of data frames corresponding to variant strings.
#'
#' @export

long_to_variant <- function(x) {

  # check input is list
  stopifnot(is.list(x))

  # work on unique elements in x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y1) {

    # check inputs
    stopifnot(is.data.frame(y1))
    stopifnot(all(names(y1) == c("gene", "pos", "n_aa", "het", "phased", "aa", "read_count")))

    # get into string format within each gene
    df_gene <- y1 |>
      group_by(gene, pos) |>
      summarise(aa = ifelse(phased[1],
                            paste(aa, collapse = "|"),
                            paste(aa, collapse = "/")),
                .groups = "drop") |>
      group_by(gene) |>
      summarise(pos = paste(pos, collapse = "_"),
                aa = paste(aa, collapse = "_"),
                .groups = "drop")

    # behave differently depending on if read count information is available
    ret <- NA
    if (all(is.na(y1$read_count))) { # no read count info

      # concatenate into one string
      ret <- df_gene |>
        mutate(variant = sprintf("%s:%s:%s", gene, pos, aa)) |>
        pull(variant) |>
        paste(collapse = ";")

    } else if (!any(is.na(y1$read_count))) { # yes read count info

      # get read counts corresponding to each position
      df_read_count <- y1 |>
        group_by(gene, pos) |>
        summarise(read_count = ifelse(phased[1],
                                      paste(read_count, collapse = "|"),
                                      paste(read_count, collapse = "/")),
                  .groups = "drop") |>
        group_by(gene) |>
        summarise(pos = paste(pos, collapse = "_"),
                  read_count = paste(read_count, collapse = "_"),
                  .groups = "drop")

      # concatenate into one string
      ret <- df_gene |>
        left_join(df_read_count, by = c("gene", "pos")) |>
        mutate(variant = sprintf("%s:%s:%s:%s", gene, pos, aa, read_count)) |>
        pull(variant) |>
        paste(collapse = ";")

    } else { # invalid read count info
      stop("")
    }

    return(ret)
  }, x_unique, SIMPLIFY = FALSE) |>
    unlist()

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Expand position strings into long form data.frames
#'
#' @description
#' Takes a vector of position strings and expands into a list of data.frames
#' containing the same information in long form.
#'
#' @param x a vector of position strings. Note, these are not internally checked for being valid position strings, it is up to the user to ensure this (see \code{?check_position_string}).
#'
#' @import dplyr
#'
#' @export

position_to_long <- function(x) {

  # basic checks
  stopifnot(is.character(x))

  # work on unique elements in x and map back at end
  x_unique <- unique(x)

  # expand each element into a list
  l_unique <- mapply(function(s1) {
    mapply(function(s2) {

      # extract codon positions to vector
      pos_vec <- strsplit(s2[2], "_")[[1]]
      if (any(grepl("-", pos_vec))) {
        pos_vec <- mapply(function(x) {
          r <- as.numeric(x)
          if (length(r) > 1) {
            r <- r[1]:r[2]
          }
          r
        }, strsplit(pos_vec, "-"), SIMPLIFY = FALSE) |>
          unlist()
      }
      pos_vec <- as.numeric(pos_vec)

      data.frame(gene = s2[1],
                 pos = pos_vec)
    }, strsplit(s1, ":"), SIMPLIFY = FALSE) |>
      dplyr::bind_rows()
  }, strsplit(x_unique, ";"), SIMPLIFY = FALSE)

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Take long form information and convert to position string
#'
#' @description
#' Takes a list of data frames in long form and converts each to position string
#' format.
#'
#' @param x a list of data frames corresponding to position strings.
#'
#' @export

long_to_position <- function(x) {

  # check input is list
  stopifnot(is.list(x))

  # work on unique elements in x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y1) {

    # check inputs
    stopifnot(is.data.frame(y1))
    stopifnot(all(names(y1) == c("gene", "pos")))

    # get into string format within each gene
    y1 |>
      group_by(gene) |>
      summarise(pos = paste(pos, collapse = "_"),
                .groups = "drop") |>
      mutate(variant = sprintf("%s:%s", gene, pos)) |>
      pull(variant) |>
      paste(collapse = ";")

  }, x_unique, SIMPLIFY = FALSE) |>
    unlist()

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Extract a position string from a variant string
#'
#' @description
#' Extract a position string from a variant string by stripping the amino acids.
#'
#' @param x a character string or vector of character strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @import dplyr
#'
#' @export

position_from_variant_string <- function(x) {

  # work on unique elements in x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y1) {
    y1 |>
      group_by(gene) |>
      reframe(pos = unique(pos))
  }, variant_to_long(x_unique), SIMPLIFY = FALSE) |>
    long_to_position()

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Subset position of a variant string
#'
#' @description
#' Given a vector of variant strings and a single position string, subsets all
#' variant strings to only the genes and codons in the position string. Retains
#' read counts at these positions if present.
#'
#' @param position_string a single position string. Note, these are not internally checked for being valid position strings, it is up to the user to ensure this (see \code{?check_position_string}).
#' @param variant_strings a variant string or vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @import dplyr
#'
#' @export

subset_position <- function(position_string, variant_strings) {

  # checks
  stopifnot(length(position_string) == 1)

  # get position string in long form
  df_position <- position_to_long(position_string)[[1]]

  # work on unique variant strings and map back at end
  v_unique <- unique(variant_strings)

  l_unique <- mapply(function(x) {
    df_diff <- anti_join(df_position, x, by = c("gene", "pos"))
    if (nrow(df_diff) == 0) {
      df_sub <- semi_join(x, df_position, by = c("gene", "pos"))
      ret <- long_to_variant(list(df_sub))
    } else {
      ret <- NA
    }
    ret
  }, variant_to_long(v_unique), SIMPLIFY = FALSE) |>
    unlist()

  # map uniques back to original
  ret <- l_unique[match(variant_strings, v_unique)]

  return(ret)
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
#' @param x a variant string or vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

order_variant_string <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y) {
    arrange(y, gene, pos, aa)
  }, variant_to_long(x_unique), SIMPLIFY = FALSE) |>
    long_to_variant()

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Reorders a position string
#'
#' @description
#' Reorders a position string in alphabetical order of genes. This can be useful
#' when checking for duplicated strings as the same information may be presented
#' in a different order.
#'
#' @param x a position string or vector of position strings. Note, these are not internally checked for being valid position strings, it is up to the user to ensure this (see \code{?check_position_string}).
#'
#' @export

order_position_string <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y) {
    arrange(y, gene, pos)
  }, position_to_long(x_unique), SIMPLIFY = FALSE) |>
    long_to_position()

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Count the number of unphased heterozygous loci
#'
#' @description
#' Count the number of unphased heterozygous loci in each variant string. Return
#' the number as a vector.
#'
#' @param x a variant string or vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

count_unphased_hets <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  h_unique <- mapply(function(y1) {
    y1 |>
      group_by(gene, pos) |>
      summarise(n = (het[1] == TRUE) & (phased[1] == FALSE),
                .groups = "drop") |>
      pull(n) |>
      sum()
  }, variant_to_long(x_unique))

  # map uniques back to original
  ret <- h_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Count the number of phased heterozygous loci
#'
#' @description
#' Count the number of phased heterozygous loci in each variant string. Return
#' the number as a vector.
#'
#' @param x a variant string or vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

count_phased_hets <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  h_unique <- mapply(function(y1) {
    y1 |>
      group_by(gene, pos) |>
      summarise(n = (het[1] == TRUE) & (phased[1] == TRUE),
                .groups = "drop") |>
      pull(n) |>
      sum()
  }, variant_to_long(x_unique))

  # map uniques back to original
  ret <- h_unique[match(x, x_unique)]

  return(ret)
}

#------------------------------------------------
#' @title Drop read counts from a variant string.
#'
#' @description
#' Takes a vector of variant strings and strips and information on read counts.
#'
#' @param x a variant string or vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

drop_read_counts <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y) {
    y$read_count <- NA
    y
  }, variant_to_long(x_unique), SIMPLIFY = FALSE) |>
    long_to_variant()

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]

  return(ret)
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
#'   contain any heterozygous calls. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#' @param comparison_strings a vector of variant strings against which the
#'   target is compared. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @import dplyr
#'
#' @export

compare_variant_string <- function(target_string, comparison_strings) {

  # checks
  stopifnot(length(target_string) == 1)
  if ((count_unphased_hets(target_string) > 0) || (count_phased_hets(target_string) > 0)) {
    stop("target string cannot contain any heterozygous loci")
  }

  # get target in long form
  df_target <- variant_to_long(target_string)[[1]] |>
    select(gene, pos, aa)

  # work on unique comparison_strings and map back at end
  c_unique <- unique(comparison_strings)

  # loop over all unique comparison strings
  n_unique <- length(c_unique)
  df_unique <- data.frame(match = rep(NA, n_unique),
                          ambiguous = NA,
                          prop = NA)
  for (i in 1:n_unique) {

    # get this comparison string in long form
    df_comparison <- variant_to_long(c_unique[i])[[1]]

    # compare target against comparison
    df_target_match <- df_target |>
      rowwise() |>
      summarise(gene = gene[1],
                pos = pos[1],
                match = any((gene == df_comparison$gene) &
                              (pos == df_comparison$pos) &
                              (aa == df_comparison$aa)),
                .groups = "drop") |>
      ungroup()

    # simple exit if not match
    if (!all(df_target_match$match)) {
      df_unique$match[i] <- FALSE
      df_unique$ambiguous[i] <- FALSE
      df_unique$prop[i] <- 0
      next
    }

    # compare comparison back against target
    df_comparison_match <- df_comparison |>
      rowwise(gene, pos, het, phased, aa, read_count) |>
      summarise(match = any((gene == df_target$gene) &
                              (pos == df_target$pos) &
                              (aa == df_target$aa)),
                .groups = "drop") |>
      group_by(gene, pos) |>
      summarise(het = het[1],
                phased = phased[1],
                read_count_total = sum(read_count),
                read_count = sum(read_count*match),
                match = any(match),
                .groups = "drop") |>
      filter(match)

    # count hets
    n_phased_het <- sum(df_comparison_match$het * df_comparison_match$phased)
    n_unphased_het <- sum(df_comparison_match$het * !df_comparison_match$phased)

    # determine if unambiguous match
    df_unique$match[i] <- TRUE
    df_unique$ambiguous[i] <- ((n_phased_het > 0) + n_unphased_het) > 1

    # get proportion if calculable
    df_unique$prop[i] <- NA
    if (!df_unique$ambiguous[i]) {
      if (any(is.na(df_comparison_match$read_count))) {
        if ((n_phased_het + n_unphased_het) == 0) {
          df_unique$prop[i] <- 1
        }
      } else {
        prop <- df_comparison_match$read_count / df_comparison_match$read_count_total
        if (all(prop == 1)) {
          df_unique$prop[i] <- 1
        } else {
          u_prop <- unique(prop[prop != 1])
          if (length(u_prop) == 1) {
            df_unique$prop[i] <- u_prop
          }
        }
      }
    }

  }

  # map uniques back to original
  ret <- df_unique[match(comparison_strings, c_unique),]

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
#' @param target_string a single position string that we want to compare. Note, these are not internally checked for being valid position strings, it is up to the user to ensure this (see \code{?check_position_string}).
#' @param comparison_strings a vector of variant strings against which the
#'   target is compared. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

compare_position_string <- function(target_string, comparison_strings) {

  # checks
  stopifnot(length(target_string) == 1)

  # get target in long form
  df_target <- position_to_long(target_string)[[1]]

  # work on unique comparison_strings and map back at end
  c_unique <- unique(comparison_strings)

  # get comparisons in long form
  list_comparison <- position_from_variant_string(c_unique) |>
    position_to_long()

  # loop over all unique comparison strings
  n_unique <- length(c_unique)
  ret_unique <- rep(NA, n_unique)
  for (i in 1:n_unique) {

    ret_unique[i] <- mapply(function(j) {
      any((df_target$gene[j] == list_comparison[[i]]$gene) &
        (df_target$pos[j] == list_comparison[[i]]$pos))
    }, 1:nrow(df_target), SIMPLIFY = FALSE) |>
      unlist() |>
      all()
  }

  # map uniques back to original
  ret <- ret_unique[match(comparison_strings, c_unique)]

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
#' @param x a vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

extract_single_locus_variants <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  l_unique <- mapply(function(y) {
    sprintf("%s:%s:%s", y$gene, y$pos, y$aa)
  }, variant_to_long(x_unique), SIMPLIFY = FALSE)

  # map uniques back to original
  ret <- l_unique[match(x, x_unique)]
}

#------------------------------------------------
#' @title Get all genotypes that are consistent with a variant string
#'
#' @description
#' For a variant string with at most one heterozygous locus we can unambiguously
#' define the genotypes that are present in this mixture. This function returns
#' all such component genotypes.
#'
#' @param x a vector of variant strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @importFrom tidyr pivot_longer
#' @import dplyr
#'
#' @export

get_component_variants <- function(x) {

  # work on unique x and map back at end
  x_unique <- unique(x)

  # focus on unambiguous
  unphased_hets <- count_unphased_hets(x_unique)
  phased_hets <- count_phased_hets(x_unique)

  # get all variants into long format and drop any read counts
  variant_list <- mapply(function(y) {
    y$read_count <- NA
    y
  }, variant_to_long(x_unique), SIMPLIFY = FALSE)

  n_unique <- length(x_unique)
  df_ret <- replicate(n_unique, NA, simplify = FALSE)
  for (i in 1:n_unique) {
    if ((unphased_hets[i] == 0) & (phased_hets[i] == 0)) {
      df_ret[[i]] <- long_to_variant(variant_list[i])

    } else if ((unphased_hets[i] == 1) && (phased_hets[i] == 0)) {

      # get amino acids into list nested within data.frame
      df_nested <- variant_list[[i]] |>
        group_by(gene, pos) |>
        summarise(aa = list(aa),
                  .groups = "drop")

      # get all possible combinations of amino acids
      combos <- expand.grid(df_nested$aa) |>
        as.matrix() |>
        t()
      colnames(combos) <- sprintf("combo_%s", 1:ncol(combos))

      # get all possible variant strings
      components <- df_nested |>
        dplyr::bind_cols(combos) |>
        select(-aa) |>
        pivot_longer(cols = starts_with("combo_"),
                     names_to = "combo",
                     values_to = "aa") |>
        group_by(combo, gene) |>
        summarise(pos = paste(pos, collapse = "_"),
                  aa = paste(aa, collapse = "_"),
                  .groups = "drop") |>
        mutate(variant = sprintf("%s:%s:%s", gene, pos, aa)) |>
        group_by(combo) |>
        summarise(variant = paste(variant, collapse = ";"),
                  .groups = "drop") |>
        pull(variant)

      df_ret[[i]] <- components

    } else if ((unphased_hets[i] == 0) && (phased_hets[i] > 0)) {

      # get amino acids into list nested within data.frame
      df_nested <- variant_list[[i]] |>
        group_by(gene, pos) |>
        summarise(aa = list(aa),
                  .groups = "drop")

      # get all phased combinations of amino acids
      combos <- df_nested$aa |>
        rbind.data.frame() |>
        as.matrix() |>
        t()
      colnames(combos) <- sprintf("combo_%s", 1:ncol(combos))
      rownames(combos) <- NULL

      # get all possible variant strings
      components <- df_nested |>
        select(-aa) |>
        dplyr::bind_cols(combos) |>
        pivot_longer(cols = starts_with("combo_"),
                     names_to = "combo",
                     values_to = "aa") |>
        group_by(combo, gene) |>
        summarise(pos = paste(pos, collapse = "_"),
                  aa = paste(aa, collapse = "_"),
                  .groups = "drop") |>
        mutate(variant = sprintf("%s:%s:%s", gene, pos, aa)) |>
        group_by(combo) |>
        summarise(variant = paste(variant, collapse = ";"),
                  .groups = "drop") |>
        pull(variant)

      df_ret[[i]] <- components

    } else {
      df_ret[[i]] <- NA

    }
  }

  # map uniques back to original
  ret <- df_ret[match(x, x_unique)]

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
