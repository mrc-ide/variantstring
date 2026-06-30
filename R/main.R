
# avoids "no visible bindings" warnings
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("gene", "pos", "aa", "het", "phased", "read_count",
                           "combo", "variant", "n_aa.x", "n_aa.y",
                           "het.x", "het.y",
                           "phased.x", "phased.y",
                           "aa.x", "aa.y",
                           "read_count.x", "read_count.y",
                           "n_aa"))
}

#------------------------------------------------
# Internal helpers (not exported)
#------------------------------------------------

# Apply a function to the unique values of x and map the results back to the
# original positions. This centralises the "work on unique entries, then expand
# back" pattern used throughout the package, both for speed (the function is
# only evaluated once per distinct value) and to avoid the repeated, error-prone
# boilerplate. `fun` must accept the vector of unique values and return an object
# of the same length that can be indexed by position (an atomic vector or list).
with_unique <- function(x, fun) {
  x_unique <- unique(x)
  fun(x_unique)[match(x, x_unique)]
}

# Expand a vector of codon position tokens into a numeric vector, resolving any
# hyphenated ranges (e.g. c("72", "74-76") -> c(72, 74, 75, 76)). Assumes any
# ranges are well formed (a single internal hyphen); validity is checked
# separately in check_variant_string() / check_position_string().
expand_positions <- function(tokens) {
  if (any(grepl("-", tokens))) {
    tokens <- unlist(lapply(strsplit(tokens, "-"), function(r) {
      r <- as.numeric(r)
      if (length(r) > 1) r[1]:r[2] else r
    }))
  }
  as.numeric(tokens)
}

# Split a gene's amino-acid string into a vector of per-locus calls, coping with
# both the underscore-separated and concise notations (e.g. "C_V_K/T" and
# "CVK/T" both give c("C", "V", "K/T")).
split_aa_loci <- function(aa_string) {
  strsplit(gsub("_", "", aa_string), "(?<=[A-Z])(?=[^/|])", perl = TRUE)[[1]]
}

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
  for (i in seq_len(n_unique)) {

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

      # ranges must contain a single internal hyphen only
      if (any(stringr::str_count(codon_pos, "-") > 1)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("codon position %s does not define a valid range", codon_pos)
        next()
      }

      # make numeric, expanding any ranges
      codon_pos <- expand_positions(codon_pos)

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
      codon_aa <- split_aa_loci(codon_aa_string)

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

  # if any invalid entries then collate all issues into the error message
  if (any(!valid)) {
    issues <- sprintf("  - %s: %s", x_unique[!valid], reason[!valid])
    stop(paste0(c("The following issues were found:", issues), collapse = "\n"),
         call. = FALSE)
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
  for (i in seq_len(n_unique)) {

    # split into genes
    s1 <- stringr::str_split(x_unique[i], ";", n = Inf, simplify = FALSE)[[1]]
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

      # ranges must contain a single internal hyphen only
      if (any(stringr::str_count(codon_pos, "-") > 1)) {
        valid[i] <- FALSE
        reason[i] <- sprintf("codon position %s does not define a valid range", codon_pos)
        next()
      }

      # make numeric, expanding any ranges
      codon_pos <- expand_positions(codon_pos)

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

  # if any invalid entries then collate all issues into the error message
  if (any(!valid)) {
    issues <- sprintf("  - %s: %s", x_unique[!valid], reason[!valid])
    stop(paste0(c("The following issues were found:", issues), collapse = "\n"),
         call. = FALSE)
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

  with_unique(x, function(x_unique) {
    lapply(strsplit(x_unique, ";"), function(genes) {
      dplyr::bind_rows(lapply(strsplit(genes, ":"), function(s2) {

        # positions (resolving any ranges) and amino acid loci
        pos_vec <- expand_positions(strsplit(s2[2], "_")[[1]])
        aa_loci <- split_aa_loci(s2[3])
        aa_list <- strsplit(aa_loci, "[/|]")
        n_aa <- lengths(aa_list)
        is_phased <- grepl("\\|", aa_loci)

        # make data.frame, one row per amino acid
        ret <- data.frame(gene = s2[1],
                          pos = rep(pos_vec, times = n_aa),
                          n_aa = rep(n_aa, times = n_aa),
                          het = rep(n_aa > 1, times = n_aa),
                          phased = rep(is_phased, times = n_aa),
                          aa = unlist(aa_list),
                          read_count = NA)

        # optionally append read counts, mirroring the amino acid format
        if (length(s2) == 4) {
          ret$read_count <- unlist(lapply(strsplit(s2[4], "_")[[1]], function(s3) {
            as.numeric(strsplit(s3, "[/|]")[[1]])
          }))
        }

        ret
      }))
    })
  })
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

  with_unique(x, function(x_unique) {
    vapply(x_unique, function(y1) {

      # check inputs
      stopifnot(is.data.frame(y1))
      stopifnot(all(names(y1) == c("gene", "pos", "n_aa", "het", "phased", "aa", "read_count")))

      # read counts must be either entirely present or entirely absent
      if (all(is.na(y1$read_count))) {
        has_reads <- FALSE
      } else if (!any(is.na(y1$read_count))) {
        has_reads <- TRUE
      } else {
        stop("read counts must be either present at every locus or absent at every locus")
      }

      # build each gene in alphabetical order (genes are sorted before any
      # manipulation, matching the rest of the package)
      gene_strings <- vapply(sort(unique(y1$gene), method = "radix"), function(g) {
        sub <- y1[y1$gene == g, , drop = FALSE]
        upos <- sort(unique(sub$pos))

        # collapse a column to one entry per locus, joining alleles with | when
        # phased and / otherwise
        collapse_locus <- function(col) {
          vapply(upos, function(p) {
            rows <- sub[sub$pos == p, , drop = FALSE]
            paste(rows[[col]], collapse = if (isTRUE(rows$phased[1])) "|" else "/")
          }, character(1))
        }

        out <- sprintf("%s:%s:%s", g, paste(upos, collapse = "_"),
                       paste(collapse_locus("aa"), collapse = "_"))
        if (has_reads) {
          out <- sprintf("%s:%s", out, paste(collapse_locus("read_count"), collapse = "_"))
        }
        out
      }, character(1))

      paste(gene_strings, collapse = ";")
    }, character(1))
  })
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

  with_unique(x, function(x_unique) {
    lapply(strsplit(x_unique, ";"), function(genes) {
      dplyr::bind_rows(lapply(strsplit(genes, ":"), function(s2) {
        data.frame(gene = s2[1],
                   pos = expand_positions(strsplit(s2[2], "_")[[1]]))
      }))
    })
  })
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

  with_unique(x, function(x_unique) {
    vapply(x_unique, function(y1) {

      # check inputs
      stopifnot(is.data.frame(y1))
      stopifnot(all(names(y1) == c("gene", "pos")))

      # build each gene in alphabetical order, joining positions with _
      gene_strings <- vapply(sort(unique(y1$gene), method = "radix"), function(g) {
        sprintf("%s:%s", g, paste(y1$pos[y1$gene == g], collapse = "_"))
      }, character(1))

      paste(gene_strings, collapse = ";")
    }, character(1))
  })
}

#------------------------------------------------
#' @title Extract a position string from a variant string
#'
#' @description
#' Extract a position string from a variant string by stripping the amino acids.
#'
#' @param x a character string or vector of character strings. Note, these are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @export

position_from_variant_string <- function(x) {

  with_unique(x, function(x_unique) {
    vapply(strsplit(x_unique, ";"), function(genes) {

      # keep only the gene name and codon position of each gene (the first two
      # colon-separated elements), dropping amino acids and read counts
      pos_only <- sub("^([^:]*:[^:]*):.*$", "\\1", genes)

      # genes are sorted alphabetically before any manipulation
      gene_names <- sub(":.*$", "", pos_only)
      paste(pos_only[order(gene_names, method = "radix")], collapse = ";")
    }, character(1))
  })
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

  with_unique(variant_strings, function(v_unique) {
    vapply(variant_to_long(v_unique), function(df) {
      keep <- (df$gene %in% df_position$gene) & (df$pos %in% df_position$pos)
      if (!any(keep)) {
        return(NA_character_)
      }
      long_to_variant(list(df[keep, , drop = FALSE]))
    }, character(1))
  })
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

  with_unique(x, function(x_unique) {
    ordered <- lapply(variant_to_long(x_unique), function(y) {
      y[order(y$gene, y$pos, y$aa, method = "radix"), , drop = FALSE]
    })
    long_to_variant(ordered)
  })
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

  with_unique(x, function(x_unique) {
    ordered <- lapply(position_to_long(x_unique), function(y) {
      y[order(y$gene, y$pos, method = "radix"), , drop = FALSE]
    })
    long_to_position(ordered)
  })
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
  count_hets(x, phased = FALSE)
}

# internal: count heterozygous loci of a given phasing in each variant string
count_hets <- function(x, phased) {
  with_unique(x, function(x_unique) {
    vapply(variant_to_long(x_unique), function(y1) {
      # one row per locus, then count hets of the requested phasing
      first <- !duplicated(paste(y1$gene, y1$pos, sep = "\r"))
      sum(y1$het[first] & (y1$phased[first] == phased))
    }, numeric(1))
  })
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
  count_hets(x, phased = TRUE)
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

  with_unique(x, function(x_unique) {
    vapply(strsplit(x_unique, ";"), function(genes) {

      # keep only the gene name, codon position and amino acid (the first three
      # colon-separated elements), dropping any read counts
      no_reads <- sub("^(([^:]*:){2}[^:]*):.*$", "\\1", genes)

      # genes are sorted alphabetically before any manipulation
      gene_names <- sub(":.*$", "", no_reads)
      paste(no_reads[order(gene_names, method = "radix")], collapse = ";")
    }, character(1))
  })
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

  # parse the target once. A heterozygous locus produces more than one row at a
  # given gene/position, so duplicated loci flag a heterozygous target
  df_target <- variant_to_long(target_string)[[1]]
  if (any(duplicated(paste(df_target$gene, df_target$pos, sep = "\r")))) {
    stop("target string cannot contain any heterozygous loci")
  }
  target_key <- paste(df_target$gene, df_target$pos, df_target$aa, sep = "\r")

  # work on unique comparison strings (parsed in one pass) and map back at end
  c_unique <- unique(comparison_strings)
  comparison_list <- variant_to_long(c_unique)

  df_unique <- data.frame(match = logical(length(c_unique)),
                          ambiguous = logical(length(c_unique)),
                          prop = NA_real_)

  for (i in seq_along(c_unique)) {

    dfc <- comparison_list[[i]]
    comp_key <- paste(dfc$gene, dfc$pos, dfc$aa, sep = "\r")

    # the target is only present if every one of its alleles is found in the
    # comparison
    if (!all(target_key %in% comp_key)) {
      df_unique$match[i] <- FALSE
      df_unique$prop[i] <- 0
      next
    }
    df_unique$match[i] <- TRUE

    # flag comparison alleles that correspond to the target, then aggregate to
    # the matching loci
    matched <- comp_key %in% target_key
    locus <- paste(dfc$gene, dfc$pos, sep = "\r")
    idx <- split(seq_along(locus), factor(locus, levels = unique(locus)))
    keep <- vapply(idx, function(ii) any(matched[ii]), logical(1))
    idx <- idx[keep]
    het_l <- vapply(idx, function(ii) dfc$het[ii][1], logical(1))
    phased_l <- vapply(idx, function(ii) dfc$phased[ii][1], logical(1))
    rc_total <- vapply(idx, function(ii) sum(dfc$read_count[ii]), numeric(1))
    rc_match <- vapply(idx, function(ii) sum(dfc$read_count[ii] * matched[ii]), numeric(1))

    n_phased_het <- sum(het_l & phased_l)
    n_unphased_het <- sum(het_l & !phased_l)

    # a match is ambiguous if it spans more than one heterozygous locus
    ambiguous <- ((n_phased_het > 0) + n_unphased_het) > 1
    df_unique$ambiguous[i] <- ambiguous

    # within-sample genotype frequency, where it can be determined
    if (!ambiguous) {
      if (any(is.na(rc_match))) {
        if ((n_phased_het + n_unphased_het) == 0) {
          df_unique$prop[i] <- 1
        }
      } else {
        prop <- rc_match / rc_total
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
  df_unique[match(comparison_strings, c_unique), ]
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

  # get target as a set of gene/position keys
  df_target <- position_to_long(target_string)[[1]]
  target_key <- paste(df_target$gene, df_target$pos, sep = "\r")

  with_unique(comparison_strings, function(c_unique) {

    # get comparison positions in long form (one parse for all)
    list_comparison <- position_to_long(position_from_variant_string(c_unique))

    # a match requires every target position to be present in the comparison
    vapply(list_comparison, function(dfc) {
      all(target_key %in% paste(dfc$gene, dfc$pos, sep = "\r"))
    }, logical(1))
  })
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

  with_unique(x, function(x_unique) {
    lapply(variant_to_long(x_unique), function(y) {
      sprintf("%s:%s:%s", y$gene, y$pos, y$aa)
    })
  })
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
#' @export

get_component_variants <- function(x) {

  with_unique(x, function(x_unique) {

    # parse once; read counts are irrelevant to component genotypes
    lapply(variant_to_long(x_unique), function(df) {
      df$read_count <- NA

      # collapse to one row per locus (in sorted gene, position order) and the
      # list of alleles observed at each
      key <- paste(df$gene, df$pos, sep = "\r")
      ord <- order(df$gene[!duplicated(key)], df$pos[!duplicated(key)], method = "radix")
      loci <- data.frame(gene = df$gene[!duplicated(key)],
                         pos = df$pos[!duplicated(key)],
                         het = df$het[!duplicated(key)],
                         phased = df$phased[!duplicated(key)])[ord, , drop = FALSE]
      allele_lists <- split(df$aa, factor(key, levels = unique(key)))[ord]

      n_unphased <- sum(loci$het & !loci$phased)
      n_phased <- sum(loci$het & loci$phased)

      # no heterozygous loci: the genotype is already fully resolved
      if (n_unphased == 0 && n_phased == 0) {
        return(long_to_variant(list(df)))
      }

      # build a matrix of allele combinations (one column per locus). With a
      # single unphased het we take the cartesian product; with phased hets each
      # combination follows the phasing, and homozygous loci are recycled.
      if (n_unphased == 1 && n_phased == 0) {
        grid <- expand.grid(allele_lists, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
      } else if (n_unphased == 0 && n_phased > 0) {
        k <- max(lengths(allele_lists))
        grid <- as.data.frame(lapply(allele_lists, function(a) if (length(a) == 1) rep(a, k) else a),
                              stringsAsFactors = FALSE)
      } else {
        # two or more heterozygous loci: component genotypes are ambiguous
        return(NA)
      }

      # turn each combination into a fully resolved variant string
      components <- lapply(seq_len(nrow(grid)), function(r) {
        data.frame(gene = loci$gene, pos = loci$pos, n_aa = 1, het = FALSE, phased = FALSE,
                   aa = unlist(grid[r, , drop = FALSE], use.names = FALSE), read_count = NA)
      })
      long_to_variant(components)
    })
  })
}

#------------------------------------------------
#' @title Overlay one variant string onto another
#'
#' @description
#' All the positions and amino-acids from the first variant are combined with
#' the second, overwriting where positions are the same.
#'
#' @details
#' Note, inputs are not internally checked for being valid variant strings, it is up to the user to ensure this (see \code{?check_variant_string}).
#'
#' @param var1 the first variant string, which overwrites in the case of shared positions.
#' @param var2 the second variant string, which is overwritten in the case of shared positions.
#'
#' @import dplyr
#'
#' @export

overlay_variant <- function(var1, var2) {

  mapply(function(x, y) {
    if (is.na(x) & is.na(y)) {
      return(NA)
    }
    if (is.na(x)) {
      return(y)
    }
    if (is.na(y)) {
      return(x)
    }
    full_join(variant_to_long(x)[[1]], variant_to_long(y)[[1]], by = join_by(gene, pos)) |>
      mutate(n_aa = coalesce(n_aa.x, n_aa.y),
             het = coalesce(het.x, het.y),
             phased = coalesce(phased.x, phased.y),
             aa = coalesce(aa.x, aa.y),
             read_count = coalesce(read_count.x, read_count.y)) |>
      select(gene, pos, n_aa, het, phased, aa, read_count) |>
      list() |>
      long_to_variant()
  }, var1, var2)

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
