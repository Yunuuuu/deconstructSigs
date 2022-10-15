#' Converts mutation list to correct input format
#'
#' Given a mutation list, outputs a data frame with counts of how frequently a
#' mutation is found within each trinucleotide context per sample ID.  Output
#' can be used as input into getTriContextFraction.
#'
#' The context sequence is taken from the BSgenome.Hsapiens.UCSC.hgX::Hsapiens
#' object. Therefore the coordinates must correspond to the human hgX assembly.
#' Default is set to the UCSC hg19 assembly, which corresponds to the GRCh37
#' assembly. If another assembly is required, it must already be present in the
#' R workspace and fed as a parameter. This method will translate chromosome
#' names from other versions of the assembly like NCBI or Ensembl. For instance,
#' the following transformation will be done: "1" -> "chr1"; "MT" -> "chrM";
#' "GL000245.1" -> "chrUn_gl000245"; etc.
#'
#' @param mut.ref Location of the mutation file that is to be converted or name
#'   of data frame in environment
#' @param sample.id Column name in the mutation file corresponding to the Sample
#'   ID
#' @param chr Column name in the mutation file corresponding to the chromosome
#' @param pos Column name in the mutation file corresponding to the mutation
#'   position
#' @param ref Column name in the mutation file corresponding to the reference
#'   base
#' @param alt Column name in the mutation file corresponding to the alternate
#'   base
#' @param bsg Only set if another genome build is required. Must be a BSgenome
#'   object.
#' @param chr.list what targetedd chromosome should be used in the analysis.
#' Default: NULL, means all chromosome
#' @param sig.type Are SBS or DBS signatures being used?
#' @return A data frame that contains sample IDs for the rows and trinucleotide
#'   contexts for the columns. Each entry is the count of how many times a
#'   mutation with that trinucleotide context is seen in the sample.
#' @examples
#' \dontrun{
#' sample.mut.ref <- readRDS(
#'   system.file("extdata", "sample.mut.ref.rds",
#'     package = "deconstructSigs"
#'   )
#' )
#' sigs.input <- mut.to.sigs.input(
#'   mut.ref = sample.mut.ref, sample.id = "Sample",
#'   chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg =
#'     BSgenome.Hsapiens.UCSC.hg19, sig.type = "SBS"
#' )
#' }
#' @export
mut.to.sigs.input <- function(mut.ref, sample.id = "Sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = NULL, chr.list = NULL, sig.type = "SBS") {
  sig.type <- match.arg(sig.type, c("SBS", "DBS"))

  # prepare data -------------------------------
  if (exists("mut.ref", mode = "list")) {
    mut.full <- mut.ref
  } else {
    if (file.exists(mut.ref)) {
      mut.full <- utils::read.table(mut.ref, sep = "\t", header = TRUE, as.is = FALSE, check.names = FALSE)
    } else {
      stop("mut.ref is neither a file nor a loaded data frame")
    }
  }

  mut <- mut.full[, c(sample.id, chr, pos, ref, alt)]
  data.table::setDT(mut)
  data.table::setnames(mut, c("sample.id", "chr", "pos", "ref", "alt"))

  if (!is.null(chr.list)) {
    mut <- mut[chr %in% chr.list]
  }

  # parse DBS -----------------------------------------
  # don't need trinucleotide context if looking for DBS
  if (identical(sig.type, "DBS")) {
    mut[, c("ref", "alt") := lapply(.SD, as.character),
      .SDcols = c("ref", "alt")
    ]
    mut <- mut[nchar(ref) == 2L & nchar(alt) == 2L]
    mut[
      , dbs_condensed := dbs_possible$dbs_condensed[
        match(paste(ref, alt, sep = ">"), dbs_possible$dbs)
      ]
    ]
    mut[, dbs_condensed := factor(dbs_condensed, unique(dbs_condensed))]
    final.df <- as.data.frame(
      table(
        samples = as.character(mut[["sample.id"]]),
        features = mut[["dbs_condensed"]]
      ),
      stringsAsFactors = FALSE,
      responseName = "counts"
    )
    data.table::setDT(final.df)
    final.df <- data.table::dcast(
      final.df,
      samples ~ features,
      drop = FALSE,
      value.var = "counts"
    )
    final.df <- column_to_rownames(final.df, "samples")
  }

  # parse SBS --------------------------------------------
  # mut.lengths <- with(mut, nchar(as.character(mut[,ref])))
  # mut.lengths <- with(mut, nchar(as.character(ref)))
  # mut <- mut[which(mut.lengths == 1),]
  # mut$mut.lengths <- nchar(as.character(mut[, ref]))
  if (identical(sig.type, "SBS")) {
    mut <- mut[
      ref %in% c("A", "T", "C", "G") &
        alt %in% c("A", "T", "C", "G") &
        !is.na(chr)
    ]

    # Fix the chromosome names (in case they come from Ensembl instead of UCSC)
    mut[, chr := GenomeInfoDb::mapSeqlevels(
      chr,
      style = "UCSC",
      best.only = TRUE, drop = TRUE
    )]
    mut[, chr := factor(chr)]
    # mut[, chr] <- factor(mut[, chr])
    # levels(mut[, chr]) <- sub("^([0-9XY])", "chr\\1", levels(mut[, chr]))
    # levels(mut[, chr]) <- sub("^MT", "chrM", levels(mut[, chr]))
    # levels(mut[, chr]) <- sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1", levels(mut[, chr]), perl = TRUE)

    # Check the genome version the user wants to use
    # If set to default, carry on happily
    if (is.null(bsg)) {
      bsg <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
    } else if (!inherits(bsg, "BSgenome")) {
      stop("The bsg parameter needs to either be set to NULL or a BSgenome object.")
    }
    # Remove any entry in chromosomes that do not exist in the BSgenome.Hsapiens.UCSC.hgX::Hsapiens object
    unknown.regions <- levels(mut[["chr"]])[
      which(!(levels(mut[["chr"]]) %in% GenomeInfoDb::seqnames(bsg)))
    ]
    if (length(unknown.regions) > 0L) {
      unknown.regions <- paste(unknown.regions, collapse = ",\ ")
      warning(
        paste("Check chr names -- not all match", attr(bsg, which = "pkgname"),
          "object:\n", unknown.regions,
          sep = " "
        )
      )
      mut <- mut[chr %in% GenomeInfoDb::seqnames(bsg)]
    }

    # Add in context
    context <- BSgenome::getSeq(
      bsg, mut[["chr"]],
      mut[["pos"]] - 1L,
      mut[["pos"]] + 1L,
      as.character = FALSE
    )
    mutcat <- Biostrings::DNAStringSet(paste0(mut[["ref"]], mut[["alt"]]))

    conflicted_idx <- Biostrings::DNAStringSet(mut[["ref"]]) !=
      Biostrings::subseq(context, 2L, 2L)
    if (any(conflicted_idx)) {
      bad_mut <- mut[
        conflicted_idx,
        paste0(sample.id, ": ", chr, " ", pos, " ",
          ref, " ", alt,
          collapse = "\n"
        )
      ]
      warning(paste("Check ref bases -- not all match context:\n ",
        bad_mut,
        sep = " "
      ))
    }

    # Reverse complement the G's and A's
    rev_idx <- mut[["ref"]] %in% c("A", "G")

    std.mutcat <- mutcat
    std.mutcat[rev_idx] <- Biostrings::complement(std.mutcat)[rev_idx]
    std.mutcat <- sub("([ACTG]$)", ">\\1", as.character(std.mutcat))

    std.context <- context
    std.context[rev_idx] <- Biostrings::reverseComplement(
      std.context
    )[rev_idx]

    # Make the tricontext -----------------------------------------
    tricontext <- unname(as.character(std.context))
    tricontext <- paste0(
      substr(tricontext, 1, 1),
      "[", std.mutcat, "]",
      substr(tricontext, 3, 3)
    )
    tricontext <- factor(tricontext, levels = generate_sbs_features())

    # Generate all possible trinucleotide contexts
    final.df <- as.data.frame(
      table(
        samples = as.character(mut[["sample.id"]]),
        features = tricontext
      ),
      stringsAsFactors = FALSE,
      responseName = "counts"
    )
    data.table::setDT(final.df)
    final.df <- data.table::dcast(
      final.df,
      samples ~ features,
      drop = FALSE,
      value.var = "counts"
    )
    final.df <- column_to_rownames(final.df, "samples")

    # previous implementation  -----------------------------
    # mut$std.mutcat <- paste0(mut[["ref"]], ">", mut[["alt"]])
    # mut$std.mutcat[rev_idx] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.mutcat[rev_idx])))) # to lowercase
    # mut$std.mutcat[rev_idx] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.mutcat[rev_idx])))) # complement
    #
    # mut$std.context <- as.character(context)
    # mut$std.context[rev_idx] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.context[rev_idx])))) # to lowercase
    # mut$std.context[rev_idx] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.context[rev_idx])))) # complement
    # mut$std.context[rev_idx] <- sapply(strsplit(mut$std.context[rev_idx], split = ""), function(str) {
    #   paste(rev(str), collapse = "")
    # }) # reverse
    # all(mut$std.context == as.character(std.context)) # [1] TRUE
    # all(mut$std.mutcat == std.mutcat) # [1] TRUE
    #
    # # Make the tricontext
    # mut$tricontext <- paste(
    #   substr(mut$std.context, 1, 1),
    #   "[", mut$std.mutcat, "]",
    #   substr(mut$std.context, 3, 3),
    #   sep = ""
    # )
    # all(mut$tricontext == tricontext) # [1] TRUE
    #
    # # Generate all possible trinucleotide contexts
    # all.tri <- generate_sbs_features()
    # samples <- as.character(unique(mut[["sample.id"]]))
    # final.matrix <- matrix(
    #     0,
    #     ncol = 96L,
    #     nrow = length(samples)
    # )
    # colnames(final.matrix) <- all.tri
    # rownames(final.matrix) <- samples
    # # print(paste("[", date(), "]", "Fill in the context matrix"))
    # for (i in samples) {
    #     tmp <- mut[which(sample.id == i), ]
    #     beep <- table(tmp$tricontext)
    #     for (l in 1:length(beep)) {
    #         trimer <- names(beep[l])
    #         if (trimer %in% all.tri) {
    #             final.matrix[i, trimer] <- beep[trimer]
    #         }
    #     }
    # }
    #
    # final.df2 <- data.frame(final.matrix, check.names = FALSE)
    # browser()
    # setequal(colnames(final.df), colnames(final.df2))
    # final.df <- final.df[, colnames(final.df2)]
    # all(colnames(final.df) == colnames(final.df2))
    # all(final.df == final.df2) # [1] TRUE
  }

  bad_samples <- rownames(final.df)[rowSums(final.df) <= 50]
  if (length(bad_samples) > 0L) {
    bad_samples <- paste(bad_samples, collapse = ", ")
    warning(paste("Some samples have fewer than 50 mutations:\n ", bad_samples, sep = " "))
  }
  final.df
}

generate_sbs_features <- function() {
  mut <- expand.grid(
    x = c("C", "T"),
    y = c("A", "C", "G", "T"),
    stringsAsFactors = FALSE
  )
  mut <- mut[mut$x != mut$y, ]
  mut <- .mapply(paste, mut, list(sep = ">"))
  mut <- unlist(mut, use.names = FALSE, recursive = FALSE)
  sbs_feature <- expand.grid(
    left = c("A", "C", "G", "T"),
    mut = paste0("[", mut, "]"),
    right = c("A", "C", "G", "T"),
    stringsAsFactors = FALSE
  )
  sbs_feature <- sbs_feature[order(sbs_feature$mut), ]
  unlist(.mapply(paste0, sbs_feature, NULL),
    use.names = FALSE, recursive = FALSE
  )
}

# all.tri <- c()
# for (i in c("A", "C", "G", "T")) {
#   for (j in c("C", "T")) {
#     for (k in c("A", "C", "G", "T")) {
#       if (j != k) {
#         for (l in c("A", "C", "G", "T")) {
#           tmp <- paste(i, "[", j, ">", k, "]", l, sep = "")
#           all.tri <- c(all.tri, tmp)
#         }
#       }
#     }
#   }
# }
# all.tri <- all.tri[order(substr(all.tri, 3, 5))]
# setequal(all.tri, generate_sbs()) # [1] TRUE
