parse_sig <- function(x, transpose = TRUE) {
  if (is.character(x)) {
    data <- data.table::fread(x, header = TRUE)
  } else if (inherits(x, "data.frame")) {
    data <- x
  } else {
    stop("x should be a data.frame or a file path (can be also a URL starting http://, file://, etc) read as a data.frame")
  }
  sig.df <- as.data.frame(
    data[, -1L],
    row.names = as.character(data[[1L]]),
    make.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (transpose) {
    sig.df <- as.data.frame(
      t(as.matrix(sig.df)),
      make.names = FALSE, 
      stringsAsFactors = FALSE
    )
  } 
  sig.df
}

check_sig <- function(signatures.ref, sig.type = "SBS") {
  sig_features <- switch (sig.type,
    SBS = generate_sbs_features(),
    DBS = unique(dbs_possible$dbs_condensed)
  )
  is_compatible <- setequal(
    colnames(signatures.ref),
    sig_features
  )
  if (!is_compatible) {
    stop(paste0("signatures.ref don't contain all ", length(sig_features), " signature features"))
  }
}

