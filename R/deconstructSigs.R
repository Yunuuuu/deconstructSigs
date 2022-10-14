#' deconstructSigs
#'
#' Takes sample information in the form of the fraction of mutations
#' in each of 96 trinucleotide contexts and identifies the weighted combination
#' of published signatures that, when summed, most closely reconstructs the
#' mutational profile.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{whichSignatures}}
#'   \item \code{\link{mut.to.sigs.input}}
#'   \item \code{\link{plotSignatures}}
#'   \item \code{\link{makePie}}
#'   }
#'
#' @docType package
#' @author Rachel Rosenthal rachel.rosenthal.14@ucl.ac.uk
#' @name deconstructSigs
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table :=
#' @importFrom data.table data.table
## usethis namespace: end
NULL

utils::globalVariables(
    c("mutation", "full_context", "mutation")
)
