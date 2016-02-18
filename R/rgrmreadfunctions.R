###
###
###
###   Purpose:   Reading GRM data
###   started:   2016/02/18 (pvr)
###
### ############################### ###

# ---- R script to read the GRM binary file ----
#' Reading Genetic Relationship Matrix from binary files
#'
#' @description
#' The function \code{ReadGRMBin} was taken from
#' http://cnsgenomics.com/software/gcta/estimate_grm.html
#' and was only slightly modified, but in essence
#' the function here does the same thing. Given three
#' files which have names starting with the same prefix,
#' a list with the diagonal elements, the offdiagonal
#' elements and the ids is returned
#'
#' @param    prefix   common prefix of input files
#' @param    AllN     flag indicating whether to read 1 line (FALSE) or n(n+1)/2 lines from file prefix.N.bin
#' @param    size     wordsize that is passed to readBin
#' @return   list with diagonal elements, offdiagonals, ids and total number of individuals
#' @export ReadGRMBin
ReadGRMBin=function(prefix, AllN=F, size=4){
  ### # the following can be replaced by cumsum
  #sum_i=function(i){
  #  return(sum(1:i))
  #}
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <-  read.table(IDFileName)
  #n=dim(id)[1]
  n <- nrow(id)
  BinFile <- file(BinFileName, "rb");
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb");
  if (AllN == T){
    N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  } else {
    N <- readBin(NFile, n=1, what=numeric(0), size=size)
  }
  #i=sapply(1:n, sum_i)
  i <- cumsum(1:n)
  ### # return list of results
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}


# --- Read binary GRM data into a vector --------------------
#' Read a GRM into a Vector
#'
#' @description
#' Given two input files
#'
#' \itemize{
#' \item a file with the flattened GRM in binary format
#' \item a file with indices
#' }
#'
#' The data from the binary file are read using readLines and
#' are stored in a vector
#'
#' @param     prefix   common prefix of input files
#' @param     size     word size passed to readBin
#' @return    vector with flattened represenation of GRM
#' @export readGRMBinToVector
readGRMBinToVector <- function(prefix, size = 4){
  ### # complete names of input files
  BinFileName <- paste0(prefix, ".bin")
  IDFileName <- paste0(prefix, ".id")
  ### # read dataframe with ids
  dfIds <- read.table(IDFileName)
  nNrRowIds <- nrow(dfIds)
  vFlatGrm <- readBin(con = BinFileName,
                      n = as.integer(nNrRowIds * (nNrRowIds + 1)/2),
                      what = numeric(0L),
                      size = size)
  close(BinFileName)
  return(vFlatGrm)

}
