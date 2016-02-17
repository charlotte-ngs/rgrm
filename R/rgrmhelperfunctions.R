## ----Flatten Matrix------------------------------------------------------
#' Flatten lower triangular part of a matrix into a vector
#'
#' \code{flattenMatrixLowerTri} uses function \code{lower.tri()} and
#' \code{upper.tri()} to flatten the lower triangular part of a matrix
#' into a vector.
#'
#' @param pMatrix           Matrix to be flattened
#' @param psOrder           Order in which elements should be stored in the flattened vector
#' @param pbIncludeDiag     whether or not to include diagonal elements of pMatrix
#' @return vFlatVecResult   flattened vector
#' @export
flattenMatrixLowerTri <- function(pMatrix, psOrder = "rowwise", pbIncludeDiag = TRUE){
  if(psOrder != "rowwise" & psOrder != "columnwise")
    stop("Unknown flattening order: ", psOrder)
  if (psOrder == "columnwise") {
    vFlatVecResult <- pMatrix[lower.tri(pMatrix, diag = pbIncludeDiag)]
  } else {
    vFlatVecResult <- pMatrix[upper.tri(pMatrix, diag = pbIncludeDiag)]
  }
  return(vFlatVecResult)
}


## ----FindGrmElement------------------------------------------------------
#' Find a single element of a symmetric matrix based on a flattened vector representation
#'
#' Given a flattened vector representation of a symmatric matrix \code{findSingleMatElement}
#' computes the row and column indices in the original matrix of a given element. The element
#' is given by the index of the flattened vector. The flattened vector is a vector that contains
#' all elements of the lower triangular part of the matrix including all diagonal elements. The
#' order in which the elements are written to the vector is going through each row and copying
#' each element from column one up to the diagonal element of that row into the flattened
#' vector.
#'
#' @param  pnElementIdx   Vector index of element to be searched for
#' @param  pvFlatVec      Flattened vector represenation of symmetric matrix
#' @return lResultList    Result list with components nRowIndex, nColIndex and nCoefficient
#' @export
findSingleMatElement <- function(pnElementIdx, pvFlatVec){
  # pnElementIdx must a be positive integer
  if (pnElementIdx < 0)
    stop(paste("Length of flattened vector must be positive, found: ", pnElementIdx))
  if (as.integer(pnElementIdx) - pnElementIdx != 0)
    stop(paste("Length of flattened vector must be an integer, found: ", pnElementIdx))
  ### # pnElementIdx cannot be larger than the length of pvFlatVec
  if (pnElementIdx > length(pvFlatVec))
    stop("Elementindex cannot be larger than length of index vector")
  ### # compute the dimensionality of the original matrx
  nMatDim <- nGetMatDim(pnLenFlatVec = length(pvFlatVec))
  ### # compute vector of indices of diagonal elements
  vDiagIdx <- cumsum(1:nMatDim)
  ### # check whether pnElementIdx is a diagonal element, if yes, then special case
  if (any(vDiagIdx == pnElementIdx)) {
    nDiagIdx <- which(vDiagIdx == pnElementIdx)
    lResultList <- list(nRowIndex=nDiagIdx, nColIndex=nDiagIdx, nCoefficient = pvFlatVec[pnElementIdx])
  } else {
    ### # off-diagonal
    nElementRow <- which(vDiagIdx > pnElementIdx)[1]
    nElementCol <- pnElementIdx - vDiagIdx[nElementRow-1]
    lResultList <- list(nRowIndex=nElementRow, nColIndex=nElementCol, nCoefficient = pvFlatVec[pnElementIdx])
  }
  return(lResultList)
}

## ----GetMatDim-----------------------------------------------------------
#' Compute matrix dimension based on the length of a flattened vector
#'
#' The flattened vector contains all elements of the lower triangular
#' part of a symmetric matrix including all diagonal elements. Because
#' any row i of a quadratic matrix contains i elements below the
#' diagonal, including the diagonal element, the flattened vector of
#' an n * n symmetric matrix contains the sum of all natural numbers
#' up to and including n. That length has a closed form and corresponds
#' to n*(n+1)/2. Given the length of the flattened vector, we have
#' a quadratic equation for n, that can be solved. We are only interested
#' in the positive solution which is n = (-1 + \sqrt(1+8l))/2
#'
#' @param    pnLenFlatVec    length of flattened vector
#' @return   nMatDimResult   dimension of symmetric matrix
#' @export
nGetMatDim <- function(pnLenFlatVec) {
  # pnLenFlatVec must a be positive integer
  if (pnLenFlatVec < 0)
    stop(paste("Length of flattened vector must be positive, found: ", pnLenFlatVec))
  if (as.integer(pnLenFlatVec) - pnLenFlatVec != 0)
    stop(paste("Length of flattened vector must be an integer, found: ", pnLenFlatVec))
  # compute the matrix dimension
  nMatDimResult <- (sqrt(1 + 8*pnLenFlatVec) - 1) / 2
  # result must be integer
  if ((as.integer(nMatDimResult) - nMatDimResult) != 0 )
    stop(paste("Resulting matrix dimension not an integer:", nMatDimResult ))
  return(nMatDimResult)
}

## ----Finding groups of elements------------------------------------------
#' Find row and column indices and matrix values for a vector of elements
#'
#' @description
#' Given a vector of element indices of a flattened vector and the flattened
#' vector itself, find row and column indices and coefficients (matrix values)
#' for the specified elements.
#'
#' @details
#' Flattening a symmetric matrix corresponds to saving the elements of the lower
#' (or upper) triangular part of the matrix in a defined order in a vector. Given
#' such a flattened representation of a matrix one might search for the exact
#' location where in the original matrix certain elements did occur.
#'
#' @param  pvElements     Vector of indices specifying elements to be searched for
#' @param  pvFlatVec      Flattened vector represenation of symmetric matrix
#' @return lResultLoL     Result list of lists with components nRowIndex, nColIndex and nCoefficient
#' @export
findMatrixElements <- function(pvElements, pvFlatVec){
  lResultLoL <- lapply(pvElements, findSingleMatElement, pvFlatVec)
  return(lResultLoL)
}

## ----CompleteTest--------------------------------------------------------
#' Test all elements of flattened vector
#'
#' The test consists of searching for all elements in the flattened
#' vector and comparing the found elements with the original elements
#' in the matrix
#'
#' @param    pvFlatMatElem   flattened vector
#' @param    pMatrix         original matrix
#' @return   dfResults       dataframe with comparison results
#' @export
runCompleteTest <- function(pvFlatMatElem, pMatrix) {

  ### # search all elements of the flattened vector
  lCompTestResults <- lapply(1:length(pvFlatMatElem), findSingleMatElement, pvFlatMatElem)

  ## ----CompareResults------------------------------------------------------
  nRowIndex <- NULL
  nColIndex <- NULL
  nCoefficient <- NULL
  nMatElement <- NULL
  bCompareResult <- NULL
  for (lCurResult in lCompTestResults){
    ### # collect items
    if (is.null(nRowIndex)){
      nRowIndex <- lCurResult$nRowIndex
      nColIndex <- lCurResult$nColIndex
      nCoefficient <- lCurResult$nCoefficient
      nMatElement <- pMatrix[lCurResult$nRowIndex,lCurResult$nColIndex]
      bCompareResult <- identical(lCurResult$nCoefficient,pMatrix[lCurResult$nRowIndex,lCurResult$nColIndex])
    } else {
      nRowIndex <- c(nRowIndex, lCurResult$nRowIndex)
      nColIndex <- c(nColIndex, lCurResult$nColIndex)
      nCoefficient <- c(nCoefficient, lCurResult$nCoefficient)
      nMatElement <- c(nMatElement, pMatrix[lCurResult$nRowIndex,lCurResult$nColIndex])
      bCompareResult <- c(bCompareResult, identical(lCurResult$nCoefficient,pMatrix[lCurResult$nRowIndex,lCurResult$nColIndex]))

    }
  }
  dfResults <- data.frame(nRowIndex,
                          nColIndex,
                          nCoefficient,
                          nMatElement,
                          bCompareResult,
                          stringsAsFactors = FALSE)

  return(dfResults)
}
