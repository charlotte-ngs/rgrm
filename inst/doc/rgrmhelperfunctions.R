## ----GenerateGRM---------------------------------------------------------
set.seed(4321)
nNrInd <- 6
mGrmP <- matrix(runif(n = nNrInd * nNrInd), ncol = nNrInd)
mGrm <- crossprod(mGrmP)/nNrInd
diag(mGrm) <- 1+diag(mGrm)
print(mGrm)

## ----FlattenGrm----------------------------------------------------------
vGrmLowTri <- vector(mode = "numeric", length = sum(1:nNrInd))
nVecIdx <- 1
for (nRowIdx in 1:nrow(mGrm)){
  for (nColIdx in 1:nRowIdx){
    vGrmLowTri[nVecIdx] <- mGrm[nRowIdx, nColIdx]
    nVecIdx <- nVecIdx + 1
  }
}
print(vGrmLowTri)

## ----LowerTriMatrix------------------------------------------------------
lower.tri(mGrm, diag = TRUE)

## ----FlattenRowWise------------------------------------------------------
mGrm[lower.tri(mGrm, diag = TRUE)]

## ----FlattenColumnWise---------------------------------------------------
vFlatMatElem <- mGrm[upper.tri(mGrm, diag = TRUE)]
print(vFlatMatElem)
identical(vFlatMatElem, vGrmLowTri)

## ----DiagonalElementsIndex-----------------------------------------------
vFlatDiagIdx <- cumsum(1:nrow(mGrm))
print(vFlatDiagIdx)

## ----ExtractDiagonalElement----------------------------------------------
vFlatMatElem[vFlatDiagIdx]

## ----ExtractOffDiagonals-------------------------------------------------
vFlatMatElem[-c(vFlatDiagIdx)]

## ----AssignElem12, echo=FALSE--------------------------------------------
nElemIdx <- 12

## ----Element12-----------------------------------------------------------
nElemIdx <- 12
vFlatMatElem[nElemIdx]

## ----ElementRow, echo=FALSE----------------------------------------------
nElementRow <- which(vFlatDiagIdx > nElemIdx)[1]

## ----GetMatRowIndex------------------------------------------------------
nElemIdx <- 12
nElementRow <- which(vFlatDiagIdx > nElemIdx)[1]
nElementCol <- nElemIdx-vFlatDiagIdx[nElementRow-1]
cat(" * Row: ", nElementRow, "\n * Col: ", nElementCol, "\n")

## ----CheckResult---------------------------------------------------------
identical(vFlatMatElem[12], mGrm[nElementRow,nElementCol])

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

## ----SearchElement12-----------------------------------------------------
lElemFound <- findSingleMatElement(pnElementIdx = 12, pvFlatVec = vFlatMatElem)
print(lElemFound)
identical(lElemFound$nCoefficient, mGrm[lElemFound$nRowIndex, lElemFound$nColIndex])

## ----ElementThreshold, echo=FALSE----------------------------------------
nElementThreshold <- 0.12

## ----SelElemGroup--------------------------------------------------------
nElementThreshold <- 0.12
vElemGroup <- which(vFlatMatElem < nElementThreshold)

## ----SearchElementGroup--------------------------------------------------
lapply(vElemGroup, findSingleMatElement, vFlatMatElem)

## ----CompleteTest--------------------------------------------------------
lCompTestResults <- lapply(1:length(vFlatMatElem), findSingleMatElement, vFlatMatElem)

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
    nMatElement <- mGrm[lCurResult$nRowIndex,lCurResult$nColIndex]
    bCompareResult <- identical(lCurResult$nCoefficient,mGrm[lCurResult$nRowIndex,lCurResult$nColIndex])
  } else {
    nRowIndex <- c(nRowIndex, lCurResult$nRowIndex)
    nColIndex <- c(nColIndex, lCurResult$nColIndex)
    nCoefficient <- c(nCoefficient, lCurResult$nCoefficient)
    nMatElement <- c(nMatElement, mGrm[lCurResult$nRowIndex,lCurResult$nColIndex])
    bCompareResult <- c(bCompareResult, identical(lCurResult$nCoefficient,mGrm[lCurResult$nRowIndex,lCurResult$nColIndex]))
    
  }
}
dfResults <- data.frame(nRowIndex, 
                        nColIndex, 
                        nCoefficient, 
                        nMatElement, 
                        bCompareResult, 
                        stringsAsFactors = FALSE)

## ----ResultTable---------------------------------------------------------
knitr::kable(dfResults)

## ----CheckForDiff--------------------------------------------------------
any(bCompareResult == FALSE)

## ------------------------------------------------------------------------
sessionInfo()

