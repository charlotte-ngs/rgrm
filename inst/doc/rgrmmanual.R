## ----LoadPackageStuff, echo=FALSE, results='hide'------------------------
if (!require(rgrm))
  devtools::load_all()

## ----GenerateMatrix, echo=FALSE------------------------------------------
set.seed(2143)
nNrInd <- 6
mGrmP <- matrix(runif(n = nNrInd * nNrInd), ncol = nNrInd)
mGrm <- crossprod(mGrmP)/nNrInd
diag(mGrm) <- 1+diag(mGrm)
print(mGrm)

## ----FlattenedGrm, echo=FALSE--------------------------------------------
vFlatMatElem <- rgrm::flattenMatrixLowerTri(pMatrix = mGrm)
print(vFlatMatElem)

## ----InterestingElements, echo=FALSE-------------------------------------
vInterestingElements <- c(8,13)

## ----ComputeElements-----------------------------------------------------
vInterestingElements <- c(8,13)
lSearchResult <- rgrm::findMatrixElements(vInterestingElements, vFlatMatElem) 
print(lSearchResult)

## ----ReadGrmBinToVector, eval=TRUE---------------------------------------
sInputPath <- system.file(package = "rgrm", file.path("extdata", "OmicKriging"))
sInputPrefix <- "ig_genotypes.grm"
lFlatMatElem <- rgrm::readGRMBinToVector(psInputPrefix = sInputPrefix, psInputPath = sInputPath)
vFlatMatElem <- lFlatMatElem$flatVec
length(vFlatMatElem)
head(vFlatMatElem)
nrow(lFlatMatElem$id)
head(lFlatMatElem$id)

## ----SmallSearchWithIds--------------------------------------------------
rgrm::findMatrixElements(1:3, vFlatMatElem, as.vector(lFlatMatElem$id[,2]))

## ------------------------------------------------------------------------
sessionInfo()

