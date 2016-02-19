###
###
###
###   Purpose:   Prepare test data for reading functions
###   started:   2016/02/18 (pvr)
###
### ##################################################### ###

### # create a subdirectory for the data from package OmicKriging
dir.create(path = file.path("inst","extdata", "OmicKriging"), recursive = TRUE)
### # copy the grm.bin file into the created subdirectory
file.copy(from = system.file(package = "OmicKriging", "doc/vignette_data/ig_genotypes.grm.bin"),
          to = file.path("inst","extdata", "OmicKriging", "ig_genotypes.grm.bin"))
### # copy the id file
file.copy(from = system.file(package = "OmicKriging", "doc/vignette_data/ig_genotypes.grm.id"),
          to = file.path("inst","extdata", "OmicKriging", "ig_genotypes.grm.id"))

