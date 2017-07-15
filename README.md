# lof-enrichment

## Cloning and Set Up
Use `git clone --recursive https://github.com/rivas-lab/lof-enrichment.git` to
clone this repo and automatically pull in the submodules as well. If you use
the typical `git clone` and need to pull in the submodules, change into the
repo and run `git submodule update --init --recursive`.

After cloning, change into the repo and use `git lfs pull` in case some files 
tracked with LFS were not downloaded. If you have just installed git LFS for
the first time, you may need to run `git lfs install` from the command line
first.

## Shiny App
First install the GeneticsDesign package in `data` using 
`R CMD install data/GeneticsDesign_1.45.0.tar.gz`. Start RStudio and open
`shinyapp/ui.R` and click the "Run App" button above the file contents. You can
use the arrow next to the "Run App" button to choose whether to launch the app
in RStudio or your browser. Note that the app relies on data that is tracked
with git LFS, so be it might make sense to make sure all the LFS data has been
pulled in using `git lfs pull`.
