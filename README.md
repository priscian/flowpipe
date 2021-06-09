# flowpipe
Flow & mass cytometry analysis pipeline

## Preliminaries
The *flowpipe* R package is fairly easy to set up. In an R session:
```
install.packages("remotes") # If necessary.
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true") # https://github.com/r-lib/remotes#environment-variables
remotes::install_github("priscian/flowpipe")
library(flowpipe)

## Once the package has been installed as described above, all you need to use it is:
library(flowpipe)
```

Well, that's not quite *all* you need. I'm working on a vignette to describe using *flowpipe* that includes details about how to (fairly readily) adapt it for your own data. Using *flowpipe* does require some R/bioinformatics knowledge (esp. to make use of *flowpipe*'s code injection points for dynamic evaluation), but it's not a deep dive, and the effort could produce substantial reductions in time spent on general analyses for investigators.

Meanwhile, [here](inst/docs/cyto-2021-poster_20210520.pdf)'s our poster from [CYTO 2021](https://cytoconference.org/) describing some of *flowpipe*'s features.
