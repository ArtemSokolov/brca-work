FROM rocker/tidyverse:latest

RUN R -e 'install.packages( "synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org") )'
RUN R -e 'devtools::install_github( "ArtemSokolov/synExtra" )'
