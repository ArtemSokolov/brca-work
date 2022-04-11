FROM rocker/tidyverse:latest

RUN R -e 'install.packages("languageserver")'
RUN R -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN R -e 'devtools::install_github( "ArtemSokolov/synExtra" )'

# USER gitpod
# RUN echo 'export PS1="\e[01;34m\w\e[0m$ "' >> $HOME/.bashrc
