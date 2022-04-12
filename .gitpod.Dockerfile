FROM rocker/tidyverse:latest

RUN R -e 'install.packages("languageserver")'
RUN R -e 'install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))'
RUN R -e 'remotes::install_github("ArtemSokolov/synExtra")'

RUN R -e 'install.packages(c("ggrepel", "cowplot"))'

### Gitpod user ###
# '-l': see https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#user
RUN useradd -l -u 33333 -G sudo -md /home/gitpod -s /bin/bash -p gitpod gitpod \
    # passwordless sudo for users in the 'sudo' group
    && sed -i.bkp -e 's/%sudo\s\+ALL=(ALL\(:ALL\)\?)\s\+ALL/%sudo ALL=NOPASSWD:ALL/g' /etc/sudoers
ENV HOME=/home/gitpod
WORKDIR $HOME
USER gitpod

RUN echo 'export PS1="\e[01;34m\w\e[0m$ "' >> $HOME/.bashrc
