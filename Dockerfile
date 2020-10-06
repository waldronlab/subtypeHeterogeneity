FROM bioconductor/bioconductor_docker:RELEASE_3_11

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN apt-get update
RUN apt-get -y install libharfbuzz-dev libfribidi-dev  
RUN Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); BiocManager::install(update = TRUE, ask = FALSE)"
RUN Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); devtools::install('.', dependencies = TRUE, repos = BiocManager::repositories())"
