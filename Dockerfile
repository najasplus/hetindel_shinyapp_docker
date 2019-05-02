FROM rocker/r-ver:3.6.0

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget


# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown', 'knitr', 'readr'), repos='http://cran.rstudio.com/')" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("sangerseqR", "Biostrings", "BiocGenerics", "msa")' && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    chown shiny:shiny /var/lib/shiny-server



#Copy configuration files into the Docker image
COPY shiny-server.conf /etc/shiny-server.conf 
COPY /app /srv/shiny-server/
COPY shiny-server.sh /usr/bin/shiny-server.sh


#Shiny awailable at port 80
EXPOSE 80

CMD ["/usr/bin/shiny-server.sh"]


