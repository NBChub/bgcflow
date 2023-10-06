#!/bin/bash

# Install devtools package from a specific CRAN mirror
R -e 'install.packages("devtools", repos="https://mirrors.dotsrc.org/cran/")'

# Install thacklr from GitHub
R -e 'devtools::install_github("thackl/thacklr")'

# Install gggenomes from GitHub
R -e 'devtools::install_github("thackl/gggenomes")'
