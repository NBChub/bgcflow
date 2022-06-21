#!/bin/bash

# exit on the first error encountered
set -e

# create log folder
mkdir -p workflow/report/logs/bigslice
LOG='workflow/report/logs/bigslice/install_bigslice.log'

# get the bigslice models
(cd resources && download_bigslice_hmmdb && rm bigslice_models.tar.gz) 2>> $LOG
bigslice --version &>> $LOG