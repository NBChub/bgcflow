#!/bin/bash

# exit on the first error encountered
set -e

# create log folder
LOG_DIR="logs/bigslice"
mkdir -p logs/bigslice
LOG="$LOG_DIR/install_bigslice.log"

# get the bigslice models
wget -O resources/bigslice_models.tar.gz https://github.com/medema-group/bigslice/releases/download/v1.0.0/bigslice-models.2020-04-27.tar.gz 2> $LOG
mkdir -p $CONDA_PREFIX/bin/bigslice-models 2>> $LOG
tar -xvf resources/bigslice_models.tar.gz -C $CONDA_PREFIX/bin/bigslice-models 2>> $LOG
rm resources/bigslice_models.tar.gz 2>> $LOG
bigslice --version &>> $LOG
