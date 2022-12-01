#!/bin/bash

# exit on the first error encountered
set -e

# create log folder
mkdir -p workflow/report/logs/arts

# get the latest arts release
rm -rf resources/arts 2> workflow/report/logs/arts/arts_setup.log
git clone git@github.com:NBChub/arts_v3.git resources/arts &> workflow/report/logs/arts/arts_setup.log
mkdir -p data/interim/arts/tmp/ 2>> workflow/report/logs/arts/arts_setup.log

# set up reference database
wget -P resources/ https://bitbucket.org/ziemertlab/arts/get/master.tar.gz -nc 2>> workflow/report/logs/arts/arts_setup.log
(cd resources/ && tar -xvzf master.tar.gz ziemertlab-arts-b4789c6b3a88/reference) 2>> workflow/report/logs/arts/arts_setup.log
mv resources/ziemertlab-arts-*/reference resources/arts/ 2>> workflow/report/logs/arts/arts_setup.log
(cd resources/ && tar -xvzf master.tar.gz ziemertlab-arts-b4789c6b3a88/astral) 2>> workflow/report/logs/arts/arts_setup.log
mv resources/ziemertlab-arts-*/astral resources/arts/ 2>> workflow/report/logs/arts/arts_setup.log
rm resources/master.tar.gz

# set up environment variable
echo "ARTS_RESULTS=data/interim/arts/tmp/" > resources/arts/.env
echo "ARTS_UPLOAD=data/interim/arts/tmp/" >> resources/arts/.env
echo "ARTS_RUN=data/interim/arts/tmp/" >> resources/arts/.env
echo "ARTS_CPU=8" >> resources/arts/.env
echo "ARTS_WEBPORT=80" >> resources/arts/.env

# install binaries
(cd resources/arts && tar -xzf linux64_bins.tar.gz -C  $CONDA_PREFIX/bin && hash -r) 2>> workflow/report/logs/arts/arts_setup.log
