#!/bin/bash

# exit on the first error encountered
set -e

# create log folder
mkdir -p workflow/report/logs/arts

# get the latest arts release
wget -P resources/ https://bitbucket.org/ziemertlab/arts/get/master.tar.gz 2>> workflow/report/logs/arts/arts_setup.log
(cd resources/ && tar -xvzf master.tar.gz) 2>> workflow/report/logs/arts/arts_setup.log
mv resources/ziemertlab-arts-* resources/arts 2>> workflow/report/logs/arts/arts_setup.log
mkdir -p data/interim/arts/tmp/ 2>> workflow/report/logs/arts/arts_setup.log
echo "ARTS_RESULTS=data/interim/arts/tmp/" > resources/arts/.env
echo "ARTS_UPLOAD=data/interim/arts/tmp/" >> resources/arts/.env
echo "ARTS_RUN=data/interim/arts/tmp/" >> resources/arts/.env
echo "ARTS_CPU=8" >> resources/arts/.env
echo "ARTS_WEBPORT=80" >> resources/arts/.env

# get ranger-dtl binaries
wget -P resources/arts http://compbio.mit.edu/ranger-dtl/ranger-dtl-linux.tar.gz 2>> workflow/report/logs/arts/arts_setup.log

# extract & add executables to path
(cd resources/arts && tar -xvzf ranger-dtl-linux.tar.gz && chmod +x ranger-dtl-*.linux) 2>> workflow/report/logs/arts/arts_setup.log
ln resources/arts/ranger-dtl-U.linux $CONDA_PREFIX/bin/ranger-dtl-U.linux 2>> workflow/report/logs/arts/arts_setup.log
ln resources/arts/ranger-dtl-D.linux $CONDA_PREFIX/bin/ranger-dtl-D.linux 2>> workflow/report/logs/arts/arts_setup.log

# clean up
rm resources/arts/ranger-dtl-linux.tar.gz resources/master.tar.gz