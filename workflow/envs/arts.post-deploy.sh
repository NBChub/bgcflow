#!/bin/bash

# exit on the first error encountered
set -e

mkdir -p resources

# create log folder
ARTS_LOG_DIR="logs/arts"
mkdir -p $ARTS_LOG_DIR
echo "Installing ARTS - saving log files at $ARTS_LOG_DIR" &> $ARTS_LOG_DIR/arts_setup.log

# get the latest arts release
ARTS_REPO="https://github.com/NBChub/arts_v3"
ARTS_RELEASE_VERSION="0.1.0"
echo "STEP 1 - Getting ARTS version $ARTS_RELEASE_VERSION from $ARTS_REPO" &>> $ARTS_LOG_DIR/arts_setup.log
rm -rf resources/arts 2>> $ARTS_LOG_DIR/arts_setup.log
wget -O resources/arts.zip $ARTS_REPO/archive/refs/tags/$ARTS_RELEASE_VERSION.zip -nc -q 2>> $ARTS_LOG_DIR/arts_setup.log
(cd resources && unzip arts.zip && mv arts_v3-$ARTS_RELEASE_VERSION arts) &>> $ARTS_LOG_DIR/arts_setup.log
rm resources/arts.zip
mkdir -p data/interim/arts/tmp/ 2>> $ARTS_LOG_DIR/arts_setup.log

# set up reference database
REF_DB="https://bitbucket.org/ziemertlab/arts/get/master.tar.gz"
echo "STEP 2 - Setting up reference databases from $REF_DB" &>> $ARTS_LOG_DIR/arts_setup.log
wget -P resources/ $REF_DB -nc -q 2>> $ARTS_LOG_DIR/arts_setup.log
(cd resources/ && tar -xvzf master.tar.gz ziemertlab-arts-b4789c6b3a88/reference) &>> $ARTS_LOG_DIR/arts_setup.log
mv resources/ziemertlab-arts-*/reference resources/arts/ 2>> $ARTS_LOG_DIR/arts_setup.log
(cd resources/ && tar -xvzf master.tar.gz ziemertlab-arts-b4789c6b3a88/astral) &>> $ARTS_LOG_DIR/arts_setup.log
mv resources/ziemertlab-arts-*/astral resources/arts/ 2>> $ARTS_LOG_DIR/arts_setup.log
rm resources/master.tar.gz

# set up environment variable
echo "STEP 3 - Setting up environment variables" &>> $ARTS_LOG_DIR/arts_setup.log
echo "ARTS_RESULTS=data/interim/arts/tmp/" > resources/arts/.env
echo "ARTS_UPLOAD=data/interim/arts/tmp/" >> resources/arts/.env
echo "ARTS_RUN=data/interim/arts/tmp/" >> resources/arts/.env
echo "ARTS_CPU=8" >> resources/arts/.env
echo "ARTS_WEBPORT=80" >> resources/arts/.env

# install binaries
echo "STEP 4 - Setting up binaries" &>> $ARTS_LOG_DIR/arts_setup.log
(cd resources/arts && tar -xzf linux64_bins.tar.gz -C  $CONDA_PREFIX/bin && hash -r) 2>> $ARTS_LOG_DIR/arts_setup.log

echo "Done!" &>> $ARTS_LOG_DIR/arts_setup.log
