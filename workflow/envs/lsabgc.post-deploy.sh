set -e

resource_dir="resources"
output_lsabgc="$resource_dir/lsaBGC"
repository="https://github.com/Kalan-Lab/lsaBGC"
version="1.40.0"
release="$repository/archive/refs/tags/v$version.tar.gz"
log="logs/lsabgc/install.log"

mkdir -p $(dirname $log)
mkdir -p $resource_dir

echo "Downloading lsaBGC version $version..." > $log
(cd $resource_dir && wget $release -nc) &>> $log

echo "Extracting lsaBGC..." >> $log
(cd $resource_dir && tar -xzf v$version.tar.gz) &>> $log
(cd $resource_dir && mv lsaBGC-$version lsaBGC) &>> $log
(cd $resource_dir && rm v$version.tar.gz) &>> $log

echo "Installing lsaBGC..." >> $log
(cd $output_lsabgc && python setup.py install) >> $log 2>> $log
(cd $output_lsabgc && pip install -e .) >> $log &>> $log
