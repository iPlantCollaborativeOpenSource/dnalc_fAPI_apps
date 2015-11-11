# Binaries and PATH
tar -zxf bin.tgz
tar -zxf FastQC.tgz
export PATH=$PWD/bin:$PWD/FastQC:$PATH
# Directories
export CWD=${PWD}
export TMPDIR="$CWD/tmp"
mkdir $TMPDIR

# Should probably use the $NSLOTS environment variable to calculate thread info if possible
THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)
if [[ $THREADS -eq 0 ]];
then
  THREADS=4
fi
export THREADS

# Convenience functions
echoerr() { echo -e "$@" 1>&2; }
