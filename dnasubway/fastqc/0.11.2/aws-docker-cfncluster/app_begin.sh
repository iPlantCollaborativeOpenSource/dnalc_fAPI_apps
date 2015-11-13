# Binaries and PATH
tar -zxf FastQC.tgz
export PATH=$PWD/FastQC:$PATH

# Should probably use the $NSLOTS environment variable to calculate thread info if possible
THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)
if [[ $THREADS -eq 0 ]];
then
  THREADS=4
fi
export THREADS

# Directories
export CWD=${PWD}
export TMPDIR="$CWD/tmp"
mkdir -p $TMPDIR

# Convenience functions
echoerr() { echo -e "$@" 1>&2; }

# Docker
DOCKER_APP_IMAGE=iplantc/dnasub_ls4_env
HOST_SCRATCH=/home
. $PWD/docker_begin.sh
