export CWD=${PWD}
export PATH="${CWD}/bin:${PATH}"
export TMPDIR="$CWD/tmp"
mkdir -p $TMPDIR

# Should probably use the $NSLOTS environment variable to calculate thread info if possible
THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)
if [[ $THREADS -eq 0 ]];
then
  THREADS=4
fi
export THREADS

# Convenience functions
echoerr() { echo -e "$@\n" 1>&2; }

# Docker
DOCKER_APP_IMAGE=iplantc/dnasub_apps:latest
HOST_SCRATCH=/home
. $PWD/docker_begin.sh
