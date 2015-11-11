# Functions common to all cufflinks apps
#
# Cannot use Agave macros in here

# These are not in the Github repo as they are too large
# but are implicitly stored on the storageSystem
tar zxf ./bin.tgz
export PATH="$PATH:$PWD/bin"

tar zxf ./R-3.2.0.tgz
export PATH="$PATH:$PWD/R/bin"
export R_LOCAL=$PWD/R

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
