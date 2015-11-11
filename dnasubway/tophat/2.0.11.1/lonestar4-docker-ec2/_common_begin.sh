
tar -xzf bin.tgz
chmod -R a+x bin/*
export CWD=${PWD}
export PATH="${CWD}/bin:${PATH}"
export TMPDIR="$CWD/tmp"
mkdir -p $TMPDIR

# Automatically find out how many cores the node
# has and create this many threads
THREADS=`cat /proc/cpuinfo | grep processor | wc -l`
# If the above trick has failed, fall back to 4 threads
if [[ $THREADS == 0 ]]; then
    THREADS=4
fi
export THREADS

echoerr() { echo -e "$@" 1>&2; }
