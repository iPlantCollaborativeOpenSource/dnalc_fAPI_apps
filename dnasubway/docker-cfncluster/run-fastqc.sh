# FastQC utilities for DNA Subway
# Ported to use Docker image iplantc/dnasub_apps
# FastQC is at /opt/FastQC/
# Scripts at /opt/scripts/fastqc

# Set up
. $PWD/app_begin.sh

#if [ "$[THREADS]" -gt "${AGAVE_JOB_PROCESSORS_PER_NODE}" ];
#then
#  export THREADS=${AGAVE_JOB_PROCESSORS_PER_NODE}
#fi

# inputs
SEQ1="${input}"
mkdir -p fastqc_out
$DOCKER_APP_RUN /opt/FastQC/fastqc --quiet -f fastq -o fastqc_out "$SEQ1"

# Cleanup
rm "$SEQ1"

# Teardown
. $PWD/app_end.sh
