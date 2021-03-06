# Cufflinks apps DNA Subway
# Ported to use Docker image iplantc/dnasub_apps
# Scripts at /opt/scripts/cufflinks
# Binaries at /opt/bin

# This is a specialized implementation of the Docker app
# design pattern for Agave apps where I use /dev/shm as volume
# to work around the fact that sqlite3 cannot work on the
# NFS mounted directory we are using as a working directory
# for processing a cuffdiff analysis

OUTPUT_DIR=cuffdiff_out

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

## -> NO USER-SERVICABLE PARTS INSIDE
if [ -z "${DOCKER_APP_IMAGE}" ];
then
    echo "Error: DOCKER_APP_IMAGE was not defined"
fi
if [ -z "$HOST_SCRATCH" ];
then
    echo "Warning: HOST_SCRATCH was not defined"
fi

#TTY="-t"
TTY=""
MYUID=$(id -u $USER)
STAMP=$(date +%s)
UUID=$(uuidgen -r)
HOSTNAME=$(uname -n)
HOST_OPTS="-u=$MYUID -m 24G"
# Create unique but human comprehensible name for each container
DOCKER_APP_CONTAINER="appc-${UUID}"
PERSIST=172800

# Make a dedicated directory in /dev/shm for the app container
# If you don't do this, multiple containers accessing /dev/shm
# can cause a race condition, leading to app failure, node crash,
# or docker subsystem malfunction
SHM="/dev/shm/appc-${UUID}"
mkdir -p $SHM
if [ ! -d $SHM ];
then
  echoerr "Couldn't create $SHM. Exiting."
  ${AGAVE_JOB_CALLBACK_FAILURE}
  exit 1
fi

# Logging
echo "hostname: $HOSTNAME" >> docker.log
echo "username: $USER" >> docker.log
echo "userid: $MYUID" >> docker.log
echo "container_name: $DOCKER_APP_CONTAINER" >> docker.log

# Refresh the image
docker pull ${DOCKER_APP_IMAGE} | grep "Status" >> docker.log

# Launch the app container, persist up to 24h
DOCKER_APP_CREATE="docker run ${HOST_OPTS} -d -v `pwd`:${HOST_SCRATCH}:rw -w ${HOST_SCRATCH} -v ${SHM}:/home/$OUTPUT_DIR --name ${DOCKER_APP_CONTAINER} ${DOCKER_APP_IMAGE} sleep ${PERSIST}"

# Export Container ID
export DOCKER_APP_CONTAINER=$( ${DOCKER_APP_CREATE} | cut -c1-12)
echo "container_id: $DOCKER_APP_CONTAINER" >> docker.log
# Append this in front of every command
export DOCKER_APP_RUN="docker exec -i ${TTY} $DOCKER_APP_CONTAINER"
## <- NO USER-SERVICABLE PARTS INSIDE

JOB=${jobName}

# GTF files to merge
QUERY1=${query1}
QUERY2=${query2}
QUERY3=${query3}
QUERY4=${query4}
QUERY5=${query5}
QUERY6=${query6}
QUERY7=${query7}
QUERY8=${query8}
QUERY9=${query9}
QUERY10=${query10}
QUERY11=${query11}
QUERY12=${query12}

# Up to ten samples, up to four replicates each
# sam1_1,sam1_2,sam1_3,sam1_4...sam10_1,sam10_2,sam10_3,sam10_4
SAM1_F1=${sam1_f1}
SAM1_F2=${sam1_f2}
SAM1_F3=${sam1_f3}
SAM1_F4=${sam1_f4}
SAM2_F1=${sam2_f1}
SAM2_F2=${sam2_f2}
SAM2_F3=${sam2_f3}
SAM2_F4=${sam2_f4}
SAM3_F1=${sam3_f1}
SAM3_F2=${sam3_f2}
SAM3_F3=${sam3_f3}
SAM3_F4=${sam3_f4}
SAM4_F1=${sam4_f1}
SAM4_F2=${sam4_f2}
SAM4_F3=${sam4_f3}
SAM4_F4=${sam4_f4}
SAM5_F1=${sam5_f1}
SAM5_F2=${sam5_f2}
SAM5_F3=${sam5_f3}
SAM5_F4=${sam5_f4}
SAM6_F1=${sam6_f1}
SAM6_F2=${sam6_f2}
SAM6_F3=${sam6_f3}
SAM6_F4=${sam6_f4}
SAM7_F1=${sam7_f1}
SAM7_F2=${sam7_f2}
SAM7_F3=${sam7_f3}
SAM7_F4=${sam7_f4}
SAM8_F1=${sam8_f1}
SAM8_F2=${sam8_f2}
SAM8_F3=${sam8_f3}
SAM8_F4=${sam8_f4}
SAM9_F1=${sam9_f1}
SAM9_F2=${sam9_f2}
SAM9_F3=${sam9_f3}
SAM9_F4=${sam9_f4}
SAM10_F1=${sam10_f1}
SAM10_F2=${sam10_f2}
SAM10_F3=${sam10_f3}
SAM10_F4=${sam10_f4}

# Supplemental files
REFGTF=${ref_gtf}
REFSEQ=${ref_seq}

# --mask-file optional
MASK=${mask_gtf}
SKIPCUFFMERGE=${skipCuffmerge}
USEGTF=${refGTF}

MANIFEST=gtf_to_merge.txt
touch $MANIFEST
if [[ -n $QUERY1  ]]; then echo $QUERY1  >> $MANIFEST; fi
if [[ -n $QUERY2  ]]; then echo $QUERY2  >> $MANIFEST; fi
if [[ -n $QUERY3  ]]; then echo $QUERY3  >> $MANIFEST; fi
if [[ -n $QUERY4  ]]; then echo $QUERY4  >> $MANIFEST; fi
if [[ -n $QUERY5  ]]; then echo $QUERY5  >> $MANIFEST; fi
if [[ -n $QUERY6  ]]; then echo $QUERY6  >> $MANIFEST; fi
if [[ -n $QUERY7  ]]; then echo $QUERY7  >> $MANIFEST; fi
if [[ -n $QUERY8  ]]; then echo $QUERY8  >> $MANIFEST; fi
if [[ -n $QUERY9  ]]; then echo $QUERY9  >> $MANIFEST; fi
if [[ -n $QUERY10 ]]; then echo $QUERY10 >> $MANIFEST; fi
if [[ -n $QUERY11 ]]; then echo $QUERY11 >> $MANIFEST; fi
if [[ -n $QUERY12 ]]; then echo $QUERY12 >> $MANIFEST; fi

if [[ $SKIPCUFFMERGE -eq 1 ]]; then
    echoerr "Skipping Cuffmerge"
else
    echoerr "Inspecting list of GTF files to merge..."
    lines=$(cat $MANIFEST | wc -l)
    unique=$(sort -u $MANIFEST | wc -l)

    if ! [[ $unique -eq $lines ]]; then
        echoerr "Error: Each GTF file to be merged must have a unique filename"
        ${AGAVE_JOB_CALLBACK_FAILURE}
        exit 1
    fi

    if [[ $unique -lt 2 ]]; then
        echoerr "Error: at least two GTF files are required for merging"
        ${AGAVE_JOB_CALLBACK_FAILURE}
        exit 1
    fi

    if [[ $USEGTF -eq 1 ]]; then
        INCLUDE_REF_GTF=" -g $REFGTF "
        echoerr "Reference GTF being used at Cuffmerge"
    else
        INCLUDE_REF_GTF=""
        echoerr "No reference GTF being used at Cuffmerge"
    fi

    ARGS="-p $THREADS -o cuffmerge_out $INCLUDE_REF_GTF -s $REFSEQ $MANIFEST";

    echoerr "Command: cuffmerge ${ARGS}"
    $DOCKER_APP_RUN cuffmerge $ARGS

    MERGEDONE=$(ls -alh cuffmerge_out)
    echoerr "DONE!
    $MERGEDONE
    "

    ANNOTATION='./cuffmerge_out/merged.gtf'
    if [[ ! -e $ANNOTATION ]]; then
        echoerr "No point going on, cuffmerge failed"
        ${AGAVE_JOB_CALLBACK_FAILURE}
        exit 1
    fi
fi

if [[ -z $ANNOTATION ]]; then
    ANNOTATION=$REFGTF
fi

# not supported/recommended
#poissonDispersion=${poissonDispersion}

# Mandatory parameters
# --min-alignment-count 10
MINALIGNMENTCOUNT=${minAlignmentCount}
# --FDR 0.05
FDR=${fdr}
# Force Replace spaces with empty characters
LABELS=${labels}
LABELS=${LABELS//\ /}

# --library-type fr-unstranded
LIBRARYTYPE=${libraryType}

# Optional
# --frag-len-mean 200
FRAGLENMEAN=${fragLenMean}
# --frag-len-std-dev 80
FRAGLENSTDEV=${fragLenStdev}

# Fetch alignment files
if [[ -e "$SAM1_F1" ]]; then SAM1_F=$SAM1_F1; echo "$SAM1_F found"
else echoerr "$SAM1_F is missing! Abort!"; exit 1; fi

if [[ -n $SAM1_F2 ]]; then SAM1_F="$SAM1_F,$SAM1_F2"; fi
if [[ -n $SAM1_F3 ]]; then SAM1_F="$SAM1_F,$SAM1_F3"; fi
if [[ -n $SAM1_F4 ]]; then SAM1_F="$SAM1_F,$SAM1_F4"; fi
echoerr "SAM1 files $SAM1_F"

if [[ -n $SAM2_F1 ]]; then SAM2_F="$SAM2_F1"        ; fi
if [[ -n $SAM2_F2 ]]; then SAM2_F="$SAM2_F,$SAM2_F2"; fi
if [[ -n $SAM2_F3 ]]; then SAM2_F="$SAM2_F,$SAM2_F3"; fi
if [[ -n $SAM2_F4 ]]; then SAM1_F="$SAM2_F,$SAM2_F4"; fi
echoerr "SAM2 files $SAM2_F"

if [[ -n $SAM3_F1 ]]; then SAM3_F="$SAM3_F1"        ; fi
if [[ -n $SAM3_F2 ]]; then SAM3_F="$SAM3_F,$SAM3_F2"; fi
if [[ -n $SAM3_F3 ]]; then SAM3_F="$SAM3_F,$SAM3_F3"; fi
if [[ -n $SAM3_F4 ]]; then SAM3_F="$SAM3_F,$SAM3_F4"; fi
echoerr "SAM3 Files $SAM3_F"

if [[ -n $SAM4_F1 ]]; then SAM4_F="$SAM4_F1"        ; fi
if [[ -n $SAM4_F2 ]]; then SAM4_F="$SAM4_F,$SAM4_F2"; fi
if [[ -n $SAM4_F3 ]]; then SAM4_F="$SAM4_F,$SAM4_F3"; fi
if [[ -n $SAM4_F4 ]]; then SAM4_F="$SAM4_F,$SAM4_F4"; fi

if [[ -n $SAM5_F1 ]]; then SAM5_F="$SAM5_F1"        ; fi
if [[ -n $SAM5_F2 ]]; then SAM5_F="$SAM5_F,$SAM5_F2"; fi
if [[ -n $SAM5_F3 ]]; then SAM5_F="$SAM5_F,$SAM5_F3"; fi
if [[ -n $SAM5_F4 ]]; then SAM5_F="$SAM5_F,$SAM5_F4"; fi

if [[ -n $SAM6_F1 ]]; then SAM6_F="$SAM6_F1"        ; fi
if [[ -n $SAM6_F2 ]]; then SAM6_F="$SAM6_F,$SAM6_F2"; fi
if [[ -n $SAM6_F3 ]]; then SAM6_F="$SAM6_F,$SAM6_F3"; fi
if [[ -n $SAM6_F4 ]]; then SAM6_F="$SAM6_F,$SAM6_F4"; fi

if [[ -n $SAM7_F1 ]]; then SAM7_F="$SAM7_F1"        ; fi
if [[ -n $SAM7_F2 ]]; then SAM7_F="$SAM7_F,$SAM7_F2"; fi
if [[ -n $SAM7_F3 ]]; then SAM7_F="$SAM7_F,$SAM7_F3"; fi
if [[ -n $SAM7_F4 ]]; then SAM7_F="$SAM7_F,$SAM7_F4"; fi

if [[ -n $SAM8_F1 ]]; then SAM8_F="$SAM8_F1"        ; fi
if [[ -n $SAM8_F2 ]]; then SAM8_F="$SAM8_F,$SAM8_F2"; fi
if [[ -n $SAM8_F3 ]]; then SAM8_F="$SAM8_F,$SAM8_F3"; fi
if [[ -n $SAM8_F4 ]]; then SAM8_F="$SAM8_F,$SAM8_F4"; fi

if [[ -n $SAM9_F1 ]]; then SAM9_F="$SAM9_F1"        ; fi
if [[ -n $SAM9_F2 ]]; then SAM9_F="$SAM9_F,$SAM9_F2"; fi
if [[ -n $SAM9_F3 ]]; then SAM9_F="$SAM9_F,$SAM9_F3"; fi
if [[ -n $SAM9_F4 ]]; then SAM9_F="$SAM9_F,$SAM9_F4"; fi

if [[ -n $SAM10_F1 ]]; then SAM10_F="$SAM10_F1"         ; fi
if [[ -n $SAM10_F2 ]]; then SAM10_F="$SAM10_F,$SAM10_F2"; fi
if [[ -n $SAM10_F3 ]]; then SAM10_F="$SAM10_F,$SAM10_F3"; fi
if [[ -n $SAM10_F4 ]]; then SAM10_F="$SAM10_F,$SAM10_F4"; fi

# Initialize OPTIONS with mandatory parameters
OPTIONS="--no-update-check --num-threads ${THREADS} --output-dir ${OUTPUT_DIR}"
OPTIONS="$OPTIONS --library-type ${LIBRARYTYPE} --min-alignment-count ${MINALIGNMENTCOUNT} --labels ${LABELS}"

# Flag  OPTIONS
if   [[ "${treatAsTimeSeries}"  -eq 1 ]]; then OPTIONS="${OPTIONS} --time-series"              ; fi
if   [[ "${multiReadCorrect}"   -eq 1 ]]; then OPTIONS="${OPTIONS} --multi-read-correct"       ; fi
if   [[ "${upperQuartileNorm}"  -eq 1 ]]; then OPTIONS="${OPTIONS} --upper-quartile-norm"      ; fi
if   [[ "${totalHitsNorm}"      -eq 1 ]]; then OPTIONS="${OPTIONS} --total-hits-norm"
elif [[ "${compatibleHitsNorm}" -eq 1 ]]; then OPTIONS="${OPTIONS} --compatible-hits-norm"     ; fi
if   [[ "${poissonDispersion}"  -eq 1 ]]; then OPTIONS="${OPTIONS} --dispersion-method poisson"; fi

# Parameter OPTIONS
# JMF note: these next few line can't be working as intended.
if [[ -n "${FRAGLENMEAN}" ]]; then
    OPTIONS="${OPTIONS} --frag-len-mean ${FRAGLENMEAN}"
fi
if [[ -n "${FRAGLENSTDEV}" ]]; then
    OPTIONS="${OPTIONS} --frag-len-std-dev ${FRAGLENSTDEV}"
fi

if [[ -n "$fdr" ]]     ; then OPTIONS="${OPTIONS} --FDR $fdr"    ; fi
if [[ -n "$REFSEQ_F" ]]; then OPTIONS="${OPTIONS} -b ${REFSEQ_F}"; fi

# and we pull the trigger...
SAMS="${SAM1_F} ${SAM2_F} ${SAM3_F} ${SAM4_F} ${SAM5_F} ${SAM6_F} ${SAM7_F} ${SAM8_F} ${SAM9_F} ${SAM10_F}"
echoerr "$SAMS"
echoerr "$(ls -al *.bam *.gtf)"
echoerr "Command: cuffdiff ${OPTIONS} ${ANNOTATION} $SAMS"

$DOCKER_APP_RUN cuffdiff ${OPTIONS} ${ANNOTATION} $SAMS 2>cuffdiff.stderr

if ! [[ -d $OUTPUT_DIR ]]; then
  echoerr "Oh no, cuffdiff failed!  I quit'"
  ${AGAVE_JOB_CALLBACK_FAILURE}
  exit 1
fi

echoerr "Preparing global R plots..."
  $DOCKER_APP_RUN perl /opt/scripts/cufflinks/cuffdiff_R_plots.pl $LABELS
echoerr "Sorting output data..."

# Get the biomart annotations
wget "https://agave.iplantc.org/files/v2/download/jfonner/system/data.iplantcollaborative.org/shared/iplant_DNA_subway/genomes/${species}/${species}.txt"

if [[ -e ${species}.txt ]]; then
        echoerr "BioMart annotations exist, grabbing them, then running cuffdiff_sort.pl"
        $DOCKER_APP_RUN perl /opt/scripts/cufflinks/cuffdiff_sort.pl $PWD $LABELS ${species}
else
        echoerr "BioMart annotations do not exist, running cuffdiff_sort.pl without them"
        $DOCKER_APP_RUN perl /opt/scripts/cufflinks/cuffdiff_sort.pl $PWD $LABELS
fi

echoerr "Done!"

# Copy files out of /dev/shm/app-UUID
# Must be done in container context since files aren't
# -technically- visible outside it. If this doesn't work
# we can copy straight from the special /dev/shm directory
# but that breaks the abstraction even more
#
mkdir cuffdiff_out_temp
$DOCKER_APP_RUN cp -R cuffdiff_out/* cuffdiff_out_temp

# Clean up
#rm -frv *.gtf *.fa* *.txt *.bam annotations cuffdiff_out/*db
#rm -rfv tmp bin bin.tgz
#docker rm -f ${DOCKER_APP_CONTAINER}
#rm -rfv $SHM

# Rename the temp directory
# For now, rename to cuffdiff_out_1 until I am sure it's OK to blow away
# any emphemral cuffdiff_out that gets generated by the script
#
#mv cuffdiff_out_temp cuffdiff_out_1

# Now, we're done
