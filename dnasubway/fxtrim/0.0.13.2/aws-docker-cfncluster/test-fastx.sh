#!/bin/bash
#$ -V
#$ -cwd # Start job in submission directory
#$ -N test-fastx # Job Name
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12
#$ -q development       # Queue name normal
#$ -A iPlant-Collabs
#$ -l h_rt=01:00:00     # Run time (hh:mm:ss) 
#$ -M jcarson@tacc.utexas.edu   # Address for email notification
#$ -m be        # Email at Begin and End of job

module purge
module load TACC
module load fastx_toolkit/0.0.13.2

# for DNA Subway w/ Agave

# inputs (for testing)
SEQ1="ra1_rep1_1.fastq.gz"

# Collection of fastx QC utilities for DNA Subway

JOB="tophat"

# fastx trimmer
#FIRST=${first}
#LAST=${last}
#FTRIM=
#if [[ -n $FIRST ]] || [[ -n $LAST ]]; then
#    FTRIM=1
#fi

# fastx quality trimmer
TQUAL=20
MINLEN=20
QTRIM=
if [[ -n $TQUAL ]] || [[ -n $MINLEN ]]; then
    QTRIM=1
fi

# fastx quality filter
MQUAL=20
MPERCENT=50
QFILT=
if [[ -n $MQUAL ]] || [[ -n $MPERCENT ]]; then
    QFILT=1
fi

INDEXING=0

tar zxf FastQC.tgz
tar zxf bin.tgz

export PATH=$PWD/FastQC:$PATH

# A little function to print to STDERR
echoerr() { echo -e "$@\n" 1>&2; }

# output directory
OUTDIR=./fastx_out 
mkdir $OUTDIR

infile=$(basename "$SEQ1")
outfile=${infile/\.*/}
outfile="${OUTDIR}/$outfile";

# Is our fastq file zipped?
# either way, we'll compress the output file..
zipper=gzip
if [[ "$infile" =~ ".gz" ]]; then
    echoerr "Decompressing $infile with $zipper"
    $zipper -d "$infile"
    infile=${infile//.gz/}
fi
if [[ "$infile" =~ ".bz2" ]]; then
    zipper=bzip2
    echoerr "Decompressing $infile with $zipper"
    $zipper -d "$infile"
    infile=${infile//.bz2/}
fi
basename="$infile";

if [[ -n "$INDEXING" && ( $INDEXING == "1" || $INDEXING == "true" ) ]]; then
    echo indexing...

    indexed_file="${outfile}-idx.fastq"
    perl bin/reindex_fq.pl $infile > $indexed_file

echo $indexed_file

    sleep 1
    mv $indexed_file $infile
echo $infile

fi

ARGS=

# are we Sanger or not?
Q=$(bin/check_qual_score.pl "$infile") 

if [[ -n $Q ]]; then
    echoerr "Quality scaling is $Q
";
    ARGS=$Q;
fi


#if [[ -n $FTRIM ]]; then
#    if [[ -n $FIRST ]]; then
#	LARGS="$ARGS -f $FIRST";
#    fi
#    if [[ -n $LAST ]]; then
#	LARGS="$LARGS -l $LAST";
#    fi
#
#    outfile="${outfile}.trim";
#    echoerr "Running fastx_trimmer $LARGS -i $infile -o ${outfile}..."
#    fastx_trimmer $LARGS -i $infile -o $outfile
#    echoerr "Done!
#        "
#    infile=$outfile
#
#    LARGS=''
#fi

if [[ -n $QTRIM ]]; then 
    if [[ -n $TQUAL ]]; then
	LARGS="$ARGS -t $TQUAL"
    fi
    if [[ -n $MINLEN ]]; then
	LARGS="$LARGS -l $MINLEN"
    fi


    outfile="${outfile}.qtrim"
    echoerr "Running fastq_quality_trimmer $LARGS -i $infile -o ${outfile}..."
    fastq_quality_trimmer $LARGS -i "$infile" -o "$outfile"
    echoerr "Done!
        "
    infile="$outfile"

    LARGS=''
fi
    
if [[ -n $QFILT ]]; then
    if [[ -n $MQUAL ]]; then
        LARGS="$ARGS -q $MQUAL"
    fi
    if [[ -n $MPERCENT ]]; then
        LARGS="$LARGS -p $MPERCENT"
    fi

    outfile="${outfile}.qfilt"
    echoerr "Running fastq_quality_filter $LARGS -i $infile -o ${outfile}..."
    fastq_quality_filter $LARGS -i $infile -o $outfile
    echoerr "Done!
        "
fi

mv $outfile ${outfile}.fq
outfile="${outfile}.fq"

base=${basename/\.*/}
final_outfile="${OUTDIR}/${base}-${JOB}.fastq"
cp "$outfile" "$final_outfile"

echoerr "\n*** OUTFILE is $final_outfile\n\n"
cd $OUTDIR

file_to_check=$(basename "$final_outfile");

if  [[ -z $file_to_check  ]]; then
    echoerr "hey, the outfile $file_to_check is empty!\n";
    exit 1;
fi

if  ! [[ -e "$file_to_check"  ]]; then
    echoerr "Hey, the outfile $file_to_check does not exist\n";
    exit 1;
fi


# Since we have trimmed, let's remove the length assertion
# from the sequence headers
echoerr "Removing explit read-lengths from header"
perl -i -pe 's/Length=\d+//i' "$file_to_check"

../FastQC/fastqc --quiet -f fastq "$file_to_check"

if [[ -n $zipper ]]; then
    echoerr "compressing output file with $zipper"
    $zipper "$file_to_check"
fi

cd ..

rm -f ./${basename}*
rm -fr ${OUTDIR}/*filt* ${OUTDIR}/*trim* ${OUTDIR}/*fastqc
rm -rf ./bin ./FastQC


