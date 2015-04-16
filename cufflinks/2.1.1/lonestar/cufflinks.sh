# convenience function
echoerr() { echo -e "$@" 1>&2; }

# Hard coded
output_dir=./cufflinks_out

# JMF Note: It is not recommended practice to make the internal variable the
# same as the Agave variable. It can lead to some unexpected substitutions 
# (especially if you modify the value internally)

# Conditionals
multiReadCorrect=${multiReadCorrect}
upperQuartileNorm=${upperQuartileNorm}
totalHitsNorm=${totalHitsNorm}
compatibleHitsNorm=${compatibleHitsNorm}

# Mandatory
query1=${query1}
minIsoformFraction=${minIsoformFraction}
preMrnaFraction=${preMrnaFraction}
maxIntronLength=${maxIntronLength}
smallAnchorFraction=${smallAnchorFraction}
overhangTolerance=${overhangTolerance}
maxBundleLength=${maxBundleLength}
minIntronLength=${minIntronLength}
trim3avgcovThresh=${trim3avgcovThresh}
trim3dropoffFrac=${trim3dropoffFrac}
minFragsPerTransfrag=${minFragsPerTransfrag}
JOB=${jobName}

# These three options only have affect when used with -g/--GTF-guide (RABT assembly)
intronOverhangTolerance=${intronOverhangTolerance}
overhangTolerance3=${overhangTolerance3}
noFauxReads=${noFauxReads}

libraryType=${libraryType}

LABEL=${label}
# Force replace spaces with empty characters
LABEL=${LABEL//\ /}

# Optional file inputs
GTF=${gtf} #bool
GUIDE=${guide_gtf} #bool
#MASK=${MASK_GTF}
BIAS=${BIAS_FASTA}

# Annotation will only be used if GTF or GUIDE is set to true
ANNOTATION=${annotation}

# Create local temp directory
export TMPDIR="$PWD/tmp"
mkdir $TMPDIR
export OMP_NUM_THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)

tar xzf bin.tgz
export PATH="${PATH}:$PWD/bin"

# Comprised of all mandatory or non-dynamically-computed options
OPTIONS="--no-update-check -o $output_dir -p $OMP_NUM_THREADS --min-isoform-fraction ${minIsoformFraction} --pre-mrna-fraction ${preMrnaFraction} --max-intron-length ${maxIntronLength} --small-anchor-fraction ${smallAnchorFraction} --min-frags-per-transfrag ${minFragsPerTransfrag} --overhang-tolerance ${overhangTolerance} --max-bundle-length ${maxBundleLength} --min-intron-length ${minIntronLength} --trim-3-avgcov-thresh ${trim3avgcovThresh} --trim-3-dropoff-frac ${trim3dropoffFrac}"

# Fetch files
# Input
query1_F=${query1}

# Reference GTF (optional)
#if [[ -n $ANNOTATION ]]; then
if [[ $GTF ]]; then
    ANNO_F=${ANNOTATION}
    OPTIONS="${OPTIONS} --GTF ${ANNO_F}"
fi
# Mask GTF (optional)
if [[ -n $MASK ]]; then
    MASK_F=${MASK}
    OPTIONS="${OPTIONS} --mask-file ${MASK_F}"
fi
# Guide GTF (optional)
if [[ $GUIDE ]]; then
    ANNO_F=${ANNOTATION}
    OPTIONS="${OPTIONS} --GTF-guide ${ANNO_F}"
fi
# Bias Fasta (optional)
if [[ -n $BIAS ]]; then
    BIAS_F=${BIAS}
    OPTIONS="${OPTIONS} --frag-bias-correct ${BIAS_F}"
fi

# Conditional or optional params
if [[ $multiReadCorrect ]];                     then OPTIONS="${OPTIONS} --multi-read-correct"; fi
if [[ $upperQuartileNorm ]];                    then OPTIONS="${OPTIONS} --upper-quartile-norm"; fi
if [[ $compatibleHitsNorm && $GTF ]];           then OPTIONS="${OPTIONS} --compatible-hits-norm"
elif [[ "$totalHitsNorm" ]];                    then OPTIONS="${OPTIONS} --total-hits-norm"; fi
if [[ -n $libraryType ]];                       then OPTIONS="${OPTIONS} --library-type ${libraryType}"; fi
if [[ -n $LABEL ]];                             then OPTIONS="${OPTIONS} --label ${LABEL}"; fi
if [[ -n $overhangTolerance3 && $GUIDE ]];      then OPTIONS="${OPTIONS} --3-overhang-tolerance ${overhangTolerance3}"; fi
if [[ -n $intronOverhangTolerance && $GUIDE ]]; then OPTIONS="${OPTIONS} --intron-overhang-tolerance ${intronOverhangTolerance}"; fi
if [[ $noFauxReads && $GUIDE ]];                then OPTIONS="${OPTIONS} --no-faux-reads"; fi

echoerr "CUFFLINKS OPTIONS: $OPTIONS $query1_F\n"

cufflinks $OPTIONS $query1_F 2>cufflinks.stderr
name=${query1_F/\.*/}
echoerr "My base name is $name"
mv ${output_dir}/transcripts.gtf ${output_dir}/${name}-${JOB}.gtf

# convert the transcripts file to bam for gbrowse upload
gtf="${output_dir}/${name}-${JOB}.gtf"
sam="${output_dir}/${name}-${JOB}.sam"
bam="${output_dir}/${name}-${JOB}"
gtf_to_sam $gtf $sam
samtools faidx $BIAS_F
samtools view -b -h -t $BIAS_F.fai -o $bam.unsorted.bam $sam
samtools sort $bam.unsorted.bam $bam
samtools index $bam.bam

rm -rf $bam.unsorted.bam bin tmp *.fa* *.bam *.gtf $sam
