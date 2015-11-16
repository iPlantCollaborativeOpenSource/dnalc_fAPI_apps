# Cufflinks apps DNA Subway
# Ported to use Docker image iplantc/dnasub_apps
# Scripts at /opt/scripts/cufflinks
# Binaries at /opt/bin

# Set up
. $PWD/app_begin.sh

# Hard coded
output_dir=cufflinks_out

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

# Comprised of all mandatory or non-dynamically-computed options
OPTIONS="--no-update-check -o $output_dir -p $THREADS --min-isoform-fraction ${minIsoformFraction} --pre-mrna-fraction ${preMrnaFraction} --max-intron-length ${maxIntronLength} --small-anchor-fraction ${smallAnchorFraction} --min-frags-per-transfrag ${minFragsPerTransfrag} --overhang-tolerance ${overhangTolerance} --max-bundle-length ${maxBundleLength} --min-intron-length ${minIntronLength} --trim-3-avgcov-thresh ${trim3avgcovThresh} --trim-3-dropoff-frac ${trim3dropoffFrac}"

# Inputs
query1_F=${query1}

# Inputs + parameters
# Reference GTF (optional)
if [[ $GTF -eq 1 ]]; then
 ANNO_F=${ANNOTATION}
 OPTIONS="${OPTIONS} --GTF ${ANNO_F}"
fi
# Mask GTF (optional)
if [[ -n $MASK ]]; then
 MASK_F=${MASK}
 OPTIONS="${OPTIONS} --mask-file ${MASK_F}"
fi
# Guide GTF (optional)
if [[ $GUIDE -eq 1 ]]; then
 ANNO_F=${ANNOTATION}
 OPTIONS="${OPTIONS} --GTF-guide ${ANNO_F}"
fi
# Bias Fasta (optional)
if [[ -n $BIAS ]]; then
 BIAS_F=${BIAS}
 OPTIONS="${OPTIONS} --frag-bias-correct ${BIAS_F}"
fi

# Conditional or optional params
if [[ $multiReadCorrect -eq 1 ]]; then OPTIONS="${OPTIONS} --multi-read-correct"; fi
if [[ $upperQuartileNorm -eq 1 ]]; then OPTIONS="${OPTIONS} --upper-quartile-norm"; fi
if [[ $compatibleHitsNorm -eq 1 && $GTF -eq 1 ]]; then OPTIONS="${OPTIONS} --compatible-hits-norm"
elif [[ $totalHitsNorm -eq 1 ]]; then OPTIONS="${OPTIONS} --total-hits-norm"; fi
if [[ -n $libraryType ]]; then OPTIONS="${OPTIONS} --library-type ${libraryType}"; fi
if [[ -n $LABEL ]]; then OPTIONS="${OPTIONS} --label ${LABEL}"; fi
if [[ -n $overhangTolerance3 && $GUIDE -eq 1 ]]; then OPTIONS="${OPTIONS} --3-overhang-tolerance ${overhangTolerance3}"; fi
if [[ -n $intronOverhangTolerance && $GUIDE -eq 1 ]]; then OPTIONS="${OPTIONS} --intron-overhang-tolerance ${intronOverhangTolerance}"; fi
if [[ $noFauxReads -eq 1 && $GUIDE -eq 1 ]]; then OPTIONS="${OPTIONS} --no-faux-reads"; fi

echoerr "Command: cufflinks $OPTIONS $query1_F"

$DOCKER_APP_RUN cufflinks ${OPTIONS} ${query1_F} 2>cufflinks.stderr

name=${query1_F/\.*/}
echoerr "My base name is $name"
mv ${output_dir}/transcripts.gtf ${output_dir}/${name}-${JOB}.gtf

# Convert the transcripts file to bam for gbrowse upload
gtf="${output_dir}/${name}-${JOB}.gtf"
sam="${output_dir}/${name}-${JOB}.sam"
bam="${output_dir}/${name}-${JOB}"

$DOCKER_APP_RUN gtf_to_sam $gtf $sam
$DOCKER_APP_RUN samtools faidx $BIAS_F
$DOCKER_APP_RUN samtools view -b -h -t $BIAS_F.fai -o $bam.unsorted.bam $sam
$DOCKER_APP_RUN samtools sort $bam.unsorted.bam $bam
$DOCKER_APP_RUN samtools index $bam.bam

# Teardown
rm -rf $bam.unsorted.bam *.fa* *.bam *.gtf $sam
. $PWD/app_end.sh

#if [ ! -s $bam.bam ];
#then
#  echoerr "$bam.bam had zero length. This can happen if the parameters are not set correctly."
#  ${AGAVE_JOB_CALLBACK_FAILURE}
#fi
