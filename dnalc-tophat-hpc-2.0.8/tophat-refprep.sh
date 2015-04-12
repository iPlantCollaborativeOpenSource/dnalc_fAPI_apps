#----
# Script for creating tar bundle of reference fasta and its index files for tophat/bowtie
# Creates index files
# Bundles everything up
#----

ref=${referenceFasta}
cleanup=${cleanupParameter}
echo "Clean up parameter is now set to:"
echo ${cleanup}

# is input compressed? if so, uncompress and update with new filename
ref_ext="${ref##*.}"
if [ ${ref_ext} = "gz" ]; then gunzip ${ref}; ref=${ref%.*}; fi

# standardize extension to .fa
rb="${ref%%.*}"   # reference basename
mv ${ref} ${rb}.fa
ref=${rb}.fa

gzip_output=1
if [ ${gzip_output} = 1 ]; then compress_flag="z"; else compress_flag=""; fi

# Bowtie is bundled with this application
tar xzf bin.tgz
chmod -R a+x bin/*
export CWD=${PWD}
export PATH="${CWD}/bin:${PATH}"

samtools faidx ${ref}
bowtie2-build -q ${ref} ${ref}

java -Xmx3g -jar $TACC_PICARD_DIR/CreateSequenceDictionary.jar R= ${ref} O= ${rb}.dict
chmod a+rw ${ref}*
tar -c${compress_flag}f ${rb}.tar ${ref} ${ref}.fai ${rb}.dict ${ref}.1.bt2 ${ref}.2.bt2 ${ref}.3.bt2 ${ref}.4.bt2 ${ref}.rev.1.bt2 ${ref}.rev.2.bt2

if [ ${cleanupParameter} ]; then echo "Cleaning up input and intermediate files"
	rm ${ref}; rm ${ref}.fai; rm ${rb}.dict
	rm -rf *.bt2
	rm -rf bin
else
	echo "Skipping clean up"
fi