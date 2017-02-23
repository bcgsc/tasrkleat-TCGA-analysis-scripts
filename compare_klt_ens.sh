# these are the only four features available in klt.gz and ucsc.gz
FEATURES="exon CDS start_codon stop_codon"
OUTPUT=__diff_klt_ens
mkdir -p ${OUTPUT}

# First extract features from both files. (NOTE: ens.gz uses 1 instead of chr1)
for gtf_prefix in klt ens; do
    for fea in ${FEATURES}; do
        zcat ${gtf_prefix}.gz \
            | \grep -P "^(?:chr)?1\t.*\t${fea}\t" \
            | sed 's/^chr//g' \
            | awk '{print $1,$3,$4,$5,$7}' \
            | tr ' ' '\t' \
            | sort -k 1n -k 3n -k 4n \
            > ${OUTPUT}/${gtf_prefix}_${fea}.tsv &
    done
done
wait

# Then compare them, respectively
for fea in ${FEATURES}; do
    echo "If the two are the same, the returncode should be 0"
    echo -n "comparing ${fea}: "
    sdiff ${OUTPUT}/klt_${fea}.tsv ${OUTPUT}/ens_${fea}.tsv > ${OUTPUT}/sdiff_${fea}.tsv
    echo $?
done

# don't cleanup for further inspection if wanted
