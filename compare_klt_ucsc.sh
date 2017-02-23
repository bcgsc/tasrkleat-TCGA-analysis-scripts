# these are the only four features available in klt.gz and ucsc.gz
FEATURES="exon CDS start_codon stop_codon"
OUTPUT=__diff_klt_ucsc
mkdir -p ${OUTPUT}

# First extract features from both files
for gtf_prefix in klt ucsc; do
    for fea in ${FEATURES}; do
        zcat ${gtf_prefix}.gz \
            | \grep -P "^chr1\t.*\t${fea}\t" \
            | sed 's/^chr//g' \
            | awk '{print $1,$3,$4,$5,$7}' \
            | tr ' ' '\t' \
            | sort -k 1n -k 3n -k 4n \
            > ${OUTPUT}/${gtf_prefix}_${fea}.tsv &
    done
done
wait

# Then compare them, respectively, which should print all 0
for fea in ${FEATURES}; do
    diff ${OUTPUT}/klt_${fea}.tsv ${OUTPUT}/ucsc_${fea}.tsv
    echo $?
done

# cleanup
rm -fvr ${OUTPUT}
