# these are the only four features available in klt.gz and ucsc.gz
FEATURES="exon CDS start_codon stop_codon"
OUTPUT=__diff
mkdir -p ${OUTPUT}

echo "# First, extracting features into separate files..."
# NOTE: ens.gz uses 1 while the other two use chr1)
for gtf_prefix in klt ucsc ens; do
    for fea in ${FEATURES}; do
        zcat ${gtf_prefix}.gtf.gz \
            | \grep -P "^(?:chr)?1\t.*\t${fea}\t" \
            | sed 's/^chr//g' \
            | awk '{print $1,$3,$4,$5,$7}' \
            | tr ' ' '\t' \
            | sort -k 1n -k 3n -k 4n \
            > ${OUTPUT}/${gtf_prefix}_${fea}.tsv &
    done
done
wait

echo 
echo "# Comparing klt.gz & ucsc.gz..."
echo "It should print out all 0s"
for fea in ${FEATURES}; do
    echo -n "comparing ${fea}: "
    diff ${OUTPUT}/klt_${fea}.tsv ${OUTPUT}/ucsc_${fea}.tsv
    echo $?
done

echo
echo "# Comparing klt.gz & ens.gz..."
echo "If the two are the same, the returncode should be 0"
for fea in ${FEATURES}; do
    echo -n "comparing ${fea}: "
    sdiff ${OUTPUT}/klt_${fea}.tsv ${OUTPUT}/ens_${fea}.tsv > ${OUTPUT}/sdiff_${fea}.log
    echo $?
done
