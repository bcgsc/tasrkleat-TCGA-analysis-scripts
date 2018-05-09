for bam in $(\grep 'reads2genome\/cba.sorted.bam$' ls-output-local.txt); do
    output_dir=$(echo $(dirname ${bam}) | sed 's/tasrkleat-TCGA-results/tasrkleat-TCGA-results-post/' )
    output=${output_dir}/genomecov.bed.csv
    # echo $output
    mkdir -p ${output_dir}
    cmd="bedtools genomecov -bga \
             -ibam ${bam} \
             -g tasrkleat-TCGA-results/tasrkleat-static/hg19.fa \
             > ${output}"
    echo ${cmd}
done


