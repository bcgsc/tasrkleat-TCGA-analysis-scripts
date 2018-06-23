# UHRC1 
for i in $(\ls -v ./UHR/C1/tasrkleat-results/kleat/postproc-styleA-polyA-confidence/*.csv); do
    echo "python compare_to_polyA_Seq.py ${i} kleat"
done


for i in $(\ls -v ./UHR/C1/tasrkleat-results/kleat/postproc-styleB-tbr-tuning/*.csv); do
    echo "python compare_to_polyA_Seq.py ${i} kleat"
done


# concatenate comparison results
# cat ./UHR/C1/tasrkleat-results/kleat/postproc-styleA-polyA-confidence/*vs-polyA-Seq.csv  > UHR/C1/tasrkleat-results/kleat/styleA-filtered-benchmark.txt
# cat ./UHR/C1/tasrkleat-results/kleat/postproc-styleB-tbr-tuning/*vs-polyA-Seq.csv  > UHR/C1/tasrkleat-results/kleat/styleB-filtered-benchmark.txt

