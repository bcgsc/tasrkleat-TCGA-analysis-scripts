# UHRC1 
for i in $(\ls -v ./UHR/C1/tasrkleat-results/kleat/postproc-styleA-polyA-confidence/*.csv); do
    echo "python compare_to_polyA_Seq.py ${i} kleat"
done


for i in $(\ls -v ./UHR/C1/tasrkleat-results/kleat/postproc-styleB-tbr-tuning/*.csv); do
    echo "python compare_to_polyA_Seq.py ${i} kleat"
done
