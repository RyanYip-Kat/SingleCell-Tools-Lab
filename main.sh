#python ldsc/munge_sumstats.py --sumstats  ldsc-seg/body_BMIz.sumstats.gz --merge-alleles ldsc-seg/w_hm3.snplist  --out UKBB_BMI
cts_name=body
python ldsc/ldsc.py \
	--h2-cts ldsc-seg/UKBB_BMI.sumstats.gz \
	--ref-ld-chr ldsc-seg/1000G_EUR_Phase3_baseline/baseline. \
	--out BMI_${cts_name} \
	--ref-ld-chr-cts ldsc-seg/Corces_ATAC.ldcts \
	--w-ld-chr ldsc-seg/weights_hm3_no_hla/weights.
