#!/usr/bin/bash
#PBS -N FUSION
#PBS -l nodes=1:ppn=1,mem=15G
#PBS -l walltime=24:00:00
#PBS -p +1023
#PBS -t 1-1078%1
#PBS -j oe
#PBS -o /public/home/biostat06/tmp/FUSION/FUSION.out
#PBS -q batch

bash
let k=0
PROJ_PATH=/public/home/biostat06/TWAS/
FUSION=${PROJ_PATH}fusion_twas-master/FUSION.assoc_test.R
var_names=${PROJ_PATH}49_tissues.txt

for pheno in $(cat ${var_names})
do

mkdir ${PROJ_PATH}FUSION_RESULT/${pheno}

for chr in `seq 1 22`
do

let k=${k}+1
if [ ${k} -eq ${PBS_ARRAYID} ]
then

TIMELOG=/public/home/biostat06/tmp/FUSION/${pheno}_${chr}.txt
/usr/bin/time -v -o ${TIMELOG} \
Rscript ${FUSION} \
--sumstats /public/home/biostat06/TWAS/Preeclampsia_gwas_summary/outcome.txt \
--weights /public/home/biostat06/TWAS/WEIGHTS2/${pheno}.pos \
--weights_dir /public/home/biostat06/TWAS/WEIGHTS2/ \
--ref_ld_chr /public/home/biostat06/TWAS/LDREF/1000G.EUR. \
--chr ${chr} \
--out ${PROJ_PATH}FUSION_RESULT/${pheno}/chr${chr}.dat

fi
done
done