#!/usr/bin/bash
#PBS -N COLOC
#PBS -l nodes=1:ppn=1,mem=5G
#PBS -l walltime=30:00:00
#PBS -p +1023
#PBS -t 1-49%20
#PBS -j oe
#PBS -o /public/home/biostat06/tmp/COLOC/COLOC.out
#PBS -q batch

bash
let k=0
PROJ_PATH=/public/home/biostat06/TWAS/
COLOC=${PROJ_PATH}COLOC_RUN_CODE/01_coloc_split_gene.r
var_names=${PROJ_PATH}49_tissues_coloc.txt

for pheno in $(cat ${var_names})
do

let k=${k}+1
if [ ${k} -eq ${PBS_ARRAYID} ]
then

TIMELOG=/public/home/biostat06/tmp/COLOC/${pheno}.txt
/usr/bin/time -v -o ${TIMELOG} \
Rscript ${COLOC} --tissue ${pheno}

fi
done