#!/bin/bash -x


for i in $(seq 1 1 104); 
do 
	
	GENE_LIST="Short_gene_list_${i}.txt"
	GENE_LIST_DIR=/group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_inputs/short_gene_lists
	echo $GENE_LIST
	qsub -v FILE=$GENE_LIST_DIR/$GENE_LIST /group/nicolae-lab/users/mstein3/FcReceptor_June2017/gemma_DE_scripts/Run_Gemma_Fc_DE.pbs 
	 
done

