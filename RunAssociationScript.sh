#!/bin/bash -x

# Script to loop through cytokine names to start jobs simultaneously 
# FILE = list of cytokines, one per line
 
for i in $(seq 1 1 21); 
do 
	echo $i 
	qsub -v FILE=$i StartGemmaRun.pbs
done


	

