#!/bin/bash -l

#$ -l tmem=10G
#$ -l h_vmem=10G
#$ -l h_rt=48:00:0
#$ -wd /SAN/ugi/plant_genom/jiajucui/
#$ -V
#$ -N run46
#$ -t 1-46
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/

#variables
i=$SGE_TASK_ID
samplename=$(cat /SAN/ugi/plant_genom/jiajucui/phylogeny/samplelist/tailocin46.txt | sed -n $i'p')
refPswithtailocin=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta
#/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/tailocin_46/
source ~/miniconda3/bin/activate phylogeny_snp 

##use bwa to realign the unmapped bam file
bwa aln -t 2 -l 1024 ${refPswithtailocin} -b /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2024_413/${samplename}_after_removal_mappedAt.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_realign.sai
#_after_removal_mappedAt.bam
bwa samse -r @RG\\tID:${samplename}\\tSM:${samplename} -f /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.sam ${refPswithtailocin} /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_realign.sai /SAN/ugi/plant_genom/jiajucui/4_mapping_to_A_thaliana/2024_413/${samplename}_after_removal_mappedAt.bam
#Keep the mapped reads and create a compressed BAM file using samtools
samtools view -@ 2 -F 4 -Sbh -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.sam

#Sort the BAM file by chromosome and position using samtools
samtools sort -@ 2 -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.sorted.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.bam

samtools view -q 20 /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.sorted.bam -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}.mapped_to_Pseudomonas.dd.q20.bam


#rm sai sam bam and unsort unq20 bams:
#keep for other refs
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_realign.sai
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.sam
rm /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_removalAt_mapped_to_ps.sorted.bam

samtools markdup /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}.mapped_to_Pseudomonas.dd.q20.bam /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam

samtools flagstat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}.mapped_to_Pseudomonas.dd.q20.markeddup.bam > /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46/${samplename}_Ps_q20_flagstats.log
