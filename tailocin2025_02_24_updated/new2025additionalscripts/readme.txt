Details of the strategy please refer to the notion note 'week 26 redo tailocin'

all folders :tailocin_46  tailocin_46_fastq  tailocin_extract  tailocin_modern85 shfortailocin


in shfortailocin:

./forinputandtools:
2024_tailocin46.sh  installtools.sh  lnfastq.sh

./logs:

./others/other_tailocininfo:
formatted_regions.sh  formatted_regions.txt  otherbams  tailocin_genes.txt

./step1_extractreadsfrom46h:
extractfastqR1R2for10nsqs.sh	pipetoextractreadsandcutfastq.sh	readme

./step2_assembly:

pipeofassemblyandMSAtree.sh  

step3_combinemodern_tree/:
cumulative_positions.sh   fullinfoMSA.sh  mergehandm_calculatecoveredlist_maskexceptlongest.sh  modern85_extractfasta.sh  readme
cumulative_positions.tsv  handm_tree.sh   mergehandm_calculatecoveredlist.sh                    oldvcf

./tools:
clustalo  minimap2

prepare inputs:
1. get the reads with new ref:
/SAN/ugi/plant_genom/jiajucui/datapreprocessing/sh_scripts/2024_225/4thfastq/2024_tailocin46.sh
or cp here 2024_tailocin46.sh

2. installtools.sh
then first ln -s the original collapsed fastq of the 46 samples and rename the 109_* to  109.*, withlnfastq.sh

3. use formatted_regions.sh to convert tailocin_genes.txt to formatted_regions.txt 
/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/others/other_tailocininfo/formatted_regions.sh


step1: extract reads and generate fastq:
run the pipetoextractreadsandcutfastq.sh, the results are in trimmed_fastq folder 

step2:  generate assemblies:
pipeofassemblyandMSAtree.sh
give assemblies, contig_count and msa and also tree, using spades, with two method for 143 and for HBPL, using minimap to generate pal and build MSA with clustalo omega.


step3: combine modern and tree building: 
#
##bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/modern85_extractfasta.sh
#then build the tree qsub handm_tree.sh 

Step4: tape measure
step1_extracttapemeasure.shstep2_retrievetape_formodernNs.shstep3_1_2024_tailocin46tape.shstep3_1_generate55tapeMSA.shstep3_2_extractreads.shstep3_2_extractreadsfor10R1R2.shstep3_3_fish_historical_all.sh

calculatelength.sh #calculate the length of 55 tape measure (they have 3 different lengths)treebuildtape.sh #build the trees (omega and omega&fullinfo (no Ns))
fullinfoMSA.sh #to build the full info MSA for tree building





caution:
besides PL0066 and PL0222, who have 0 contig, HB0863 and 27. have no mapped reads after minimap2, need to be removed in msa building.

evidence:
'(r2t) [jiajucui@pchuckle mappings]$ ls -lh | grep ' 24 J'
-rw-r--r-- 1 jiajucui miscbio   24 Jun 30 00:26 27.ESP_1975_mapped_sorted.bam.bai
-rw-r--r-- 1 jiajucui miscbio   24 Jun 30 00:26 HB0863_mapped_sorted.bam.bai
-rw-r--r-- 1 jiajucui miscbio   24 Jun 30 00:26 PL0066_mapped_sorted.bam.bai
-rw-r--r-- 1 jiajucui miscbio   24 Jun 30 00:26 PL0222_mapped_sorted.bam.bai
(r2t) [jiajucui@pchuckle mappings]$ less HB0863_mapped
HB0863_mapped.bam             HB0863_mapped_sorted.bam      
HB0863_mapped.sam             HB0863_mapped_sorted.bam.bai  
(r2t) [jiajucui@pchuckle mappings]$ less HB0863_mapped_sorted.bam
"HB0863_mapped_sorted.bam" may be a binary file.  See it anyway? 
(r2t) [jiajucui@pchuckle mappings]$ samtools view HB0863_mapped_sorted.bam
NODE_1_length_576_cov_3.792706	4	*	0	0	*	*	0	0	CTCATGGCCCTGCTCGGTCAGGATCGACTCGGGGTGAAACTGTACCCCTTCGACGTTCAGTGTCTTGTGACGCAGGCCCATGATCTCGTCAACCGAGCCGTCGTCATGTGCTGTCCATGCAGTCACTTCCAGAGCGTCAGGCAGTGTATCCAGCTTGACCACCAGCGAGTGATAACGCGTGACGACCAGCGGATGGTTCAAGCCTTCGAACACGCCCTGGTCTTCATGGATCACCGGACTGGTCTTGCCGTGCATGACCTGACGCGCACGGACCACGTCACCACCAAATGCCTGACCAATGGACTGATGCCCGAGGCAGACACCCAGAATCGGCAGTTTGCCCGCGAAGTGATTGATCACTTCCAGCGAGACGCCAGCTTCATTGGCAGTGCAAGGACCGGGCGAAACAACGATGCGCTCGGGGTTCAGGGCTTCGATCTGCGCGATGGTCAGCTCGTCATTGCGGATGACCTTCACGTCGGCGCCCAGTTCGCCCAGGTACTGCACGACGTTGTAGGTGAAGGAATCATAGTTATCGATCATCAGCAGCATTTTGCTATTAACCTTTTGAATTTA	*	rl:i:0
NODE_2_length_329_cov_2.905109	4	*	0	0	*	*	0	0	TTCAAGGCCGTACACACGCCAGACATGCCGCTCAAACTGCTCACTACCGGTGTGGATATGACCATCACCTGCGAATTCAACAATGGCAAGACCTACGTTCTGTCCGGCGCCTATCTGGTCGAAGAACCTGTCAGCAAGGCCGATGACGGCACCATCGACTTGCGATTCGATGGCAGTCAGGGGAGCTGGCAATGAGCGAAGTCATCGAACTGGTGAGCCCGATCGAGGCTCATGGTGAAACCGTCTCGCAACTGACGTTTCGCCGCCCCACCGCACAGGAAGCGCGCGCCATCAAGGCCCTGCCTTACCGGATCGACAAAAACGAGGAC	*	rl:i:0
NODE_3_length_319_cov_2.621212	4	*	0	0	*	*	0	0	ATTCGATCAACGAGAAGTGGGCGCCGATCAAAACCTTCTTCTCGGGTCTGTTCACGGACATCAAGCCGTTCATTGACCCGATCCTGAACTGGTTCGGTATCGGCTCCGATGACAAGTCTCTGCTGCAAAAAGCCACTCAAAAGCTCAACGAGTTCGCTGAAGAACGGCGTGTGGATAACACCGGTCCCGGTGGTGGCAAAGGGGCTTTTCTGACCGCCGACGCTGTGCAGGTCAGCCAGTTGAAACAGCAGCAGATCAATCAGGCAATGGGCATACCGGCAACCAGCCAGTTGCTGAGCGCGCCTAACCTGCCGGCCCC	*	rl:i:0
(r2t) [jiajucui@pchuckle mappings]$ samtools view 27.ESP_1975_mapped_sorted.bam
NODE_1_length_246_cov_1.020942	4	*	0	0	*	*	0	0	CAAACCCGGTCAGCTCGGCGCGTGGCACATTGGCCTTGGCGCCTTCAGTGACAATCGCCGCAAGACCACTGGCGCTTTCCGGAAGCCGTTCACTCAGGTCAAGAATGTCACGGCTCATTGCCTGGAATTGCTGCGGGGTTTCAAACGTCACCGACCGCTTCACACCGGCCATGCTTGTCTCGAAACCCATCGCGGCCTTGACCGCATCGATCAGCGGCTACGCCAGGGTGCTTTCGCCGATCATTT	*	rl:i:0
NODE_2_length_207_cov_1.118421	4	*	0	0	*	*	0	0	TTCGATCACATCGTCGGCCACCGGTACATCGGCACGTTGCAGCGGTTGAACCCGTTGCTGATCCAACGAAGGAACGACCAGAATCGGTGTCGTGCTGACCGCCACCGGCATACTGGCCACAACCCTCGCCACCTTGACGAGCAGCGCGTCCTGCACCAGATCCGCCGTGGCCTGAGCAGTCACCCCGGTATCGAGACCTCCGCCCTG	*	rl:i:0
(r2t) [jiajucui@pchuckle mappings]$ samtools view PL0066_mapped_sorted.bam
(r2t) [jiajucui@pchuckle mappings]$ 
'
flag 4 means unmapped


#then found ERROR: Sequence 75.LTU_1894_S30 contains not enough characters (18045)
#ERROR: Sequence HB0814 contains not enough characters (18055), the length of ref is 18057, check the bam reads and found HB0814 'NODE_16_length_528_cov_4.803383 16      tailocin        1815    60      6S295M1D227M ' and 'NODE_4_length_1175_cov_4.098214 0       tailocin        15610   60      75M2D1062M38S', there are in total three Ds, deletion.

samtools mpileup -uf $tailocin_fasta $sorted_bam | bcftools call -c --skip-variants indels --ploidy 1 -Oz -o  $vcf_out                                  

#need -V indels/ (bcftools call -c --skip-variants indels  ) to remove indels or HB0814 and 75. will be short than 18057 since there are some D in the bam reads

