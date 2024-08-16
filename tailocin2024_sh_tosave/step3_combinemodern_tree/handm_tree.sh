#!/bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=4:30:0
#$ -S /bin/bash
#$ -N m87otree
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/

# Remember first PL0042 cant use --meta, check readme
# and PL0066 and PL0222 have zero contig, need to be excluded in MSA

# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
houtput_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract"
moutput_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_modern85"
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/tools/
#vim ../../tailocin_modern85/haplotype_selected/all_
#all_HTF_samples.fasta                   all_modern76fa_markexceptlongest.fasta  all_TFA_samples.fasta                   
#(base) [jiajucui@pchuckle step3_combinemodern_tree]$ vim ../../tailocin_modern85/haplotype_selected/all_
#filtered_HTF_samples.fasta
new_tree_dir="${houtput_dir}/tree"

#bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/modern85_extractfastafromvcf.sh
#cat and remove the Ns, since we believe clustalo omega could do a good job to align them to the right place, may be Ns could mislead since the ref hyplotypes have differen length.
cat $houtput_dir/haplotype_selected/all_TFA_samples.fasta    $moutput_dir/haplotype_selected/all_TFA_samples.fasta  | sed 's/N//g'> $houtput_dir/tree/handm_all_TFA_samples.fasta
TFAfinal_fasta=$houtput_dir/tree/handm_all_TFA_samples_clustalo.fasta
$tools/clustalo -i $houtput_dir/tree/handm_all_TFA_samples.fasta -o $TFAfinal_fasta --force
#--force: Overwrites the output file if it already exists.
cat $houtput_dir/haplotype_selected/all_HTF_samples.fasta    $moutput_dir/haplotype_selected/all_HTF_samples.fasta | sed 's/N//g'> $houtput_dir/tree/handm_all_HTF_samples.fasta
HTFfinal_fasta=$houtput_dir/tree/handm_all_HTF_samples_clustalo.fasta
$tools/clustalo -i $houtput_dir/tree/handm_all_HTF_samples.fasta --force -o $HTFfinal_fasta

cat $houtput_dir/haplotype_selected/filtered_TFA_samples.fasta    $moutput_dir/haplotype_selected/filtered_TFA_samples.fasta | sed 's/N//g'> $houtput_dir/tree/handm_filtered_TFA_samples.fasta
TFAfilter_fasta=$houtput_dir/tree/handm_filtered_TFA_samples_clustalo.fasta
$tools/clustalo -i $houtput_dir/tree/handm_filtered_TFA_samples.fasta --force -o $TFAfilter_fasta

cat $houtput_dir/haplotype_selected/filtered_HTF_samples.fasta    $moutput_dir/haplotype_selected/filtered_HTF_samples.fasta | sed 's/N//g'> $houtput_dir/tree/handm_filtered_HTF_samples.fasta
HTFfilter_fasta=$houtput_dir/tree/handm_filtered_HTF_samples_clustalo.fasta
$tools/clustalo -i $houtput_dir/tree/handm_filtered_HTF_samples.fasta --force -o $HTFfilter_fasta

cd $new_tree_dir
iqtree -s $HTFfinal_fasta -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/tailocin_HTF"
iqtree -s $TFAfinal_fasta -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/tailocin_TFA"
iqtree -s $HTFfilter_fasta -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/filter_HTF"
iqtree -s $TFAfilter_fasta -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/filter_TFA"

mkdir -p $new_tree_dir/fullinfo/
cd $new_tree_dir/fullinfo/

TFAfinal_fasta2=$houtput_dir/tree/fullinfo/handm_all_TFA_samples_clustalo.fasta
HTFfinal_fasta2=$houtput_dir/tree/fullinfo/handm_all_HTF_samples_clustalo.fasta
TFAfilter_fasta2=$houtput_dir/tree/fullinfo/handm_filtered_TFA_samples_clustalo.fasta
HTFfilter_fasta2=$houtput_dir/tree/fullinfo/handm_filtered_HTF_samples_clustalo.fasta
# Step 5: Build a phylogenetic tree
bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/step3_combinemodern_tree/fullinfoMSA.sh $HTFfinal_fasta $HTFfinal_fasta2
bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/step3_combinemodern_tree/fullinfoMSA.sh $HTFfilter_fasta $HTFfilter_fasta2
bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/step3_combinemodern_tree/fullinfoMSA.sh $TFAfinal_fasta $TFAfinal_fasta2
bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/step3_combinemodern_tree/fullinfoMSA.sh $TFAfilter_fasta $TFAfilter_fasta2




iqtree -s $HTFfinal_fasta2 -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/fullinfo/ftailocin_HTF"
iqtree -s $TFAfinal_fasta2 -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/fullinfo/ftailocin_TFA"
iqtree -s $HTFfilter_fasta2 -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/fullinfo/ffilter_HTF"
iqtree -s $TFAfilter_fasta2 -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/fullinfo/ffilter_TFA"


echo "Analysis complete. Check the directories for results."

