# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/tools/
omega_55tapemeasure=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/omega55.fasta
fasta55=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/alltapemeasureforfish.fasta
$tools/clustalo -i $fasta55 --force -o $omega_55tapemeasure
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/tree
cd /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/tree
new_tree_dir=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/tree
iqtree -s $omega_55tapemeasure -nt AUTO -bb 1000 -alrt 1000 -pre "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/tree/tape55"



#full info tree 
mkdir -p $new_tree_dir/fullinfo/
cd $new_tree_dir/fullinfo/
fullinfo55=$new_tree_dir/fullinfo/fullinfo_omega55.fasta
bash /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/step3_combinemodern_tree/fullinfoMSA.sh $omega_55tapemeasure $fullinfo55

iqtree -s $fullinfo55 -nt AUTO -bb 1000 -alrt 1000 -pre "/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure/step3/tree/fullinfo/fulltape55"


            
