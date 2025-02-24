tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools

reference_genome=/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps/Pseudomonas.OTU5_ref.fasta
sixgene=/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/output.fasta
paf_file=/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/Pseudomonas.OTU5_ref_sixgene.paf
$tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$sixgene"> "$paf_file"
