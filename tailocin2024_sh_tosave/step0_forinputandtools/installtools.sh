#seqtk
#just git and make
#under jiajucui/tools
cd ./tools
#SPAdes
    wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
    tar -xzf SPAdes-4.0.0-Linux.tar.gz
    cd SPAdes-4.0.0-Linux/bin/
#requires python 3.8 failed 
#try in env phylogeny_snp with python  Python 3.10.6 
conda install -c bioconda spades 
#version is 3.15.5
#  environment location: /SAN/ugi/aMetagenomics/jiajucui/tmpbigfile/miniconda3newdirtostore/envs/phylogeny_snp




#minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
#r2t env


#clustalo
wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64
#rename 
#UNIX and Mac users should rename the downloaded file to clustalo and place in the location of their choice. This file may need to be made executable e.g.: chmod u+x clustalo
