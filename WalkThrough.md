# Genome annotation pipeline
This pipeline trains and runs Augustus gene prediction for a new eukaryotic species.\
The codes are intended to be running on NCI Gadi, it will need adjusting if running on other servers.
### Required inputs:
- Trinity assembly assembled from RNA-seq of your species (fasta file)
- Soft-masked genome assembly of your species (fasta file)
### Outputs:
- Gene annotation (gff3)
- Predicted peptides (fasta)
- Optimised Augustus species-specific parameters for your species
## 1. Setting up the working directory
**Step 1.**  Create a working directory
```
export workingdir=/path/to/working_directory
mkdir ${workingdir}
```
**Step 2.** Download the latest Uniprot database fasta file, and the Augustus working directory
```
cd ${workingdir}
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
git clone
```
