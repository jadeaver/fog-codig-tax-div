# fog-codig-tax-div

### This repo contains scripts utilized in analyses described by Deaver et al. 2021

### The data that support the findings of this study are openly available in MG-RAST at https://www.mg-rast.org/linkin.cgi?project=mgp88050, project ID mgp88050 and accession numbers mgm4829552.3 to mgm4829571.3.

### bash-code 
The code described in this file intends to pull 16S rRNA gene sequences from SILVA and RefSeq database files for desired microorganisms specified to the species level. The purpose of pulling these representative sequences is to create a pyhlogenetic tree for the bacterial and archaeal species present in ACoD samples described in this study.

### Inputs:
The Org_name.txt file contains the species-level names of microorganisms for which representative 16S rRNA gene sequences were retrieved. The database files must be downloaded separately from  https://www.arb-silva.de and https://www.ncbi.nlm.nih.gov/. 

### Tree Input: 
The final compilation of representative 16S rRNA gene sequences are provided in the tree_input.fasta file. This file was submitted for tree creation using a custom workflow within the NGPhylogeny tool (https://ngphylogeny.fr/) that utilized the tools MUSCLE, BMGE, FastTree, and Newick Display. 

### Tree Output:
The output tree for this study is provided, FastTree_output_tree.nhx.

### Renaming Tree nodes: 
Node names were renamed to otu# (otu1, otu2, . . . , otu1417). A mapping file was provided as input with the previous node name and an associated otu#.  Renaming the nodes was necessary to ensure node names matched taxa names in the phyloseq object. 

### R-code.R
This script takes the OTU/metadata/taxonomy tables and phylogenetic tree to create a phyloseq object (**phyloseq**: https://joey711.github.io/phyloseq/index.html). Weighted UniFrac distances are calculated and used to generate NMDS plots and to input into ANOSIM and Mantel tests (**vegan**:  https://rdrr.io/rforge/vegan/). Indicator species analysis is performed (**indicspecies**: https://rdrr.io/cran/indicspecies/)

### Inputs:
Located in "data" folder
unrooted-tree: FastTree_otu.nhx
otu-table: otu_T.xlsx
taxonomy-table: taxa_T.xlsx
metadata: metadata.xlsx
