This repository contains custom scripts used for the analyses and results visualization presented in the manuscript:
MichalikA., Castillo Franco D., Deng J., Prus-Frankowska M., Stroiñski A., £ukasik P. (preprint-2025): "The tiniest genomes shrink much further: extreme reductive evolution in planthopper symbionts" bioRxiv https://doi.org/10.1101/2025.05.16.654412
Specifically, the repository contains:
1. A script for the annotation of fast-evolving ancient planthopper bacterial symbionts Sulcia and Vidania (SymCap.py)
This Python3 script uses a set of curated protein and rRNA/tmRNA/ncRNA reference alignments to conduct HMMer searches against a set of contigs representing symbiont genomes, and outputs annotations as GFF files, gene contents tables, and gene alignments.
Dependencies:
- HMMER 3.2
- Getorf
- rnammer
- tRNAscan-SE
- mafft
- align_nucl_by_codon.py

Inputs:
- A folder with fasta files of genomic contigs for annotation, one sequence per file;
- A folder with amino acid alignments of all protein-coding genes to be searched;
- A folder with nucleotide alignments of all all rRNA, tmRNA and ncRNA genes to be searched;
- Table with information about genes to be searched, such as name, EC number, categorization, etc.

Before running the script, edit lines 31-61 to provide correct paths to the above folders/files, and additional parameters such as the Translation Table ( #4: Spiroplasma, Hodgkinia, Nasuia; # 11: Sulcia, Vidania, most bacteria; #  5: Invertebrate mitogenome).
All files used for the final annotations of Sulcia and Vidania in the manuscript are available in the folder SymCap_files
Run the script: ./SymCap.py


2. A Processing script for the visualization of the symbiont genome contents table (Genome_content_visualisation.pyde)
This Processing v. 3.5.4 (Python mode) script was developed for the visualization of the results of gene contents comparisons using the SymCap script. It inputs a table showing the presence or absence of functional genes, pseudogenes, and putative pseudogenes in each of the genomes, and draws a figure where these genes are organized in a tabular form, with a different color used for each category.
As input for the Genome_content_visualisation script, you need to provide the output file from the SymCap script: Genome_content_table.csv. You may choose to sort or trim the table, and the script can be modified easily to enable the visualization of similar data tables.

3. A script for visualizing the contribution of symbionts to amino acid biosynthesis pathways (Nutritional_complementarity.pyde)
This Processing v3.5.4 (Python mode) script was developed to visualize symbiont contributions to amino acid biosynthesis pathways. It generates a grid-style figure with color-coded pie charts representing the symbionts involved in amino acid biosynthesis.
Before running the script, you need to provide .txt files containing:
- A list of species
- A list of symbionts
- A list of genes
- A list of colors for symbiont identification
- Merged data from all sources
Examples of all required files are available in the folder Nutritional_complementarity_files.

4. A script for comparison of genome organization (Promer_style_genome_comparison.pyde)
This Processing v3.5.4 (Python mode) script was developed for comparing the genome organization of Sulcia and Vidania symbionts.
Before running the script, you need to provide:
- A list of genes and gene categories
- GFF files of annotated genomes
- A GFF file of the reference genome
The script uses annotated genomes in .pro format as input. To convert GFF files into .pro files, use the process_gff.py script:
Run the script: ./process_gff.py
After that, you can run the Promer_style_genome_comparison.pyde script. As input, you need to provide:
- .pro file containing the target genomes
-  .pro file of the reference genome



