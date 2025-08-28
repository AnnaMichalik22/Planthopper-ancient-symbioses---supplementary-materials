# 🧬 Symbiont Genome Analysis Scripts

Custom scripts for analyses and visualization used in the manuscript:  

**Michalik A., Castillo Franco D., Deng J., Prus-Frankowska M., Stroiñski A., Łukasik P. (preprint-2025):**  
*The tiniest genomes shrink much further: extreme reductive evolution in planthopper symbionts*  
[bioRxiv, 2025.05.16.654412](https://doi.org/10.1101/2025.05.16.654412)  

---

## 📂 Repository Overview

This repository provides scripts to annotate, compare, and visualize the genomes of ancient bacterial symbionts (*Sulcia* and *Vidania*).

**Included scripts:**
1. **SymCap.py** – Genome annotation of fast-evolving symbionts  
2. **Promer_style_genome_comparison.pyde** – Genome organization comparison  
3. **Nutritional_complementarity.pyde** – Visualization of symbiont contributions to amino acid biosynthesis  
4. **Genome_content_visualisation.pyde** – Visualization of genome content tables  

---

## ⚙️ Dependencies

Make sure the following are installed before running the scripts:

- [HMMER 3.2](http://hmmer.org/)  
- [Getorf](https://www.ebi.ac.uk/Tools/sfc/getorf/)  
- [RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/)  
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)  
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)  
- `align_nucl_by_codon.py` (included)  
- [Processing v3.5.4](https://processing.org/) (Python mode)  

---

## 📜 Scripts

### 1. **Genome Annotation – `SymCap.py`**
Annotates symbiont genomes (*Sulcia* and *Vidania*) using curated protein and RNA reference alignments.  

**Inputs:**
- FASTA files of contigs (one sequence per file)  
- Amino acid alignments (protein-coding genes)  
- Nucleotide alignments (rRNA, tmRNA, ncRNA)  
- Gene metadata table (name, EC number, category, etc.)  

**Usage:**
```bash
# Edit lines 31–61 in SymCap.py to configure paths and parameters
./SymCap.py
**Output**:

- GFF annotation files
- Gene content tables
- Gene alignments

👉 Example files are in SymCap_files/

### 2. **Genome Comparison – `Promer_style_genome_comparison.pyde`**
A **Processing v3.5.4 (Python mode)** script for comparing genome organization of *Sulcia* and *Vidania* symbionts.  
This was developed after we found that **PROMER** (a standard comparison tool) struggled to detect extremely divergent but homologous regions.

**Required inputs:**
- Functional gene list (grouped into categories)  
- Gene positions in the query genome (`.pro` format, generated via `process_gff.py`)  
- GFF files of annotated genomes  

**Convert GFF → PRO:**
```bash
./process_gff.py [gff_file] > [pro_file]
# Example:
./process_gff.py VFDICMUL.gff > VFDICMUL.pro

### 3. **Nutritional Complementarity – `Nutritional_complementarity.pyde`**
A **Processing v3.5.4 (Python mode)** script developed to visualize symbiont contributions to amino acid biosynthesis pathways.  

**Output:**  
- A **grid-style figure** where:  
  - Each row = a species  
  - Each column = a gene  
  - Each cell = a pie chart showing the symbionts encoding that gene  
  - Pie chart size = proportional to the number of symbionts containing the gene  

**Required inputs (.txt files):**
- Species list  
- Symbiont list  
- Gene list  
- Color scheme for symbiont identification  
- Merged dataset (list of all possible metagenome–gene combinations with info on the symbiont genome where the gene is present)  

👉 Example input files are available in `Nutritional_complementarity_files/`.  

### 4. **Genome Content Visualization – `Genome_content_visualisation.pyde`**
A **Processing v3.5.4 (Python mode)** script for visualizing the results of gene content comparisons produced by the **SymCap** script.  

**Input:**
- `Genome_content_table.csv` (output from `SymCap.py`)  

**Features:**
- Generates a **tabular figure** showing presence/absence of:  
  - Functional genes  
  - Pseudogenes  
  - Putative pseudogenes  
- Uses distinct colors for each gene category  
- Can be easily customized to sort/trim the table or adapted to visualize similar datasets  

As input for the Genome_content_visualisation script, you need to provide the output file from the SymCap script: Genome_content_table.csv. 
You may choose to sort or trim the table, and the script can be modified easily to enable the visualization of similar data tables.
