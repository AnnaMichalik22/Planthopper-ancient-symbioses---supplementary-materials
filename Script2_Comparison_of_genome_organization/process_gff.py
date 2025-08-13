#! /usr/bin/env python3

import sys, re

if len(sys.argv) != 3:
	sys.exit('This script processes a gff file - like the one produced by symCAP script -\n'
	         'and the second tab-delimited file with info on assignment of genes to pre-defned categories,\n'
	         'and returns to stdout a table to be used as an input of genome_plotting.\n'
	         'It direct the output to STDOUT - you probably want to redirect it to a file instead.\n'
	         'Usage: ./process_gff_20250307.py <gff_file> <categories_file>\n')
Script, Gff_file, Categories_tsv_file = sys.argv


def read_tsv_to_dict(filename):
    # Reads a tab-delimited two-column file into a dictionary.
    # Returns: A dictionary containing the data from the file, or None if an error occurs.
    result_dict = {}
    try:
        with open(filename, 'r', encoding='utf-8') as infile:
            for line in infile:
                line = line.strip()
                if line: # check for blank lines.
                    key, value = line.split('\t', 1) #split only once
                    result_dict[key] = value
        return result_dict
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except ValueError:
        print(f"Error: Invalid format in file '{filename}'.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None

def extract_gene_name(gene_string):
    # Extracts the gene name from a string containing gene information.
    # Args: A string such as "ID=gene_001;locus_tag=SMDICMUL_001;gbkey=Gene;gene_biotype=protein_coding;gene=lipB;name=lipB;"
    # Returns: The gene name if found, or 'unnamed_gene' if not found.

    match = re.search(r"gene=([^;]+)", gene_string)
    if match:
        return match.group(1)
        
    # If "gene=" tag is not found, try the "name=" tag
    match = re.search(r"name=([^;]+)", gene_string)
    if match:
        return match.group(1)

    return 'unnamed_gene'  # Gene name not found



Genes_categorized = read_tsv_to_dict(Categories_tsv_file)

GFF = open(Gff_file, 'r')
"""
Example input:
##gff-version 3
##sequence-region  SMDICMUL 1 142701
SMDICMUL	HMMER_3.1b2	gene	1	645	.	+	0	ID=gene_001;locus_tag=SMDICMUL_001;gbkey=Gene;gene_biotype=protein_coding;gene=lipB;name=lipB;
SMDICMUL	HMMER_3.1b2	CDS	1	645	.	+	0	ID=product_001;Parent=gene_001;locus_tag=SMDICMUL_001;gbkey=CDS;gene=lipB;product=Octanoyltransferase;transl_table=11;codon_start=1
SMDICMUL	HMMER_3.1b2	gene	674	3217	.	+	0	ID=gene_002;locus_tag=SMDICMUL_002;gbkey=Gene;gene_biotype=protein_coding;gene=valS;name=valS;
SMDICMUL	HMMER_3.1b2	CDS	674	3217	.	+	0	ID=product_002;Parent=gene_002;locus_tag=SMDICMUL_002;gbkey=CDS;gene=valS;product=Valine--tRNA ligase;transl_table=11;codon_start=1
SMDICMUL	rnammer_1.2	gene	62387	62491	.	-	0	ID=gene_072;locus_tag=SMDICMUL_072;gbkey=Gene;gene_biotype=rRNA;gene=gene;name=gene;
SMDICMUL	rnammer_1.2	rRNA	62387	62491	.	-	0	ID=product_072;Parent=gene_072;locus_tag=SMDICMUL_072;gbkey=rRNA;product=product;note=inference:RNAmmer 1.2
"""


Genome_list = []
Gene_table = []
for line in GFF:   
   Line = line.split()
   
   if Line[0] == "##sequence-region":
      Genome_list.append(Line[1])
      Gene_table.append(["genome_entry",Line[1],Line[3]])
   elif Line[0] in Genome_list:
      gene_name = extract_gene_name(Line[8])
      
      if (Line[2] == "CDS"):
         Gene_category = 'Other'
         if gene_name in Genes_categorized:
             Gene_table.append([Line[0], gene_name, Genes_categorized[gene_name], Line[3], Line[4], Line[6]])
             #['SMDICMUL', 'valS', 'genetic_information_processing', 674, 3217, '+']
      
      elif (Line[2] == "rRNA"):
         Gene_table.append([Line[0], gene_name, 'rRNA', Line[3], Line[4], Line[6]])
         
for entry in Gene_table:
   for item in entry:
      print(item, ",", sep='',end='') 
   print('\n',end='') 
   
