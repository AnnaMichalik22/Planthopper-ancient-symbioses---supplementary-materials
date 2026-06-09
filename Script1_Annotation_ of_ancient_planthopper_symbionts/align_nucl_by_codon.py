#! /usr/bin/env python3

import sys, os

if len(sys.argv) != 4:
	sys.exit('------------------------ align_nucleotide_by_codon.py v.3.0     Piotr Lukasik      21st March 2024 ------------------------\n'
	         'This script translates nucleotide sequences, aligns resulting protein sequences using mafft, and reverse-translates the alignment to nucleotide space.\n'
	         'Usage: ./align_nucl_by_codon.py <input_fasta_file> <output_fasta_file> <translation_table> \n')

Script, Input_fasta, Output_fasta, Translation_table = sys.argv

Output_prot_fasta = Output_fasta[:-6] + "_prot.fasta"
Output_nucl_fasta = Output_fasta[:-6] + "_nucl.fasta"



def ImportFasta(fasta_file):
   FASTA = open(fasta_file, 'r')
   Seq_list = []
   Sequence = ''
   Seq_heading = ''
   for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
      if line.startswith('>'):
         if Sequence != '':   # Saves the existing Seq_heading and Sequence to a list before overwriting them
            Seq_list.append([Seq_heading, Sequence])
            #print(Seq_list[-1])
         Sequence = ''
         Seq_heading = line.strip().strip(">") # Takes the whole name as heading
      else:
         Sequence = Sequence + line.strip().upper()
   Seq_list.append([Seq_heading, Sequence.strip('')]) # Saves the final sequence (Seq_heading and Sequence) to a list
   #print(Seq_list[-1])
   FASTA.close()
   return(Seq_list)




os.system('transeq -frame 1 -table %s -sequence %s -outseq ./seq_temporary.fasta' % (Translation_table, Input_fasta))
os.system("mafft --maxiterate 1000 --thread 2 --quiet --localpair ./seq_temporary.fasta > %s" % Output_prot_fasta)

FASTA = ImportFasta(Input_fasta)
for sequence in FASTA:
    no_gaps = ""
    for nucl in sequence[1]:
        if nucl in ['A', 'C', 'G', 'T', 'N']:
            no_gaps += nucl
    sequence[1] = no_gaps

PROTFASTA = ImportFasta(Output_prot_fasta)

for i in range(len(FASTA)):
    FASTA[i].append(PROTFASTA[i][1])

OUTPUT = open(Output_nucl_fasta, 'w')

for sequence in FASTA:

    #print("\n>%s" % sequence[0])
    nucl_pos = 0
    nucl_aligned = ""
    for k in range(len(sequence[2])):
        if sequence[2][k] == "-":
            nucl_aligned += "---"
            #print("-\t---")
        else:
            nucl_aligned += sequence[1][nucl_pos:nucl_pos+3]
            nucl_pos += 3
            #print("%s\t%s" % (sequence[2][k], sequence[1][nucl_pos:nucl_pos+3]))
            
    print(">%s" % sequence[0], file=OUTPUT)
    print(nucl_aligned, file=OUTPUT)
    

OUTPUT.close()

os.system("rm ./seq_temporary.fasta")

print("Job finished. Successfully completed alignment of %d sequences." % len(FASTA))
