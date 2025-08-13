#! /usr/bin/env python3

#### symCAP_V.py version 20240622
### Developed by Piotr Łukasik, with modifications by Diego Castillo France.
### Last updated: 22nd June 2024
### Before running, ensure that all paths specified in the top 60ish lines are correct!


### This script annotates all genomes in the specified folder, using:
### 1) a list of aligned protein sequences of all query genes
### 2) a list of aligned sequences of rRNA, tmRNA and ncRNA genes
### It classifies all genes based on relative length to the query as "Functional", "Putative" or "Pseudogene".
### It outputs GFF files, alignments of all genes classified as "Functional", and a table with information on classification of each gene in each of the genomes.


### Dependencies:
    ### EMBOSS tools - getorf
    ### HMMER
    ### rnammer
    ### tRNAscan-SE
    ### mafft
    ### Piotr's script "align_nucl_by_codon.py"

import datetime, time, sys, re, glob, os
from operator import itemgetter
start_time = time.time()


# ------- Parameters to be specified before each run ------- #

  ### Where output folders etc. will be created
work_dir = "/home/anna.michalik/Sulcia_annotation/final_annotation/"

  ### Folder with genomes to annotate. Any fasta files in that dir will be regarded as genomes for annotation
genomes_for_annotation_dir = "/home/anna.michalik/Sulcia_annotation/final_genomes/"

  ### Folder containing alignments of all protein-coding genes to be annotated. File names: gene_prot.fasta
#protein_ref_dir = "/Users/Piotr/Documents/01_Cicadas/Annotation/annotation_references/protein_alignments_allcics/"
protein_ref_dir = "/home/anna.michalik/Sulcia_annotation/gene_references/proteins24/"

  ### Folder containing alignments of all rRNA, tmRNA and ncRNA genes to be searched. File names: gene.fasta
rRNA_ref_dir = "/home/anna.michalik/Sulcia_annotation/gene_references/RNAs/"

  ### Table with detailed information about every gene: name, EC number, etc.
gene_detail_table = "/home/anna.michalik/Sulcia_annotation/Sulcia_gene_list_final.txt"

evalue_for_hmmsearch = "1e-5"  ### Used for searching ORF databases using hmmsearch [Block 4]
Threshold_functional = 0.6
Threshold_putative = 0.9
Upper_threshold_functional = 1.5
Size_reference = "RANSCY"
Translation_table = 11	### As defined by getorf / NCBI!
						###  4: Spiroplasma, Hodgkinia, Nasuia; 
						### 11: Sulcia, Vidania, most bacteria;
						###  4: Mitogenome


#### These are the cutoffs for the reporting of regions without annotation, and ORFs within these gaps
gap_length_nucleotide_threshold = 300  ### Gaps in annotation above this limit will be reported
ORF_min_length = 120     ### Only ORFs above this length limit will be reported
ORF_overlap_limit = 60  ### Number of nucleotides that the ORF is allowed to overlap with adjacent annotations



#Blocks_to_run = ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B11", "B12", "B15", "B16", "B18", "B19", "B20", "B21"]
Blocks_to_run = ["B1", "B2", "B3", "B4", "B5", "B6", "B11", "B12","B15", "B16", "B18", "B19", "B20", "B21"]



 # Block 1. List genomes, save genome names to genome_list
 # Block 2. For each genome in genome_list, do six-frame translation and export all ORFs to work_dir/genomes/GenomeX/GenomeX_orfs.fasta
 # Block 3. Read protein list of Reference1, produce list of lists - Prot_seq_list
 # Block 4. Export all proteins from Prot_seq_list, search reference2, create HMM
 # Block 5. For each gene, for each genome, open hmm generated in the previous step, parse it, classify genes, updates start codons, reclassify, and and add entries to Prot_seq_list
 # Block 6. Print annotation summary for each gene
 # Block 7. Saves the Prot_seq_list table to file -- 'Protein_table.sav'
 # Block 8. Reopens 'Protein_table', reading the output as Prot_seq_list
 # Block 9.  -----------
 # Block 10. -----------
 # Block 11. Searches genomes for rRNA genes, generates list "rRNA_gene_list"
 # Block 12. Searches genomes for tRNA genes, generates list "tRNA_gene_list"
 # Block 13. -----------
 # Block 14. -----------
 # Block 15. Exports Prot_seq_list, tRNA_gene_list and rRNA_gene_list as table that can be illustrated in Processing
 # Block 16. Exports and aligns all functional gene copies
 # Block 17. -----------
 # Block 18. For each genome, combines "Prot_seq_list", "rRNA_gene_list", "tRNA_gene_list" - and sorts them, assigning IDs
 # Block 19. For all CDS, uses the details of TETULN annotation to complete , "rRNA_gene_list", "tRNA_gene_list" - and sorts them, assigning IDs
 # Block 20. Exports GFF
 # Block 21. Checks genomes for long sequence stretches without annotations


 

# ------- Modules, definitions ------- #
if Translation_table == 5:  ### Mitogenome
    start_codon_list = ['ATA', 'ATT', 'ATC', 'ATG', 'GTG', 'TTG']
    stop_codon_list = ['TAA', 'TAG']
if Translation_table == 11:   ### Sulcia, bacteria
    start_codon_list = ['ATA', 'ATT', 'ATC', 'ATG', 'GTG', 'TTG', 'CTG']
    stop_codon_list = ['TAA', 'TAG', 'TGA']
elif Translation_table == 4:   ### Spiroplasma, Hodgkinia, Nasuia
    start_codon_list = ['ATA', 'ATT', 'ATG', 'ATC', 'GTG', 'TTG', 'CTG']
    stop_codon_list = ['TAA', 'TAG']

Compl = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N':'N'}
min_rRNA_gene_length = {'23s_rRNA' : 1800, '16s_rRNA': 1000, '5s_rRNA': 80, 'rnpB': 180, 'ssrA': 250} ### prevents erroneous rRNA annotation when hit is short
max_rRNA_gene_length = {'23s_rRNA' : 3100, '16s_rRNA': 1600, '5s_rRNA': 120} ### prevents erroneous superlong rRNA annotation, which happened in one case



		# Imports a multifasta as a list of lists, where secondary lists consist of two
		# elements: heading and sequence. Important: as a heading, takes only the second part of the seq name is used
def ImportFastaRefs(fasta_file):
   FASTA = open(fasta_file, 'r')
   Seq_list = []
   Sequence = ''
   Seq_heading = ''
   for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
      if line.startswith('>'):
         if Sequence != '':   # Saves the existing Seq_heading and Sequence to a list before overwriting them
            Seq_list.append([Seq_heading, Sequence])
         Sequence = ''
         Seq_heading = line.strip().split()[1] # As a heading, takes only the second part of the seq name
      else:
         Sequence = Sequence + line.strip().upper()
   Seq_list.append([Seq_heading, Sequence]) # Saves the final sequence (Seq_heading and Sequence) to a list
   FASTA.close()
   return(Seq_list)


		# Imports a multifasta as a list of lists, where secondary lists consist of two
		# elements: heading and sequence. Takes the whole sequence name as heading
def ImportFasta(fasta_file):
   FASTA = open(fasta_file, 'r')
   Seq_list = []
   Sequence = ''
   Seq_heading = ''
   for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
      if line.startswith('>'):
         if Sequence != '':   # Saves the existing Seq_heading and Sequence to a list before overwriting them
            Seq_list.append([Seq_heading, Sequence])
         Sequence = ''
         Seq_heading = line.strip().strip(">") # Takes the whole name as heading
      else:
         Sequence = Sequence + line.strip().upper()
   Seq_list.append([Seq_heading, Sequence.strip('')]) # Saves the final sequence (Seq_heading and Sequence) to a list
   FASTA.close()
   return(Seq_list)





# ------- Specifying variables essential for further steps ------- #
Prot_seq_list = []
rRNA_gene_list = []
tRNA_gene_list = []






# Block 1. List genomes, save genome names to genome_list
#### This block lists all fasta file in directory "genomes_for_annotation_dir" provided at the beginning of script,
#### and creates a list "genome_list" : ["Genome1", "Genome2", ....]
#### Block 1 is essential for most of consecutive steps.

if "B1" in Blocks_to_run:
    print("\n\n######################## Executing Block 1 #########################################")
    ### Verifying if work directory specified exists, and if not, creating it!
    if work_dir.endswith("/"):
        work_dir = work_dir[:-1]
    if not glob.glob("%s" % work_dir):
        os.system("mkdir %s" % work_dir)

    if not glob.glob("%s" % genomes_for_annotation_dir):
        sys.exit("ERROR! Exiting ... path to directory with genomes for annotation is incorrect!")
        
    os.chdir(genomes_for_annotation_dir)
    print("Listing genomes in the provided directory...\n\nThe following genomes will be annotated: ")
    genome_list = []
    for file in glob.glob("./*.fasta"):
        print("    " + file.split('/')[1])
        genome_list.append(file.split('/')[1].split(".")[0])
    if genome_list == []:
        print("No genomes found. Aborting.")
        
    print("")
    os.chdir(work_dir)
    
    print("######################## Block 1 executed successfully! ############################\n\n")




# Block 2. For each genome in genome_list, do six-frame translation and export all ORFs to work_dir/genomes/GenomeX/GenomeX_orfs.fasta
#### Uses EMBOSS's getorf
#### Exports all ORFs to work_dir/genomes/GenomeX/GenomeX_orfs.fasta
#### Settings for Hodgkinia analysis -> min peptide length = 10 aa; translation_table = 4. 
#### Exports as ORFs any regions between STOP codons; START codons do not matter at this stage!
# Requires: genome_list

if "B2" in Blocks_to_run:
    print("\n\n######################## Executing Block 2 #########################################")
    
    minsize = 30					### getorf parameter; minimum length of ORF (nucleotides)  
    table = Translation_table		### getorf parameter; translation table
    
    print("For all genomes, doing six-frame translation (table %s) and exporting orfs above %d nucleotides from all genomes.... " % (table, minsize))
    
    if not glob.glob("%s/annotation/" % work_dir):
        os.system("mkdir %s/annotation/" % work_dir)
    if not glob.glob("%s/annotation/files" % work_dir):
        os.system("mkdir %s/annotation/files" % work_dir)
    if not glob.glob("%s/annotation/files/genomes" % work_dir):
        os.system("mkdir %s/annotation/files/genomes" % work_dir)
    for genome in genome_list:
        genome_path = "%s/annotation/files/genomes/%s" % (work_dir, genome)
        if not glob.glob(genome_path):
            os.system("mkdir %s" % genome_path)
        if not glob.glob("%s/%s_orfs.fasta" % (genome_path, genome)):
            os.system("getorf -find 0 -minsize %s -table %s  -outseq %s/%s_orfs.fasta -sequence %s%s.fasta" % (minsize, table, genome_path, genome, genomes_for_annotation_dir, genome))
            print("    from genome " + genome + ", using getorf...........DONE!")
        else:
            print("    from genome " + genome + " ORFs had been exported previously")

    print("######################## Block 2 executed successfully! ############################\n\n")






# Block 3. Read protein list of Reference1, produce list of lists - Prot_seq_list - crucial for consecutive steps
#### Prot_seq_list: [["Prot_name_1", "MHKIAIK..."], ["Prot_name_1", "LMKIKIAVV..."]...] - for now at least.
# Requires: genome_list

if "B3" in Blocks_to_run:
    print("\n\n######################## Executing Block 3 #########################################")
    print("Reading list of protein references, making hmms...")
    
    for filename in glob.glob("%s/*_prot.fasta" % protein_ref_dir):
        gene_name = filename.strip().split("/")[-1].split("_")[0]
        
        PROT_FASTA = ImportFasta(filename)
        ref_pos = 0
        for i in range(len(PROT_FASTA)):
            if Size_reference in PROT_FASTA[i][0]:
                ref_pos = i
                break
                
        new_sequence = ""
        seq = PROT_FASTA[ref_pos]
        for aa in seq[1]:
            if aa != "-":
                new_sequence += aa
        seq[1] = new_sequence
        seq[0] = gene_name
        Prot_seq_list.append(seq)
                
        if not glob.glob("%s/%s_prot.hmm" % (protein_ref_dir, gene_name)):
            os.system("hmmbuild --singlemx --amino -o temp %s/%s_prot.hmm %s%s_prot.fasta" % (protein_ref_dir, gene_name, protein_ref_dir, gene_name))
                
    print("DONE!")
    print("Listing %s protein-coding genes that will be searched for in genomes to be annotated:" % len(Prot_seq_list))
    for Protein in Prot_seq_list:
        print(Protein[0], end=', ')

    Gene_range = range(len(Prot_seq_list))

    print("\n######################## Block 3 executed successfully! ############################\n\n")





# Block 4. Export all proteins from Prot_seq_list, search reference2, create HMM
#### Export all proteins from Prot_seq_list as single-fasta files, saving each to a dedicated folder in a new folder "genes"
#### Like this: entry "rpoB" ---> ./genes/rpoB/rpoB.fasta
#### Then, use jackhmmer to find that protein in reference2, and create HMM profile
#### Next, use hmmsearch to search for that HMM in list of ORFs for each genome from a list.
####

if "B4" in Blocks_to_run:
    print("\n\n######################## Executing Block 4 #####################################")
    print("****** Searching for individual proteins from 'Reference annotation' in genomes to be annotated ******")
    
    if not glob.glob("%s/annotation/" % work_dir):
        os.system("mkdir %s/annotation/" % work_dir)
    if not glob.glob("%s/annotation/files" % work_dir):
        os.system("mkdir %s/annotation/files" % work_dir)
    if not glob.glob("%s/annotation/files/genes" % work_dir):
        os.system("mkdir %s/annotation/files/genes" % work_dir)    
    print("Results will be listed in Block 5")
    
    for prot_no in Gene_range:
        Protein = Prot_seq_list[prot_no]
        print("\nProtein: " + Protein[0])
        gene_path = "%s/annotation/files/genes/%s/" % (work_dir, Protein[0])   
        if not glob.glob(gene_path):
            os.system("mkdir %s" % gene_path)
            

        for genome in genome_list:
            if not glob.glob("%s%s_%s.report" % (gene_path, Protein[0], genome)):
                os.system("hmmsearch -E %s --domtblout %s%s_%s.report %s%s_prot.hmm %s/annotation/files/genomes/%s/%s_orfs.fasta > %s%s.hmmsearch.report" % (evalue_for_hmmsearch, gene_path, Protein[0], genome, protein_ref_dir, Protein[0], work_dir, genome, genome, gene_path, genome))
                print("Searching genome %s for %s" % (genome, Protein[0]))
            else:
                print("Genome %s has been searched for protein %s previously" % (genome, Protein[0]))

            
    print("######################## Block 4 executed successfully! ############################\n\n")




# Block 5. For each gene, for each genome, open hmm, parse it, and add entries to Prot_seq_list
#### For each gene, for each genome, opens hmmsearch report and extracts specific values for each hit:
#### [Target_seq_ID, Evalue, Env_start, Env_end, ORF_start_in_genome, ORF_end_in_genome, Orientation]
#### Then, saves them as rows in Hit_table 
# Requires: genome_list created in "Block1"
# Requires: Prot_seq_list created in "Block3"


if "B5" in Blocks_to_run:
    print("\n\n######################## Executing Block 5 #####################################")
    print("Parsing HMMER results for each genome, and each gene; extracting hits; classifying them; updating borders and re-classifying")


    for i in range(0,len(genome_list)):                  ### i - genome coords in genome_list; in Prot_seq_list, this needs to be +2'ed
        genome = genome_list[i]
        print("\n\n### Genome %s" % genome)
        print("Importing genome %s and computing the complement........." % genome, end ="")
        GENOME = ImportFasta("%s%s.fasta" % (genomes_for_annotation_dir, genome_list[i]))
        GENOME_compl = ''
        for nucl in GENOME[0][1]:
            GENOME_compl += Compl[nucl]
        GENOME.append(["", GENOME_compl])
        print("DONE! Genome length: %d bp, complement length: %d bp" % (len(GENOME[0][1]), len(GENOME[1][1])))
        print("Importing ORF list for genome %s........." % genome, end = "")
        ORF_list = ImportFasta("%s/annotation/files/genomes/%s/%s_orfs.fasta" % (work_dir, genome, genome))

        print("DONE!\nNow scrolling through genes.............\n")

        
        ORFs_used = {} ### This distionary contains assignments of all ORFs - and will hopefully help avoid naming the same gene twice!
        
        
        for k in Gene_range:
            Protein = Prot_seq_list[k]
            gene_path = "%s/annotation/files/genes/%s/" % (work_dir, Protein[0])    

            print("Protein: %s" % Protein[0])

            HITS = open(gene_path + Protein[0] + "_" + genome + ".report", 'r')
            Hit_list = []
            for line in HITS:   
                if not line.startswith('#'):
                    Line =  line.split()
                    Hit_list.append(Line[0])
                    if Line[0] not in ORFs_used:
                        ORFs_used[Line[0]] = [[Line[3].split("_")[0], Line[6]]]
                    else:
                        ORFs_used[Line[0]].append([Line[3].split("_")[0], Line[6]])
            HITS.close()
                    
            Hit_table = []
            
            # It's all super easy if Hit_table is empty
            if len(Hit_list) == 0:
                appendix = ["Absent", 0, 0, 0, 0, 0, 0, 0]
                print("        --> No apparent hits. " + Protein[0] + " appears to be absent in " + genome + "\n")
        
        
            else:            
                # That is going to be a pain, but I am going to align all hits....
            
                ali_name = "%s/annotation/files/genes/%s/ref_vs_%s.align" % (work_dir, Protein[0], genome)
                if not glob.glob(ali_name):
                    sequence_DUMP = "%s/annotation/files/genes/%s/%s_seqs_for_alignment.fasta" % (work_dir, Protein[0], genome)
                    TEMP = open(sequence_DUMP, 'w')
                    print(">%s\n%s" % (Protein[0], Protein[1]), file = TEMP)
                    for ORF in ORF_list:
                        if ORF[0].split()[0] in Hit_list:
                            print(">%s\n%s" % (ORF[0], ORF[1]), file = TEMP)
                    TEMP.close()
                    os.system("mafft --op 10 --maxiterate 1000 --thread 2 --quiet --localpair %s > %s" % (sequence_DUMP, ali_name))
                    os.system("rm %s" % sequence_DUMP)
                
                
                ### Cool. Now, importing the alignment and analyzing alignment quality for each sequence
                print(ali_name)
                ALI = ImportFasta(ali_name)
                
                for seq in ALI:
                    aa_match_len = 0          ### no aas aligning between the query gene and the reference
                    aa_start = 0              ### aa position in the query gene that corresponds to the start of the reference
                    aa_end =  0               ### aa position in the query gene that corresponds to the end of the reference
                    ORF_len = 0               ### no aas in the ORF 
                    for aa_no in range(len(seq[1])):
                        if ALI[0][1][aa_no] != "-" and seq[1][aa_no] != "-":
                            aa_match_len += 1
                            aa_end = ORF_len +1
                        if aa_match_len == 1:
                            aa_start = ORF_len
                        if seq[1][aa_no] != "-":
                            ORF_len += 1
                            ### Essentially, for each pos in alignment other than gap: 1) extend ORF_len...
                            
                    aa_len_from_start = ORF_len - aa_start     ### length of the query gene from aa_start all the way to the stop codon
                        
                    seq.append(aa_start)
                    seq.append(aa_end)
                    seq.append(aa_match_len)
                    seq.append(aa_len_from_start)
                    
                    ## seq = ["TETULN_rpoB", "MKIRP...", 2000, 3000, 200, 220]
                        
 
                for ORF_name in Hit_list:
                    for seq in ALI:
                        if seq[0].split()[0] == ORF_name:
                            start_in_genome = int(seq[0].split()[1].strip("["))
                            end_in_genome = int(seq[0].split()[3].strip("]"))
                            orientation = 'F'
                            if start_in_genome > end_in_genome:
                                orientation = 'R'
                            
                            Hit_table.append([ORF_name, seq[-2], seq[-4], seq[-3], start_in_genome, end_in_genome, orientation, seq[-1]])


                print("   Protein: %s, length %d" % (Protein[0], ALI[1][-1]))
                print("    [Target_seq_ID, Match_len, start_in_ORF, end_in_ORF, start_in_genome, end_in_genome, Orientation, len_from_start]:")
                #print("Hit_list: %s, \nHit_table: %s" % (Hit_list, Hit_table))
                for row in Hit_table:
                    print("    ", end = "")
                    for i in range(len(row)):
                        print(row[i], end = '\t')
                    print('')
                        
                reference_len = ALI[0][-2]
                appendix = []

                if len(Hit_table) > 0:   ### Look at first (best) hit. If long enough to be classified as PUTATIVE or FUNCTIONAL, don't look at others.
                    if  Hit_table[0][6] == 'F':
                        Hit_start = int(Hit_table[0][4])+(int(Hit_table[0][2]))*3
                    elif Hit_table[0][6] == 'R':
                        Hit_start = int(Hit_table[0][4])-(int(Hit_table[0][2]))*3
                    Hit_end = int(Hit_table[0][5])
                    Hit_len = Hit_table[0][1]       ### Match length
                    ORF_len_from_start = Hit_table[0][7]    ### 
                    if Hit_len/reference_len > Threshold_functional and ORF_len_from_start/reference_len < Upper_threshold_functional:
                        appendix = ['Functional', 1, int(Hit_table[0][4]), int(Hit_table[0][5]), Hit_start, Hit_end, Hit_table[0][6]]
                        print("        --> Single HMM hit, match_len %.2f percent of query; ORF_len_from_start %.2f percent of query. Classified as a FUNCTIONAL protein.\n" % (Hit_len/reference_len*100, ORF_len_from_start/reference_len*100))

                    elif Hit_len/reference_len > Threshold_putative or ORF_len_from_start/reference_len > Upper_threshold_functional:
                        appendix = ['Putative', 1, int(Hit_table[0][4]), int(Hit_table[0][5]), Hit_start, Hit_end, Hit_table[0][6]]
                        print("        --> Single HMM hit, match_len %.2f percent of query; ORF_len_from_start %.2f percent of query. Classified as a PUTATIVE protein.\n" % (Hit_len/reference_len*100, ORF_len_from_start/reference_len*100))
                
                    else:
                        # OK, the length of the first ORF was below the 'Putative' threshold.
                        # So, we are dealing with an apparent pseudogene. Let's get the hit borders right.
                        # Hit_start should be OK, but I need to update Hit_end position using env values.
                        if len(Hit_table) == 1:
                            appendix = ['Pseudogene', 1, int(Hit_table[0][4]), int(Hit_table[0][5]), Hit_start, Hit_end, Hit_table[0][6]]
                            print("        --> Single HMM hit, match_len %.2f percent of query; ORF_len_from_start %.2f percent of query. Classified as a PSEUDOGENE.\n" % (Hit_len/reference_len*100, ORF_len_from_start/reference_len*100))
                
                        elif len(Hit_table) > 1:
                            print("        --> For the top HMM hit, match_len is %.2f percent of query; ORF_len_from_start %.2f percent of query. Classified as a PSEUDOGENE - but still working on it" % (Hit_len/reference_len*100, ORF_len_from_start/reference_len*100))
                            if  Hit_table[0][6] == 'F':
                                Hit_end = int(Hit_table[0][4]) + int(Hit_table[0][3])*3 -1
                                for i in range(1, len(Hit_table)):
                                    Addl_hit_start = int(Hit_table[i][4])+(int(Hit_table[i][2]))*3
                                    Addl_hit_end = int(Hit_table[i][4]) + int(Hit_table[i][3])*3 - 1
                                    if abs(Addl_hit_end - Hit_start)+1 < 3.6*reference_len and abs(Hit_end - Addl_hit_start) < 3.6*reference_len:
                                        Hit_start = min(Hit_start, Addl_hit_start)
                                        Hit_end = max(Hit_end, Addl_hit_end)

                            elif  Hit_table[0][6] == 'R':
                                Hit_end = int(Hit_table[0][4]) - int(Hit_table[0][3])*3 + 1
                                for i in range(1, len(Hit_table)):
                                    Addl_hit_start = int(Hit_table[i][4])-(int(Hit_table[i][2]))*3
                                    Addl_hit_end = int(Hit_table[i][4]) - int(Hit_table[i][3])*3 + 1
                                    if abs(Addl_hit_end - Hit_start)+1 < 3.6*reference_len and abs(Hit_end - Addl_hit_start) < 3.6*reference_len:
                                        Hit_start = max(Hit_start, Addl_hit_start)
                                        Hit_end = min(Hit_end, Addl_hit_end)
                    
                            appendix = ['Pseudogene', len(Hit_table), int(Hit_table[0][4]), int(Hit_table[0][5]), Hit_start, Hit_end, Hit_table[0][6]]
                            print("        --> PSEUDOGENE borders updated: start at %d; end at %d; length - %.2f percent of query\n" % (Hit_start, Hit_end, (abs(Hit_end-Hit_start)+1)/3/reference_len*100))
            
            
            
            
            
            
            #############
            #############
            #############
            ### In this section, I am updating gene borders
            
            if appendix[0] in ['Functional', 'Putative']:
                change_aa_match_len = 0
                
                if appendix[6] == 'F':
                    #First, check if the three nucleotides immediately after ORF correspond to STOP codon - and if so, update gene end
                    #print(GENOME[0][1][appendix[3]:appendix[3]+3])
                    if GENOME[0][1][appendix[3]:appendix[3]+3] in stop_codon_list:
                        appendix[5] = appendix[3]+3
                    else:
                        print("!!!!!!!!!!!!!!!!!!!!")

                
                    Annotation_start = appendix[4]
                    ORF_start = appendix[2]
                    Updated_start = 0
                    if GENOME[0][1][Annotation_start-1:Annotation_start+2] in start_codon_list:
                        print("   First currently annotated codon is a start codon; annotation not changed")
                    else:
                        print("   First currently annotated codon is NOT a start codon; searching for an alternative start...")
                        Test_start = Annotation_start
                        while Test_start >= ORF_start:
                            if GENOME[0][1][Test_start-1:Test_start+2] in start_codon_list:
                                print("   Start codon found at position", Test_start, "which is", int((Annotation_start-Test_start)/3), "codons BEFORE the originally found one. The annotation will be updated to start at this codon.")
                                Updated_start = Test_start
                                break
                            Test_start -= 3
                            
                        if not Updated_start:     #### If there was no alternative start codon before the original start position...
                            Test_start = Annotation_start
                            while Test_start <= appendix[3]:
                                if GENOME[0][1][Test_start-1:Test_start+2] in start_codon_list:
                                    print("   Start codon found at position", Test_start, "which is", int((Test_start-Annotation_start)/3), "codons AFTER the originally found one. The annotation will be updated to start at this codon.")
                                    Updated_start = Test_start
                                    break
                                Test_start += 3
                                change_aa_match_len -= 1
                    if Updated_start:
                        appendix[4] = Updated_start

                    #Printing the ORF nucleotide sequence with the old and new starting positions
                    seq_pos = int(ORF_start)
                    while seq_pos <= int(appendix[3])+3:
                        print(GENOME[0][1][seq_pos-1:seq_pos+2], end = "")
                        seq_pos += 3
                        if seq_pos == appendix[4]:
                            print("^", end = '')
                        elif seq_pos == Annotation_start and Annotation_start != appendix[4]:
                            print("*", end = '')
                        else:
                            print(" ", end = "")
                    print("")


                ### Now, if the protein is in the Reverse orientation, things get slightly more complicated.
                if appendix[6] == 'R':
                    
                    #First, check if the three nucleotides immediately after ORF correspond to STOP codon - and if so, update gene end
                    #print(GENOME[1][1][appendix[3]-2:appendix[3]-5:-1])
                    if GENOME[1][1][appendix[3]-2:appendix[3]-5:-1] in stop_codon_list:
                        appendix[5] = appendix[3]-3
                    else:
                        print("!!!!!!!!!!!!!!!!!!!!")
          
                    
                    Annotation_start = appendix[4]
                    ORF_start = appendix[2]
                    Updated_start = 0
                    
                        
                    if GENOME[1][1][Annotation_start-1:Annotation_start-4:-1] in start_codon_list:
                        print("   First currently annotated codon is a start codon; annotation not changed")
                    else:
                        print("   First currently annotated codon is NOT a start codon; searching for an alternative start...")
                        ORF_start = appendix[2]
                        Test_start = Annotation_start
                        while Test_start <= ORF_start:
                            if GENOME[1][1][Test_start-1:Test_start-4:-1] in start_codon_list:
                                print("   Start codon found at position", Test_start, "which is", int((Test_start-Annotation_start)/3), "codons BEFORE the originally found one. The annotation will be updated to start at this codon.")
                                Updated_start = Test_start
                                break
                            Test_start += 3
                        if not Updated_start:
                            Test_start = Annotation_start
                            while Test_start >= appendix[3]:
                                if GENOME[1][1][Test_start-1:Test_start-4:-1] in start_codon_list:
                                    print("   Start codon found at position", Test_start, "which is", int((Annotation_start-Test_start)/3), "codons AFTER the originally found one. The annotation will be updated to start at this codon.")
                                    Updated_start = Test_start
                                    break
                                Test_start -= 3
                                change_aa_match_len -= 1
                    if Updated_start:
                        appendix[4] = Updated_start

                    #Printing the ORF nucleotide sequence with the old and new starting positions
                    seq_pos = ORF_start
                    while seq_pos >= appendix[3]-3:
                        print(GENOME[1][1][seq_pos-1:seq_pos-4:-1], end = "")
                        seq_pos -= 3
                        if seq_pos == appendix[4]:
                            print("^", end = '')
                        elif seq_pos == Annotation_start and Annotation_start != appendix[4]:
                            print("*", end = '')
                        else:
                            print(" ", end = "")
                    print("")





            #############
            #############
            #############
            ### NoW, in this section, I am updating gene classification


            

            if appendix[0] in ["Functional", "Putative"]:
                aa_match_len = Hit_table[0][1] + change_aa_match_len
                aa_len_from_start = Hit_table[0][7]
                
                if aa_match_len/reference_len > Threshold_functional and aa_len_from_start/reference_len < Upper_threshold_functional:
                    Gene_classification = "Functional"
                elif aa_match_len/reference_len > Threshold_putative or aa_len_from_start/reference_len > Upper_threshold_functional:
                    Gene_classification = "Putative"
                else:
                    Gene_classification = "Pseudogene"
                    
                if appendix[0] == Gene_classification:
                    print("   Gene classification not changed - remains as %s\n\n" % Gene_classification.upper())
                else:
                    print("   Gene classification changed from %s to %s\n\n" % (appendix[0].upper(), Gene_classification.upper()))
                    appendix[0] = Gene_classification


            Prot_seq_list[k].append(appendix)
        
        
        ##### Finally, here I am checking whether no ORFs have been annotated under two names...
        print("\n\n\n....OK, now let's make sure that no ORFs have been annotated twice...")
        for ORF in ORFs_used:
            if len(ORFs_used[ORF]) > 1:  ### if a given gene has been classified more than once
                #print(ORF, " : ", ORFs_used[ORF])
                ORF_ambigs = {}
                for item in ORFs_used[ORF]:
                    if item[0] not in ORF_ambigs:
                        ORF_ambigs[item[0]] = float(item[1])
                    else:
                        ORF_ambigs[item[0]] = min(ORF_ambigs[item[0]], float(item[1]))
                        
                if len(ORF_ambigs) > 1:    ### so only if hits to a given ORF have been to different genes...
                    real_call = list(ORF_ambigs.keys())[0]
                    for item in ORF_ambigs:
                        if ORF_ambigs[real_call] > ORF_ambigs[item]:
                            real_call = item   ### that way, I identify the gene that ORF really belongs to........
                    print("ORF %s resembles different genes: %s" % (ORF, ORF_ambigs))
                    print("   annotation as gene %s appears to be correct" % real_call)
              
                    for item in ORF_list:
                        if ORF == item[0].strip().split()[0].strip(">"):
                            ORF_end = int(item[0].strip().split()[3].strip("]"))
                    
                    for gene in Prot_seq_list:
                        #if gene[0] in ORF_ambigs and gene[0] != real_call and gene[3] == ORF_end:
                        if gene[0] in ORF_ambigs and gene[-1][3] == ORF_end and gene[0] != real_call:
                            print("      ", gene[0], gene[-1], ".......INCORRECTLY ANNOTATED - ANNOTATION REMOVED!")
                            gene[-1] = ["Absent", 0, 0, 0, 0, 0, 0, 0]
                            #print(gene)

                        #elif gene[0] == real_call:
                        #    print("                     ", gene[0], gene[-1], ".......OK:)"
                    

    print("######################## Block 5 executed successfully! #########################\n\n")















# Block 6. Print annotation summary for each gene
#### Requires: Prot_seq_list from "Block5"

if "B6" in Blocks_to_run:
    print("\n\n######################## Executing Block 6 #####################################")

    print("****** Annotation summary ******")
    print("[Genome, Gene_classification, Number_of_hmmer_hits, Gene_start_in_genome, Gene_end_in_genome, Gene_length_bp]")
    for i in Gene_range:
        print(Prot_seq_list[i][0] + ": ", len(Prot_seq_list[i][1])*3, "bp")
        for k in range(2,len(Prot_seq_list[i])):
            print("    Genome " + genome_list[k-2] + ":", Prot_seq_list[i][k][0], Prot_seq_list[i][k][1], Prot_seq_list[i][k][4], Prot_seq_list[i][k][5], sep = '\t', end = '\t')
            if Prot_seq_list[i][k][0] == 'Absent':
                print('0')
            else:
                print(abs(Prot_seq_list[i][k][5]-Prot_seq_list[i][k][4])+1)

    print("######################## Block 6 executed successfully! #########################\n\n")







# Block 7. Saves the Prot_seq_list table to file -- 'Protein_table.sav'
#### Requires: Prot_seq_list from "Block5"

if "B7" in Blocks_to_run:
    print("\n\n######################## Executing Block 7 #####################################")

    #### Save Prot_seq_list to file:
    print("Saving the result table to file -- 'Protein_table.sav' -- to be reopened later ...")

    PROTEIN_LIST = open("%s/annotation/Protein_table.sav" % work_dir, 'w')
    for row in Prot_seq_list:
        for i in range(0, 2):
            print(row[i], end = '\t', file = PROTEIN_LIST)
        for i in range(2, len(row)):
            for k in range(len(row[i])):
                print(row[i][k], end = ' ', file = PROTEIN_LIST)
            print('\t', end = '', file = PROTEIN_LIST)
        print('\n', end = '', file = PROTEIN_LIST)
    
    print("DONE!")
    print("######################## Block 7 executed successfully! ########################\n\n")



# Block 8. Reopens 'Protein_table', reading the output as Prot_seq_list
#### Requires: Protein_table.sav

if "B8" in Blocks_to_run:
    print("\n\n######################## Executing Block 8 #####################################")
    print("Reopening the result table which had been saved to file -- 'Protein_table' ...")

    PROTEIN_LIST = open("%s/annotation/Protein_table.sav" % work_dir, 'r')
    Prot_seq_list = []
    for line in PROTEIN_LIST:
        Line = line.strip().split('\t')
        for i in range (2, len(Line)):
            Line[i] = Line[i].strip().split( )
            for k in range(1, len(Line[i])-1):
                Line[i][k] = int(Line[i][k])
        Prot_seq_list.append(Line)
   
#print(Prot_seq_list)

    print("DONE!")
    print("######################## Block 8 executed successfully! ########################\n\n")








#### Block 11. Searches genomes for rRNA genes, generates list "rRNA_gene_list"
# Requires: genome_list

if "B11" in Blocks_to_run:
    print("\n\n######################## Executing Block 11 ####################################")
    print("Scans genomes for rRNA genes using rnammer 1.2, adding them to newly created rRNA_gene_list")
    print("Then, uses nhmmer and a list of curated references to search the genomes for ncRNA and rRNA")
    print("All finds are added to the same list\n")

    if not glob.glob("%s/annotation/" % work_dir):
        os.system("mkdir %s/annotation/" % work_dir)
    if not glob.glob("%s/annotation/files" % work_dir):
        os.system("mkdir %s/annotation/files" % work_dir)
    if not glob.glob("%s/annotation/files/rRNA" % work_dir):
        os.system("mkdir %s/annotation/files/rRNA" % work_dir)

    rRNA_gene_list = [["16s_rRNA", ""], ["23s_rRNA", ""], ["5s_rRNA", ""], ["ssrA", ""], ["rnpB", ""]]
    #ncRNAs = ["ssrA", "rnpB"]
    ### Entries will look like this: [[gene_name, '', ['Functional/Absent', 'Gene_start', 'Gene_end', 'Orientation']...
    ### These are the genes that will be looked for
    
    for filename in glob.glob("%s*.fasta" % rRNA_ref_dir):
        gene_name = filename.strip().split("/")[-1].split(".")[0]
        if not glob.glob("%s%s.hmm" % (rRNA_ref_dir, gene_name)):
            os.system("hmmbuild --singlemx --dna -o temp %s%s.hmm %s%s.fasta" % (rRNA_ref_dir, gene_name, rRNA_ref_dir, gene_name))


    for genome in genome_list:
    
        rRNAs_in_genome = []
    
    
        """ This portion of the script is about searching for rRNAs using rnammer... an old tool hard or impossible to get to work. Keeping here just in case!
        print("Running rnammer on genome %s........." % genome, end = '')
        
        path_to_genome = genomes_for_annotation_dir + genome + ".fasta"
        path_to_output = "%s/annotation/files/rRNA/%s_rnammer.gff" % (work_dir, genome)
        if not glob.glob(path_to_output):
            os.system("rnammer -S bac -m lsu,ssu,tsu -gff - < %s > %s" % (path_to_genome, path_to_output))
    
        print("DONE!\nNow parsing the rnammer output...............", end = '')
        
        
        RNAMMER_OUTPUT = open(path_to_output, "r")
        for line in RNAMMER_OUTPUT:
            if not line.startswith("#"):
                rRNA_gene = line.strip().split()
                if abs(int(rRNA_gene[3])-int(rRNA_gene[4])) < max_rRNA_gene_length[rRNA_gene[8]] and abs(int(rRNA_gene[3])-int(rRNA_gene[4])) > min_rRNA_gene_length[rRNA_gene[8]]:
                    rRNAs_in_genome.append([rRNA_gene[8], rRNA_gene[3], rRNA_gene[4], rRNA_gene[6]]) ### [Gene_name, Gene_start, Gene_end, strand]
        RNAMMER_OUTPUT.close()
        """
        
        print("Searching for rRNA and ncRNA genes using references...................", end = '')

        
        for gene_name in rRNA_gene_list:
            gene = gene_name[0]
        
            appendix = ['Absent', '', '', '']
            output_path = "%s/annotation/files/rRNA/%s" % (work_dir, gene)
            if not glob.glob(output_path):
                os.system("mkdir %s" % output_path)   
            
            if not glob.glob("%s/%s.fasta" % (output_path, genome)):
                os.system("nhmmer -E 1e-10 -A %s/%s.sto %s%s.hmm %s%s.fasta > %s/%s.report" % (output_path, genome, rRNA_ref_dir, gene, genomes_for_annotation_dir, genome, output_path, genome))
            if glob.glob("%s/%s.sto" % (output_path, genome)):
                os.system("esl-reformat fasta %s/%s.sto > %s/%s.fasta" % (output_path, genome, output_path, genome))
                FASTA = ImportFasta("%s/%s.fasta" % (output_path, genome))
                ### seq name goes like this: >OKARIM1/128763-128540 [subseq from] OKARIM1
                if FASTA[0][0] != "":
                    gene_start = int(FASTA[0][0].strip().split()[0].split("/")[1].split("-")[0])
                    gene_end = int(FASTA[0][0].strip().split()[0].split("/")[1].split("-")[1])  
                    start_in_genome = min(gene_start, gene_end)
                    end_in_genome = max(gene_start, gene_end)
                    orientation = "+"
                    if gene_start > gene_end:
                        orientation = "-"
                    appendix = [gene, start_in_genome, end_in_genome, orientation]
                        
            if not appendix[0] == 'Absent' and abs(end_in_genome-start_in_genome) > min_rRNA_gene_length[appendix[0]]:
                rRNAs_in_genome.append(appendix)
        print(rRNAs_in_genome)
        
        for gene in rRNA_gene_list:
            found = 0
            for entry in rRNAs_in_genome:
                if gene[0] == entry[0]:
                    orientation = 'F'
                    if entry[3] == '-':
                        orientation = 'R'
                    gene.append(['Functional', entry[1], entry[2], orientation])
                    found = 1
                    break
            if found == 0 and gene[0]:
                gene.append(['Absent', '', '', ''])
                
        print("DONE!\n   rRNA and ncRNA genes found: ", end = '')
        for entry in rRNA_gene_list:
            if entry[-1][0] == 'Functional':
                print(entry[0], end = ", ")
        print("\n")
    
        #print("   rRNAs_in_genome: \n", rRNAs_in_genome)    
        
            
    print("The final rRNA_gene_list:\n", rRNA_gene_list, "\n")
            

    print("######################## Block 11 executed successfully! ########################\n\n")





















#### Block 12. Searches genomes for tRNA genes, generates list "tRNA_gene_list"
# Requires: genome_list

if "B12" in Blocks_to_run:
    print("\n\n######################## Executing Block 12 ####################################")
    print("Scanning genomes for tRNA genes using tRNAscan-SE v.2.0.5, adds them to newly created tRNA_gene_list\n")

    if not glob.glob("%s/annotation/" % work_dir):
        os.system("mkdir %s/annotation/" % work_dir)
    if not glob.glob("%s/annotation/files" % work_dir):
        os.system("mkdir %s/annotation/files" % work_dir)
    if glob.glob("%s/annotation/files/tRNA" % work_dir):
        os.system("rm -r %s/annotation/files/tRNA/*" % work_dir)
    if not glob.glob("%s/annotation/files/tRNA" % work_dir):
        os.system("mkdir %s/annotation/files/tRNA" % work_dir)

    tRNA_gene_list = [['AAA', ''], ['AAC', ''], ['AAG', ''], ['AAT', ''], ['ACA', ''], ['ACC', ''], ['ACG', ''], ['ACT', ''], ['AGA', ''], ['AGC', ''], ['AGG', ''], ['AGT', ''], ['ATA', ''], ['ATC', ''], ['ATG', ''], ['ATT', ''], ['CAA', ''], ['CAC', ''], ['CAG', ''], ['CAT', ''], ['CCA', ''], ['CCC', ''], ['CCG', ''], ['CCT', ''], ['CGA', ''], ['CGC', ''], ['CGG', ''], ['CGT', ''], ['CTA', ''], ['CTC', ''], ['CTG', ''], ['CTT', ''], ['GAA', ''], ['GAC', ''], ['GAG', ''], ['GAT', ''], ['GCA', ''], ['GCC', ''], ['GCG', ''], ['GCT', ''], ['GGA', ''], ['GGC', ''], ['GGG', ''], ['GGT', ''], ['GTA', ''], ['GTC', ''], ['GTG', ''], ['GTT', ''], ['TAA', ''], ['TAC', ''], ['TAG', ''], ['TAT', ''], ['TCA', ''], ['TCC', ''], ['TCG', ''], ['TCT', ''], ['TGA', ''], ['TGC', ''], ['TGG', ''], ['TGT', ''], ['TTA', ''], ['TTC', ''], ['TTG', ''], ['TTT', '']]


        
    for genome in genome_list:
        print("Running tRNAscan-SE on genome %s........." % genome, end = '')
        
        path_to_genome = genomes_for_annotation_dir + genome + ".fasta"
        path_to_output = "%s/annotation/files/tRNA/%s_tRNAscan.gff" % (work_dir, genome)
        
        if not glob.glob(path_to_output):
            os.system("tRNAscan-SE -D -B -q -o %s %s > temp" % (path_to_output, path_to_genome))
    
        print("DONE!\nNow parsing the output.....................", end = '')
        
        tRNAs_in_genome = []
        tRNAscan_OUTPUT = open(path_to_output, "r")
        tRNAscan_OUTPUT.readline() ### Skip the first three lines of the output!
        tRNAscan_OUTPUT.readline() ### Skip the first three lines of the output!
        tRNAscan_OUTPUT.readline() ### Skip the first three lines of the output!

        for line in tRNAscan_OUTPUT:
            tRNA_gene = line.strip().split()
            tRNAs_in_genome.append([tRNA_gene[5], tRNA_gene[4], tRNA_gene[2], tRNA_gene[3]]) ### [Anticodon, AA/Pseudo, Gene_start, Gene_end]
        tRNAscan_OUTPUT.close()
        

        for gene in tRNA_gene_list:
            new_gene = []
            found = 0
            for entry in tRNAs_in_genome:
                if gene[0] == entry[0]:
                    if entry[2] < entry[3]:
                        orientation = 'F'
                    else:
                        orientation = 'R'
                    new_gene.append([entry[1], entry[2], entry[3], orientation])
                    found = 1
            if found == 0:
                gene.append([['Absent', '', '', '']])
            else:
                gene.append(new_gene)

                
        print("DONE!\nGenes found:", end = '')
        for entry in tRNA_gene_list:
            if entry[-1][0][0] != 'Absent':
                print(entry[0], end = "")
                if len(entry[-1]) > 1:
                    print("(%d), " % len(entry[-1]), end = "")
                else:
                    print(", ", end = "")
        print("\n\n")
        
    ### Printing summary
    print("The number of tRNA genes with a given anticodon annotated in each of the genomes:")
    print(" ", end = '')
    for genome in genome_list:
        print(genome, end = ' ')
    print('')
    for line in tRNA_gene_list:
        present = 0
        for i in range(2, len(line)):
            if line[i][0][0] != 'Absent':
                present = 1
        if present == 1:
            print(line[0], end = ' ')
            for i in range(2, len(line)):
                if line[i][0][0] == 'Absent':
                    print(' ', end = ' ')
                else:
                    print(len(line[i]), end = ' ')
            print('')
    print("\n\n")    
    
    
    
    #print("The final tRNA_gene_list:\n", tRNA_gene_list, "\n")
            
    print("######################## Block 12 executed successfully! ########################\n\n")












#### Block 15. Outputs "prot_seq_list", "rRNA_gene_list" and "tRNA_gene_list" as table that can be used for genome contents comparison using Processing
# Requires: genome_list, prot_seq_list, rRNA_gene_list, tRNA_gene_list

if "B15" in Blocks_to_run:
    print("\n\n######################## Executing Block 15 ####################################")
    print("Outputs 'prot_seq_list', 'rRNA_gene_list' and 'tRNA_gene_list' as table that can be used for illustration using Processing.\nThe table will be saved as will be saved as '%s/annotation/gene_contents_table.csv'\n" % work_dir)

    # That's how the table will look like:
    #gene_name	gene_type	funct_category	Genome0001	Genome02	Genome0003	Genome04 ...
    #geneAAA	Protein	Methionine_biosynth	Functional	Putative	Pseudogene	Absent
    #geneAAB	tRNA	tRNA	Absent	Functional	Putative	Absent
    #...

    FEATURE_TABLE = open(gene_detail_table, 'r', encoding = "ISO-8859-1")
    REFERENCE_FEATURES = []
    
    for line in FEATURE_TABLE:
        LINE = line.strip('\n').split('\t')
        REFERENCE_FEATURES.append(LINE)

    FEATURE_TABLE.close()

    
    if not glob.glob("%s/annotation/" % work_dir):
        os.system("mkdir %s/annotation/" % work_dir)
    
    GENE_TABLE = open("%s/annotation/genome_contents_table.csv" % work_dir, 'w')
    TABLE = []

    First_row = ["gene_name", "gene_type", "funct_category"]
    for genome in genome_list:
        First_row.append(genome)
    TABLE.append(First_row)

    if len(Prot_seq_list) > 0:
        for i in Gene_range:
            gene = Prot_seq_list[i][0]
            classification = ""
            for entry in REFERENCE_FEATURES:
                if entry[0] == gene:
                    classification = entry[7]
                    break
            
            New_row = [gene, 'Protein', classification]
            for k in range(2, len(Prot_seq_list[i])):
                New_row.append(Prot_seq_list[i][k][0])
            TABLE.append(New_row)
            
    if len(rRNA_gene_list) > 0:
        for i in range(len(rRNA_gene_list)):
            if rRNA_gene_list[i][0] == "rnpB":
                New_row = [rRNA_gene_list[i][0], 'ncRNA', 'translation']
            elif rRNA_gene_list[i][0] == "ssrA":
                New_row = [rRNA_gene_list[i][0], 'tmRNA', 'translation']
            else:
                New_row = [rRNA_gene_list[i][0], 'rRNA', 'translation-rRNA']
            for k in range(2, len(rRNA_gene_list[i])):
                New_row.append(rRNA_gene_list[i][k][0])
            TABLE.append(New_row)
            
    if len(tRNA_gene_list) > 0:
        for i in range(len(tRNA_gene_list)):
            anticodon_present_anywhere = 0
            for j in range(2, len(tRNA_gene_list[i])):
                if tRNA_gene_list[i][j][0][0] != 'Absent':
                    anticodon_present_anywhere = 1            
            if anticodon_present_anywhere == 1:
                New_row = [tRNA_gene_list[i][0], 'tRNA', 'translation-tRNA']
                for j in range(2, len(tRNA_gene_list[i])):
                    if len(tRNA_gene_list[i][j]) > 1:
                        New_row.append(len(tRNA_gene_list[i][j]))
                    elif tRNA_gene_list[i][j][0][0] == 'Pseudo':
                         New_row.append('Pseudogene')
                    elif tRNA_gene_list[i][j][0][0] == 'Absent':
                         New_row.append('Absent')
                    else:
                        New_row.append('Functional')
                TABLE.append(New_row)
            
            
    
    for row in TABLE:
        for entry in row:
            print(entry, end = ',', file = GENE_TABLE)
        print('\n', end = '', file = GENE_TABLE)
    
    GENE_TABLE.close()
    print("DONE!")



    print("######################## Block 15 executed successfully! ########################\n\n")

















#### Block 16. Outputs all functional proteins and RNAs from across the genomes, aligns them
# Requires: genome_list
# Requires: prot_seq_list, rRNA_gene_list, tRNA_gene_list

if "B16" in Blocks_to_run:
    print("\n\n######################## Executing Block 16 ####################################")
    print("Outputs all functional protein and RNA sequences from across the genomes, aligns them, exports alignments to work_dir /annotation/alignment/ \n")

    # Make alignment folder
    # For each gene, export functional copies from each genome to ./annotation/files/genes/gene_name
    # Keep a list of copies of each gene
    # Align genes using mafft
    
    Compl = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N':'N'}
    if not glob.glob("%s/annotation/alignments" % work_dir):
        os.system("mkdir %s/annotation/alignments/" % work_dir)
    if not glob.glob("%s/annotation/files/genomes" % work_dir):
        os.system("mkdir %s/annotation/files/genomes" % work_dir)
    if not glob.glob("%s/annotation/files/genes" % work_dir):
        os.system("mkdir %s/annotation/files/genes" % work_dir)
    
    ### First, remove any existing genome concatenation, then create a new one:
    ### Concatenate all genome fasta files, then import concatenation
    
    if glob.glob("%s/annotation/files/genomes/concatenated_genomes.fasta" % work_dir): 
        os.system("rm %s/annotation/files/genomes/concatenated_genomes.fasta" % work_dir)
        
    for genome in genome_list:
        print(genome, "%s%s.fasta" % (genomes_for_annotation_dir, genome))   ####### TEST!!!
        os.system("cat %s%s.fasta >> %s/annotation/files/genomes/concatenated_genomes.fasta" % (genomes_for_annotation_dir, genome, work_dir))

        
    GENOMES = ImportFasta("%s/annotation/files/genomes/concatenated_genomes.fasta" % (work_dir))
    #GENOMES is a table where [i][0] is a name, [i][1] genomic sequence - and I'm about to add [i][2] - complement


    for i in range(len(GENOMES)):                  ### i - genome coords in genome_list; in Prot_seq_list, this needs to be +2'ed
        GENOME_compl = ''
        for nucl in GENOMES[i][1]:
            GENOME_compl += Compl[nucl]
        GENOMES[i].append(GENOME_compl)
        print("Imported and complemented genome %s of %d bp" % (GENOMES[i][0], len(GENOMES[i][1])))

    
    if len(Prot_seq_list) > 0:
        print("\nNow exporting and aligning all proteins from across the genomes....\n")           

        for k in Gene_range:
            
            print("Exporting and aligning protein %s..............." % Prot_seq_list[k][0], end = "")
            
            ### First, export all copies of a given gene identified as 'Functional'
            GENE_EXPORTS = open("%s/annotation/files/genes/%s/%s.fasta" % (work_dir, Prot_seq_list[k][0], Prot_seq_list[k][0]), "w")
            seq_count = 0
            for i in range(len(genome_list)):
                if Prot_seq_list[k][i+2][0] == "Functional":
                    seq_count += 1
                    print(">%s_%s" % (genome_list[i], Prot_seq_list[k][0]), file = GENE_EXPORTS)
                    
                    if Prot_seq_list[k][i+2][6] == "F":
                        print(GENOMES[i][1][int(Prot_seq_list[k][i+2][4])-1:int(Prot_seq_list[k][i+2][5])], file = GENE_EXPORTS)

                    elif Prot_seq_list[k][i+2][6] == "R":
                        seq_pos = int(Prot_seq_list[k][i+2][4])
                        while seq_pos >= int(Prot_seq_list[k][i+2][5]):
                            print(GENOMES[i][2][seq_pos-1], end = "", file = GENE_EXPORTS)
                            seq_pos -= 1
                        print("", file = GENE_EXPORTS)  
            
            GENE_EXPORTS.close()
            
            
            ### Now, align!
            if seq_count > 1:
                os.system("align_nucl_by_codon.py %s/annotation/files/genes/%s/%s.fasta %s/annotation/alignments/%s.fasta %s" % (work_dir, Prot_seq_list[k][0], Prot_seq_list[k][0], work_dir, Prot_seq_list[k][0], Translation_table))
                print("alignment of %s sequences complete!" % seq_count)
            elif seq_count == 1:
                os.system("cp %s/annotation/files/genes/%s/%s.fasta %s/annotation/alignments/%s.fasta" % (work_dir, Prot_seq_list[k][0], Prot_seq_list[k][0], work_dir, Prot_seq_list[k][0]))
                print("only one sequence identified as functional and output to alignments!")
            else:
                print("alignment aborted... no sequences identified as functional!")

    else:
        print("Protein-coding gene list is empty, skipping...")



    if len(rRNA_gene_list) > 0:
        print("len(rRNA_gene_list) ", len(rRNA_gene_list), "\n", rRNA_gene_list)
        print("\nNow exporting and aligning rRNAs from across the genomes....\n")           

        for gene in rRNA_gene_list:
            print("Exporting and aligning rRNA gene %s..............." % gene[0], end = "")

            if not glob.glob("%s/annotation/files/genes/%s" % (work_dir, gene[0])):
                os.system("mkdir %s/annotation/files/genes/%s" % (work_dir, gene[0]))


            ### First, export all copies of the gene identified as 'Functional'
            seq_count = 0
            GENE_EXPORTS = open("%s/annotation/files/genes/%s/%s.fasta" % (work_dir, gene[0], gene[0]), "w")
            for i in range(len(genome_list)):
                if gene[i+2][0] == "Functional":
                    print(">%s_%s" % (genome_list[i], gene[0]), file = GENE_EXPORTS)
                    
                    if gene[i+2][3] == "F":
                        print(GENOMES[i][1][int(gene[i+2][1])-1:int(gene[i+2][2])], file = GENE_EXPORTS)
                        
                    elif gene[i+2][3] == "R":
                        seq_pos = int(gene[i+2][2])
                        while seq_pos >= int(gene[i+2][1]):  
                            print(GENOMES[i][2][seq_pos-1], end = "", file = GENE_EXPORTS)
                            seq_pos -= 1
                        print("", file = GENE_EXPORTS)
                seq_count += 1
            
            GENE_EXPORTS.close()
            
            
            ### Now, align!
            if seq_count > 1:
                os.system("mafft --adjustdirection --maxiterate 1000 --thread 2 --quiet --localpair %s/annotation/files/genes/%s/%s.fasta > %s/annotation/alignments/%s.fasta" % (work_dir, gene[0], gene[0], work_dir, gene[0]))
                print("DONE!")
            elif seq_count == 1:
                os.system("cp %s/annotation/files/genes/%s/%s.fasta %s/annotation/alignments/%s.fasta" % (work_dir, gene[0], gene[0], work_dir, gene[0]))
                print("one sequence exported, no alignments done!")            
            else:
                print("no sequences of gene %s annotated and available for alignment!" % gene[0]) 
    else:
        print("rRNA gene list is empty, skipping...")


            
            
    print("######################## Block 16 executed successfully! ########################\n\n")            

    


















#### Block 18. Combines "gene_lists", sorts the contents
# Requires: genome_list
# Requires: Prot_seq_list, rRNA_gene_list, tRNA_gene_list

if "B18" in Blocks_to_run:
    print("\n\n######################## Executing Block 18 ####################################")
    print("Combines Protein, rRNA and tRNA gene lists, adds them to annotation.\nThe output is a triple-nested list 'all_gene_table', where first level corresponds to genome, second - all present genes, third - features of the entry.")
    
    all_gene_table = [] # list where items represent complete lists of genes for a given genome
    
    for g in range(len(genome_list)):
        all_genes_in_the_genome = [] 
        # this is a list of all genes in this genome. Each gene has the following parameters:
        # [ID_to_be_added; Protein/rRNA/tRNA/other; gene_name; Functional/Putative/Pseudogene; gene_start; gene_end; gene_orientation]
        
        if len(Prot_seq_list) > 0:
           for gene in Prot_seq_list:
               if len(gene) > 2:
                   if gene[g+2][0] != "Absent":
                       all_genes_in_the_genome.append(["", "Protein", gene[0], gene[g+2][0], int(min(gene[g+2][4], gene[g+2][5])), int(max(gene[g+2][4], gene[g+2][5])), gene[g+2][6]])
        
        
        if len(rRNA_gene_list) > 0:
           for gene in rRNA_gene_list:
               if len(gene) > 2:
                   if gene[g+2][0] != "Absent":
                       if gene[0] == "ssrA":
                           classification = "tmRNA"
                       elif gene[0] == "rnpB":
                           classification = "ncRNA"
                       else:
                           classification = "rRNA"
                       all_genes_in_the_genome.append(["", classification, gene[0], "Functional", int(min(gene[g+2][1], gene[g+2][2])), int(max(gene[g+2][1], gene[g+2][2])), gene[g+2][3]])
                       
                       
        if len(tRNA_gene_list) > 0:
            # Somewhat different structure: each line looks like this - ['GTC', '', [['Asp', '129534', '129462', 'R']], [['Pseudo', '56152', '56222', 'F']], [['Absent', '', '', '']]...]
            for gene in tRNA_gene_list:
                if gene[g+2][0][0] != "Absent":
                    for m in range(len(gene[g+2])):
                        all_genes_in_the_genome.append(["", "tRNA", gene[0], gene[g+2][m][0], int(min(gene[g+2][m][1], gene[g+2][m][2])), int(max(gene[g+2][m][1], gene[g+2][m][2])), gene[g+2][m][3]])
        

        
        
        # sorting all entries by start_position!
        all_genes_in_the_genome = sorted(all_genes_in_the_genome, key=itemgetter(4), reverse=False)
        
        
        # labeling features
        for h in range(len(all_genes_in_the_genome)):
            all_genes_in_the_genome[h][0] = "%s_%03d" % (genome_list[g], h+1)
                
        
        print("\n----------------------------\n", "Genome ", genome_list[g], sep = '')
        for gene in all_genes_in_the_genome:
            for item in gene:
                print(item, end = "\t")
            print("")            

        # output "all_genes_in_the_genome" of each genome for further processing3 analysis 
        with open("{}/annotation/all_genes_{}.csv".format(work_dir,genome_list[g]), 'w') as GENE_TABLE:
            for row in all_genes_in_the_genome:
                for entry in row:
                    print(entry, end = ',', file = GENE_TABLE)
                print('\n', end = '', file = GENE_TABLE)        
        
        all_gene_table.append(all_genes_in_the_genome)
    
    
    
    print("######################## Block 18 executed successfully! ########################\n\n")











#### Block 19. For all proteins & rRNA genes, provides details of the annotation
# Requires: genome_list
# Requires: all_gene_table

if "B19" in Blocks_to_run:
    print("\n\n######################## Executing Block 19 ####################################")
    print("Updating all_gene_table with gene details for the primary reference, from file %s\n" % gene_detail_table)
    print("Printing all updated entries...\n")
    
    
        ### Reading feature_table, generated from the reference's GB feature table file
        ### using script process_GB_feature_table.py 
    FEATURE_TABLE = open(gene_detail_table, 'r', encoding = "ISO-8859-1")
    REFERENCE_FEATURES = []
    
    for line in FEATURE_TABLE:
        LINE = line.strip('\n').split('\t')
        REFERENCE_FEATURES.append(LINE)

    FEATURE_TABLE.close()
    
    
#######################    0     1      2          3           4           5           6           7      
#REFERENCE_FEATURES = [["ID", "type", "gene", "locus_tag", "product", "EC_number", "protein_id", "other"]]

#######################      0     1                           2                   3                     4           5             6     
#all_genes_in_the_genome = [[ID; Protein/rRNA/tRNA/other; gene_name; Functional/Putative/Pseudogene; gene_start; gene_end; gene_orientation]
# -> feature

#CP008699.1	Genbank	gene	10960	11862	.	-	.	ID=gene11;Name=HCTETULN_012;gbkey=Gene;gene_biotype=protein_coding;locus_tag=HCTETULN_012
#CP008699.1	Genbank	CDS	10960	11862	.	-	0	ID=cds11;Parent=gene11;gbkey=CDS;product=hypothetical protein;protein_id=AIC63763.1;transl_table=4

#CP008699.1	Genbank	gene	95178	96149	.	+	.	ID=gene126;Name=etfA;gbkey=Gene;gene=etfA;gene_biotype=protein_coding;locus_tag=HCTETULN_127
#CP008699.1	Genbank	CDS	95178	96149	.	+	0	ID=cds112;Parent=gene126;Dbxref=NCBI_GP:AIC63864.1;Name=AIC63864.1;Note=ab initio prediction:Prodigal:2.60%2Csimilar to AA sequence:UniProtKB:P38974;gbkey=CDS;gene=etfA;product=Electron transfer flavoprotein large subunit;protein_id=AIC63864.1;transl_table=4

#CP008699.1	Genbank	gene	46580	46652	.	+	.	ID=gene58;Name=HCTETULN_059;gbkey=Gene;gene_biotype=tRNA;locus_tag=HCTETULN_059
#CP008699.1	Genbank	tRNA	46580	46652	.	+	.	ID=rna4;Parent=gene58;Note=COORDINATES: profile:Aragorn:1.2.34%2CtRNA-Trp(tca);gbkey=tRNA;product=tRNA-Trp

#CP008699.1	Genbank	gene	96677	98158	.	+	.	ID=gene128;Name=HCTETULN_129;gbkey=Gene;gene_biotype=rRNA;locus_tag=HCTETULN_129
#CP008699.1	Genbank	rRNA	96677	98158	.	+	.	ID=rna13;Parent=gene128;Note=COORDINATES: profile:RNAmmer:1.2%2CHCTETULN_128;gbkey=rRNA;product=16S ribosomal RNA

   
    for all_genes_in_the_genome in all_gene_table:
        for feature in all_genes_in_the_genome:             					# For each Present/Putative/Pseudo gene:
            attribute_gene = "ID=gene_%s;locus_tag=%s;" % (feature[0].split("_")[1], feature[0])	#     Get info to be included in the "gene" entry
            attribute2 =  "ID=product_%s;Parent=gene_%s;locus_tag=%s;" % (feature[0].split("_")[1],feature[0].split("_")[1],feature[0])		#     and in the "CDS/rRNA..." entry
        
            feature_no = 0											# For each entry in all_genes..., I need to find the position of the reference!
            for i in range(len(REFERENCE_FEATURES)):
                if REFERENCE_FEATURES[i][0] == feature[2]:
                    feature_no = i
         
            if feature[1] == "Protein" and (feature[3] == "Functional" or feature[3] == "Putative"):
            
                attribute_gene += "gbkey=Gene;gene_biotype=protein_coding;"
                attribute2 += "gbkey=CDS;"
                # attribute_gene: ID,
                # attribute_CDS: ID, locus_tag, gene, product, eC_number, tranl_table=4, codon_start=1, inference (similar to protein_id)
            
                if REFERENCE_FEATURES[feature_no][2] != '' and feature_no > 0:
                    if not REFERENCE_FEATURES[feature_no][2].startswith("xxx"):
                        attribute_gene += "gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                        attribute_gene += "name=%s;" % REFERENCE_FEATURES[feature_no][2]
                        attribute2 += "gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                    else:   
                        attribute_gene += "name=%s;" % feature[0] 
                
                if REFERENCE_FEATURES[feature_no][4] != '' and feature_no > 0:
                    attribute2 += "product=%s;" % REFERENCE_FEATURES[feature_no][4]
            
                if REFERENCE_FEATURES[feature_no][5] != '' and feature_no > 0 and str(REFERENCE_FEATURES[feature_no][5]) != '0':
                    attribute2 += "eC_number=%s;" % REFERENCE_FEATURES[feature_no][5].strip("EC:")
            
                if REFERENCE_FEATURES[feature_no][6] != '' and feature_no > 0 and str(REFERENCE_FEATURES[feature_no][6]) != '0':
                    if REFERENCE_FEATURES[feature_no][6].startswith("WP_"):
                        attribute2 += "note=inference: phmmer, similar to AA sequence: %s;" % REFERENCE_FEATURES[feature_no][6]
                    else:
                        attribute2 += "note=inference: phmmer, similar to AA sequence: UniProtKB:%s;" % REFERENCE_FEATURES[feature_no][6]
                                        
                if feature[3] == "Putative":
                    attribute_gene += "note=putative pseudogene;"
                    attribute2 += "note=putative pseudogene;"


                attribute2 += "transl_table=%d;codon_start=1" % Translation_table
            
                feature.append(attribute_gene)
                feature.append(attribute2)
            
                print(feature)
                


            if feature[1] == "Protein" and feature[3] == "Pseudogene":

                attribute_gene += "gbkey=Gene;gene_biotype=Pseudogene;"            
                # attribute_gene: ID, locus_tag, gene, gene_desc (product), inference (similar to protein_id), pseudo

                if REFERENCE_FEATURES[feature_no][2] != '' and feature_no > 0:
                    if not REFERENCE_FEATURES[feature_no][2].startswith("xxx"):
                        attribute_gene += "gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                        attribute_gene += "name=%s;" % REFERENCE_FEATURES[feature_no][2]
                    else:   
                        attribute_gene += "name=%s;" % feature[0] 
                
                if REFERENCE_FEATURES[feature_no][4] != '' and feature_no > 0:
                    attribute_gene += "gene_desc=%s;" % REFERENCE_FEATURES[feature_no][4]
                        
                if REFERENCE_FEATURES[feature_no][6] != '' and feature_no > 0 and str(REFERENCE_FEATURES[feature_no][6]) != '0':
                    if REFERENCE_FEATURES[feature_no][6].startswith("WP_"):
                        attribute_gene += "note=inference: phmmer, similar to AA sequence: %s;" % REFERENCE_FEATURES[feature_no][6]
                    else:
                        attribute_gene += "note=inference: phmmer, similar to AA sequence: UniProtKB:%s;" % REFERENCE_FEATURES[feature_no][6]

                attribute_gene += "pseudogene=allelic"
            
                feature.append(attribute_gene)
                feature.append("")
            
                print(feature)


#CP008699.1	Genbank	gene	96677	98158	.	+	.	ID=gene128;Name=HCTETULN_129;gbkey=Gene;gene_biotype=rRNA;locus_tag=HCTETULN_129
#CP008699.1	Genbank	rRNA	96677	98158	.	+	.	ID=rna13;Parent=gene128;Note=COORDINATES: profile:RNAmmer:1.2%2CHCTETULN_128;gbkey=rRNA;product=16S ribosomal RNA


            if feature[1] == "rRNA":
            
                attribute_gene += "gbkey=Gene;gene_biotype=rRNA;"            
                attribute2 += "gbkey=rRNA;"
            
                # attribute_gene: ID, locus_tag, gene
                # attribute_rRNA: ID, locus_tag, product
            
                if REFERENCE_FEATURES[feature_no][2] != '':
                    attribute_gene += "gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                    attribute_gene += "name=%s;" % REFERENCE_FEATURES[feature_no][2]
                
                if REFERENCE_FEATURES[feature_no][4] != '':
                    attribute2 += "product=%s;" % REFERENCE_FEATURES[feature_no][4]
                
                attribute2 += "note=inference:RNAmmer 1.2"
                
                feature.append(attribute_gene)
                feature.append(attribute2)
            
                print(feature)


            if feature[1] == "ncRNA":
            
                # attribute_gene: ID, locus_tag, gene
                # attribute_rRNA: ID, locus_tag, product
            
                if REFERENCE_FEATURES[feature_no][2] != '':
                    attribute_gene += "gbkey=Gene;gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                
                if REFERENCE_FEATURES[feature_no][4] != '':
                    attribute2 += "ncRNA_class=RNase_P_RNA;product=%s;" % REFERENCE_FEATURES[feature_no][4]
            
                feature.append(attribute_gene)
                feature.append(attribute2)
            
                print(feature)                


            if feature[1] == "tmRNA":
                
                attribute_gene += "gbkey=Gene;gene_biotype=tmRNA;"
                
                # attribute_gene: ID, locus_tag, gene
                # attribute_rRNA: ID, locus_tag, product
            
                if REFERENCE_FEATURES[feature_no][2] != '':
                    attribute_gene += "gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                
                if REFERENCE_FEATURES[feature_no][4] != '':
                    attribute2 += "product=%s;" % REFERENCE_FEATURES[feature_no][4]
            
                feature.append(attribute_gene)
                feature.append(attribute2)
            
                print(feature)                        


#CP008699.1	Genbank	gene	46580	46652	.	+	.	ID=gene58;Name=HCTETULN_059;gbkey=Gene;gene_biotype=tRNA;locus_tag=HCTETULN_059
#CP008699.1	Genbank	tRNA	46580	46652	.	+	.	ID=rna4;Parent=gene58;Note=COORDINATES: profile:Aragorn:1.2.34%2CtRNA-Trp(tca);gbkey=tRNA;product=tRNA-Trp
                
            
            if feature[1] == "tRNA":
            
       
                attribute2 += "gbkey=tRNA;"
                attribute_gene += "gbkey=Gene;gene_biotype=tRNA;"
            
                # attribute_gene: ID, locus_tag, gene, note (like tRNA-Gly, putative pseudogene)
                # attribute_tRNA: ID, locus_tag, product
            
                attribute_gene += "gene=%s;" % REFERENCE_FEATURES[feature_no][2]
                
                if REFERENCE_FEATURES[feature_no][2] == 'Pseudo':
                    attribute_gene += "note=%s, putative pseudogene" % REFERENCE_FEATURES[feature_no][4]

                else:
                    attribute_gene += "note=%s;" % REFERENCE_FEATURES[feature_no][4]
                
                attribute2 += "product=%s;" % REFERENCE_FEATURES[feature_no][4]

                attribute2 += "note=inference:tRNAscan-SE v.2.0.5"
            
                feature.append(attribute_gene)
                feature.append(attribute2)
            
                print(feature)
                
                
        print("\n")
                        
            
    print("######################## Block 19 executed successfully! ########################\n\n")








#### Block 20. Exports annotations as GFF
# Requires: all_genes_in_the_genome
# Requires: genome_list
### For each genome:
###     open([genome].gff)
###     print("##gff-version 3")
###     print("##sequence-region", genome, 1, len(genome_sequence), sep=" ")
###     print
###     print
###     print
###     print
###     print("##FASTA", "\n", ">genome", "\n", "genome_sequence")
###     close(GFF)

if "B20" in Blocks_to_run:
    print("\n\n######################## Executing Block 20 ####################################")
    print("Exports annotations in GFF3 format. Output directory - %s/annotation/GFFs/\n" % work_dir)
    
    
    if not glob.glob("%s/annotation/" % work_dir):
        os.system("mkdir %s/annotation/" % work_dir)
    if not glob.glob("%s/annotation/GFFs" % work_dir):
        os.system("mkdir %s/annotation/GFFs" % work_dir)


    for g in range(len(genome_list)):                  ### g - genome coords in genome_list; in Prot_seq_list, this needs to be +2'ed
        print("Exports annotation for genome %s.........." % genome_list[g], end = '')
        GFF = open("%s/annotation/GFFs/%s.gff" %  (work_dir, genome_list[g]), 'w')
        GENOME = ImportFasta("%s%s.fasta" % (genomes_for_annotation_dir, genome_list[g]))

        print("##gff-version 3", file = GFF)
        print("##sequence-region ", genome_list[g], 1, len(GENOME[0][1]), sep = ' ', file = GFF)
    
        for feature in all_gene_table[g]:      			### or essentially, in all_genes_in_the_genome

############ 0               1                 2                 3                       4          5           6                   7              8
# feature: [ID; Protein/rRNA/tRNA/other; gene_name; Functional/Putative/Pseudogene; gene_start; gene_end; gene_orientation; attributes_gene; attributes_other]

            #Column 1: "seqid" = Genome_name = genome_list[g]
            Col1 = genome_list[g]
        
            #Column 2: "source" = "hmmsearch"
            Col2 = "HMMER_3.1b2"
        
            # Column 3: "type" = "gene, rRNA or CDS"
            # leaving blank for now! 
        
            #Column 4: start
            Col4 = feature[4]
        
            #Column 5: end
            Col5 = feature[5]
    
            #Column 6: "score" = "."
            Col6 = "."
        
            #Column 7: "strand"
            if feature[6] == "F":
                Col7 = "+"
            elif feature[6] == "R":
                Col7 = "-"
            else:
                Col7 = "WTF?"
        
            #Column 8: "phase" = 0
            Col8 = 0
        
            #Column 9: "attributes" - col 7 & col 8 in gene!


            if (feature[1] == "Protein" and feature[3] == "Functional") or (feature[1] == "Protein" and feature[3] == "Putative"):
                print(Col1, "HMMER_3.1b2", "gene", Col4, Col5, Col6, Col7, Col8, feature[7], sep = "\t", file = GFF)
                print(Col1, "HMMER_3.1b2", "CDS", Col4, Col5, Col6, Col7, Col8, feature[8], sep = "\t", file = GFF)

            if feature[1] == "Protein" and feature[3] == "Pseudogene":
                print(Col1, "HMMER_3.1b2", "pseudogene", Col4, Col5, Col6, Col7, Col8, feature[7], sep = "\t", file = GFF)
                
            if feature[1] == "rRNA":
                print(Col1, "rnammer_1.2", "gene", Col4, Col5, Col6, Col7, Col8, feature[7], sep = "\t", file = GFF)
                print(Col1, "rnammer_1.2", "rRNA", Col4, Col5, Col6, Col7, Col8, feature[8], sep = "\t", file = GFF)

            if feature[1] == "tRNA":
                print(Col1, "tRNAscan-SE_1.23", "gene", Col4, Col5, Col6, Col7, Col8, feature[7], sep = "\t", file = GFF)
                print(Col1, "tRNAscan-SE_1.23", "tRNA", Col4, Col5, Col6, Col7, Col8, feature[8], sep = "\t", file = GFF)

            if feature[1] == "tmRNA":
                print(Col1, "HMMER_3.1b2", "gene", Col4, Col5, Col6, Col7, Col8, feature[7], sep = "\t", file = GFF)
                print(Col1, "HMMER_3.1b2", "tmRNA", Col4, Col5, Col6, Col7, Col8, feature[8], sep = "\t", file = GFF)

            if feature[1] == "ncRNA":
                print(Col1, "HMMER_3.1b2", "gene", Col4, Col5, Col6, Col7, Col8, feature[7], sep = "\t", file = GFF)
                print(Col1, "HMMER_3.1b2", "ncRNA", Col4, Col5, Col6, Col7, Col8, feature[8], sep = "\t", file = GFF)




        print("##FASTA\n>", genome_list[g], "\n", GENOME[0][1], sep = "", file = GFF)
        
        GFF.close()
        
        print("DONE!")
        
    print("\n######################## Block 20 executed successfully! ########################\n\n")






#### Block 21. Checking whether in any of the genomes there are large chunks without annotation. If so, they are listed



if "B21" in Blocks_to_run:
    print("\n\n######################## Executing Block 21 ####################################")
    print("Checking whether in any of the genomes there are large regions without annotation, and ORFs within these regions")
    
    # Creating an output folder for unannotated ORFs ...
    if not glob.glob("%s/annotation/files/unannotated_ORFs" % work_dir):
        os.system("mkdir %s/annotation/files/unannotated_ORFs" % work_dir)


    
    for genome_no in range(len(genome_list)):
        genome = genome_list[genome_no]
        print("\nNow checking for annotation gaps, and ORFs spanning over %s nucleotides within these gaps, in the genome " % ORF_min_length, end = "")
        
        GENOME = ImportFasta("%s%s.fasta" % (genomes_for_annotation_dir, genome))
        genome_length = len(GENOME[0][1])
        print(GENOME[0][0], "\t", GENOME[0][1][0:24], "...\tlength ", genome_length, sep = "")
        
        
        ### Reading in ORFs for that genome ...
        ORFs = ImportFasta("%s/annotation/files/genomes/%s/%s_orfs.fasta" % (work_dir, genome, genome))
        """ This is how the ORFs table should look:
        ORFs = [["VFAGAHAV_6538 [973 - 941] (REVERSE SENSE)", "KTEYVCRLLLF"],
        ["VFAGAHAV_1294 [49392 - 49457]", "LYGFLFKYTQKFKNLQKEHYFF"]]
        """
        
        pos_ORFs = []
        for ORF in ORFs:
            ORF_start = int(ORF[0].strip().split()[1][1:])
            ORF_end = int(ORF[0].strip().split()[3][:-1])
            pos_start = min(ORF_start, ORF_end)
            pos_end = max(ORF_start, ORF_end)
            pos_ORFs.append([ORF[0], pos_start, pos_end, ORF[1]])

        """ This is how the pos_ORFs table should look:
        pos_ORFs =
        [['VFAGAHAV_6538 [973 - 941] (REVERSE SENSE)', 941, 973, 'KTEYVCRLLLF'],
        ['VFAGAHAV_1294 [49392 - 49457]', 49392, 49457, 'LYGFLFKYTQKFKNLQKEHYFF']]
        """

        ### Creating a table where for every base in the genome, there is information on whether anything is annotated there:
        genome_contents = []
        for base in range(genome_length):
            genome_contents.append(0)
        for feature in all_gene_table[genome_no]:
            for base in range(feature[4], feature[5]):
                genome_contents[base] = 1
        
        ### Creating a file where longer ORFs within unannotated regions will be output:
        unannotated_ORFs = open("%s/annotation/files/unannotated_ORFs/%s_unannotated_orfs.fasta" % (work_dir, genome), "w")

        
        ### Now, identifying gaps
        ### Starting at position 0, I am checking 
        pos = 0
        gap_st = 0
        unannotated_len = 0
        while pos < genome_length:
            if genome_contents[pos] == 0:
                if unannotated_len == 0:
                   gap_st = pos
                unannotated_len += 1
            elif genome_contents[pos] == 1:
                if unannotated_len > gap_length_nucleotide_threshold:
                    gap_start = gap_st+1
                    gap_end = pos
                    print("Gap_start: %d,\tGap_end: %d,\tGap_length: %d" % (gap_start, gap_end, unannotated_len))
                    
                    #### Now, looping through the list of ORFs to verify if they may fit within the gap...
                    for ORF in pos_ORFs:
                        if ORF[1] in range (gap_start-ORF_overlap_limit, gap_end):
                            if ORF[2] in range (gap_start, gap_end+ORF_overlap_limit):
                                if ORF[2]-ORF[1] > ORF_min_length:
                                    print("              Unannotated ORF: %s, length %s" % (ORF[0], ORF[2]-ORF[1]))
                                    print(">%s\n%s" % (ORF[0], ORF[3]), file = unannotated_ORFs)
                
                unannotated_len = 0
                gap_st = 0
            pos += 1
            
        ### Now, checking the possible gap between the final annotation and the contig end --- same procedure as above ...
        if unannotated_len > gap_length_nucleotide_threshold:
            
            gap_start = gap_st+1
            gap_end = pos
            print("Gap_start: %d,\tGap_end: %d,\tGap_length: %d" % (gap_start, gap_end, unannotated_len))
                    
            #### Looping through the list of ORFs to verify if they may fit within the gap...
            for ORF in pos_ORFs:
                if ORF[1] in range (gap_start-ORF_overlap_limit, gap_end):
                    if ORF[2] in range (gap_start, gap_end+ORF_overlap_limit):
                        if ORF[2]-ORF[1] > ORF_min_length:
                            print("    Unannotated ORF in the gap: %s, length %s\n" % (ORF[0], ORF[2]-ORF[1]))
                            print(">%s\n%s" % (ORF[0], ORF[3]), file = unannotated_ORFs)
            
        print("")
        
        unannotated_ORFs.close()
        

    print("\n######################## Block 21 executed successfully! ########################\n\n")

import datetime

current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print("\nScript completed on:", current_date)
end_time = time.time()
elapsed_time_seconds = end_time - start_time
elapsed_time_minutes = elapsed_time_seconds / 60
print(f"The script took {elapsed_time_seconds} seconds ({elapsed_time_minutes} minutes) to run.")
