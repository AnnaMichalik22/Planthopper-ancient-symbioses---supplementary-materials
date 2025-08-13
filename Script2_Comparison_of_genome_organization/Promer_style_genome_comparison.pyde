add_library('pdf')


#### Basic settings

taxa = "VF"    # "SM" or "VF"
ref = "OLIH"

#target_genomes_in_order = ["MEEALB", "MEETHA", "MEESKO"]
target_genomes_in_order = ["OLIH"]

win_w = 800
bottom_line = 800
scale_y = 200  # how many bases per pixel for query genomes? The greater the value, the more compacted the figure
scale_x = 70   # it's a bit complicated, but the greater the value, the more compacted the plot.
h = 8    # the height of x-axis reference boxes
H = 4    # the height of diagnal alignment boxes

Gene_colors  = {'rRNA': [160,0,0],
                'Metabolism': [250,150,0],
                'Biosynthesis': [125,155,255],
                'Genetic_information_processing': [0,200,0],
                'Other': [160,160,160]}

taxa_ref = taxa + ref
TABLE_ref = loadTable('{}.pro'.format(taxa_ref),'csv')    # genome_entry,VFOLIH,136080,
                                                          # VFOLIH,tufA,Genetic_information_processing,10,1200,+,

TABLE_query = loadTable('VFOLIH.pro','csv')    # genome_entry,VFOLIH,136080,
                                                   # VFOLIH,tufA,Genetic_information_processing,10,1200,+,

output_name = "/Users/annam/OneDrive/Dokumenty/Ania Michalik/publications/Sul_Vid/custom scripts/A script for the comparison of genome organization/VF_OLIH_vs_Meenoplidae.pdf"

size(800, 1000)
beginRecord(PDF, output_name)
background(255)




##### Parsing reference_file                                  
dic_ref = {}
for i in range(TABLE_ref.getRowCount()):
    row =  TABLE_ref.getRow(i)
    if row.getString(1) == taxa_ref:
        ref_size = row.getString(2)
        print(ref_size)
    if row.getString(0) == taxa_ref:
        if row.getString(2) == 'rRNA':
            gene_len = float(row.getString(4)) - float(row.getString(3))
            if gene_len > 1800:
                gene_name = '23s_rRNA'
            elif gene_len > 1000:
                gene_name = '16s_rRNA'
            else:
                gene_name = '5s_rRNA'
            dic_ref[gene_name] = row.getString(2) + ',' + row.getString(3) + ',' + row.getString(4) + ',' + row.getString(5)
        else:
            dic_ref[row.getString(1)] = row.getString(2) + ',' + row.getString(3) + ',' + row.getString(4) + ',' + row.getString(5)
        
        ### this is how dic_ref looks like:    
        ### {u'rpmF': u'Genetic_information_processing,121474,121653,+', u'def': u'Genetic_information_processing,59879,60370,+'



##### Parsing query file - getting the query genome lengths ...
num_genome = len(target_genomes_in_order)
combined_query_sizes = 0        ### total size of queries
sizes_in_order = []             ### the list of sizes of subsequent genomes
y_start_of_each_genome = {}  

for i in range(TABLE_query.getRowCount()):
    row =  TABLE_query.getRow(i)
    if row.getString(0) == "genome_entry" and row.getString(1)[2:] in target_genomes_in_order:     # genome_entry,VFOLIH,136080,
        y_start_of_each_genome[row.getString(1)] = combined_query_sizes
        combined_query_sizes += row.getFloat(2)
        sizes_in_order.append(row.getFloat(2))


print(sizes_in_order)
print(combined_query_sizes)


### Draw main frame
fill(0,0,0,0)
strokeWeight(2)
stroke(0)
rect(scale_x, bottom_line, win_w-scale_x*2, -combined_query_sizes/scale_y-h)   



######## Drawing vertical lines for the genomes ...
bp_w = float(win_w-scale_x*2)/float(ref_size)
top_position = bottom_line-combined_query_sizes/scale_y-h

ver_line_num = int(float(ref_size)/5000) + 1
f = createFont("Arial",16)
strokeWeight(1)
for i in range(ver_line_num):
    stroke(150)
    line(scale_x+i*5000*bp_w, top_position, scale_x+i*5000*bp_w, bottom_line+20)    # the line extend 20 at the bottom
    textFont(f,12)
    fill(0)
    textAlign(CENTER)
    text(str(i*10)+"k",scale_x+i*10000*bp_w, bottom_line+35)   



### Draw REFERENCE gene boxes --- at the level (bottom_line)
strokeWeight(1)
stroke(0)
for item in dic_ref.keys():
    feature = dic_ref[item].split(',')    # [category, start, end, +/-]
    COLOR = Gene_colors[feature[0]]
    fill(COLOR[0],COLOR[1],COLOR[2])
    Xs = float(feature[1])
    Xe = float(feature[2])
    D = feature[3]
    
    if D == '+':
        rect(scale_x+Xs*bp_w, bottom_line-h, (Xe-Xs)*bp_w, h)    # Going from bottom_line-h to bottom_line
    if D == '-':
        rect(scale_x+Xs*bp_w, bottom_line, (Xe-Xs)*bp_w, h)      # Going from bottom_line to bottom_line+h


### draw horizontal lines separating the genomes
stroke(0)
strokeWeight(1)
Y_line = bottom_line-h
#line(scale_x, Y_line, win_w-scale_x,Y_line)

for i in range(len(sizes_in_order)):
        Y_line -= sizes_in_order[i]/scale_y
        line(scale_x, Y_line, win_w-scale_x,Y_line)


### Draw different genes
strokeWeight(0)
Y_start = bottom_line - h    # avoid overlaps at the bottom
for i in range(TABLE_query.getRowCount()):
    row = TABLE_query.getRow(i)   # e.g., "VFOLIH,tufA,Genetic_information_processing,10,1200,+,"
    if row.getString(0)[2:] in target_genomes_in_order:
        sample = row.getString(0)    # VFOLIH
        ###pos = target_genomes_in_order.index(sample)
        gene = row.getString(1)          # tufA
        category = row.getString(2)      # Genetic_information_processing
        Qs = row.getFloat(3)             # 10
        Qe = row.getFloat(4)             # 1200
        Qd = row.getString(5)            # +
        if category == 'rRNA':           #### Quite rough gene categorization, since the names are missing from gffs ...
            gene_len = Qe - Qs
            if gene_len > 1800:
                gene = '23s_rRNA'
            elif gene_len > 1000:
                gene = '16s_rRNA'
            else:
                gene = '5s_rRNA'

        Y_value = Y_start - y_start_of_each_genome[sample]/scale_y
        if gene in dic_ref.keys():    # necessary check! 
            feature_ref = dic_ref[gene].split(',')    # [category, start, end, +/-]
            COLOR = Gene_colors[feature_ref[0]]
            fill(COLOR[0],COLOR[1],COLOR[2])
            Xs = float(feature_ref[1])
            Xe = float(feature_ref[2])
            Rd = feature_ref[3]
        
            Xs_quad = scale_x+Xs*bp_w
            Xe_quad = scale_x+Xe*bp_w
            Qs_quad = Y_value-Qs/scale_y
            Qe_quad = Y_value-Qe/scale_y

            if Qd == Rd:
                quad(Xs_quad, Qs_quad+H, Xs_quad, Qs_quad-H, Xe_quad, Qe_quad-H, Xe_quad, Qe_quad+H)
            if Qd != Rd:
                quad(Xs_quad, Qe_quad+H, Xs_quad, Qe_quad-H, Xe_quad, Qs_quad-H, Xe_quad, Qs_quad+H)



### Draw genome names
textFont(f,12)
fill(0)
textAlign(LEFT)
for genome in target_genomes_in_order:
    name = taxa+genome
    
    text_y = bottom_line - h - y_start_of_each_genome[name]/scale_y
    text(name,scale_x-65, text_y-30) 




### Draw LEGEND
textFont(f,16)
fill(0)
textAlign(LEFT)
text("Annotated {} genome ({}bp)".format(taxa_ref, ref_size), scale_x, bottom_line+60)


textFont(f,12)
legend_x = scale_x
for i in Gene_colors.keys():
    COLOR = Gene_colors[i]
    letter = float(len(i))
    fill(COLOR[0],COLOR[1],COLOR[2])
    rect(legend_x, bottom_line+90, 20, 20)
    text(i, legend_x+30, bottom_line+110)
    legend_x += sqrt(letter)*42


            
endRecord()
