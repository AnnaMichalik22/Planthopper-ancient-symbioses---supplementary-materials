add_library('pdf')
output_name = "symbio_genes.pdf"

Species_list = "Species_list.txt"
Gene_list = "genes_list.txt"
Symb_list = "Symbionts_list.txt"
Symb_color_list = "Symbionts_RGB_colors.txt"
Data_all = "Data_all.txt"

Col_width = 20
Row_width = 20
Start_X = 100
Start_Y = 100
pie_scaling_factor = 8

size(1800, 1800)
beginRecord(PDF, output_name)
background(255)


def create_species_dictionary(filepath): # script reads in a text file, with one word per line, and makes it into a dictionary where the contents of every new line is a key, and the value is a consecutive integer
    list_dict = {}
    with open(filepath, 'r') as file:
        count = 0
        for line in file:
            word = line.strip()
            list_dict[word] = count
            count += 1
    return(list_dict)

def create_column_dictionary(filepath): # script reads in a text file, with four tab-delimited columns per line, and makes it into a dictionary where the first column is a key, and the value is list of three consecutive columns (integers)
    column_dict = {}
    with open(filepath, 'r') as file:
        #count = 0
        for line in file:
            parts = line.strip().split('\t')  # Split the line by tabs
            if len(parts) == 4:
                dict_key = parts[0]
                dict_values = [int(parts[1]), int(parts[2]), int(parts[3])] # convert to int
                column_dict[dict_key] = dict_values
                #column_list.append(dict_key)
                #key_dict[dict_key] = count
            #count += 1
                
    return(column_dict)

def import_vis_data(filepath): # script reads in a text file, with three tab-delimited columns per line, and makes it into a list
    data_list = []
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')  # Split the line by tabs
            if len(parts) == 3:
                data_row = [parts[0],parts[1],parts[2]]                
                data_list.append(data_row)

    return(data_list)


def draw_pie_chart(data, x, y):
    """Draws a single pie chart with scaled diameter and values only."""
    
    Symbs = data.strip().strip(",").split(",")
    total = len(Symbs)

    if total == 0:
        return
    
    diameter = sqrt(float(total))*pie_scaling_factor
    
    angle_per_category = 360 / total

    if total == 1:
        fill(color(colors[Symbs[0]][0], colors[Symbs[0]][1], colors[Symbs[0]][2]))
        circle(x, y, diameter)
        
    start_angle = 0
    for no in range(len(Symbs)):
        fill(color(colors[Symbs[no]][0], colors[Symbs[no]][1], colors[Symbs[no]][2]))
        arc(x, y, diameter, diameter, radians(start_angle), radians(start_angle + angle_per_category))
        start_angle += angle_per_category





Species_dict = create_species_dictionary(Species_list)
Gene_dict = create_species_dictionary(Gene_list)
Symb_dict = create_species_dictionary(Symb_list)
colors = create_column_dictionary(Symb_color_list)
Data_all_list = import_vis_data(Data_all)


if Gene_dict:
    z = 0
    for item in Gene_dict:
        print(z, item)
        z += 1

print(len(Gene_dict))


##### Draw grid
Species_no = len(Species_dict)
Gene_no = len(Gene_dict)

# Draw horizontal lines
stroke(0)
strokeWeight(1)
for no in range(Species_no+1):
    line(Start_X-10, Start_Y+no*Row_width, Start_X+Gene_no*Col_width, Start_Y+no*Row_width)

# Draw vertical lines
stroke(0)
strokeWeight(1)
for no in range(Gene_no+1):
    print(no)
    line(Start_X+no*Col_width, Start_Y-10, Start_X+no*Col_width, Start_Y+Species_no*Row_width)

# Write_genome_names
f = createFont("Arial",10)
textAlign(RIGHT)
strokeWeight(1)
fill(0)
for name in Species_dict.keys():
    text(name, Start_X-5, Start_Y+Row_width-4+int(Species_dict[name])*Row_width)

# Draw circles
for gene_record in Data_all_list:
    if len(gene_record) == 3:
        Species = gene_record[0]
        Gene = gene_record[1]
        data = gene_record[2]
        if Species in Species_dict.keys():
            x = Start_X + Gene_dict[Gene]*Col_width+Col_width/2
            y = Start_Y + Species_dict[Species]*Row_width+Row_width/2
            draw_pie_chart(data, x, y)

# Make legend
Legend_top_Y = Start_Y+Species_no*Row_width+50
for Symb in Symb_dict.keys():
    print(Symb)
    fill(color(colors[Symb][0], colors[Symb][1], colors[Symb][2]))
    Legend_Y = Legend_top_Y + int(Symb_dict[Symb])*Row_width
    circle(100, Legend_Y, 10)
    textAlign(LEFT)
    strokeWeight(1)
    fill(0)
    text(Symb, 120, Legend_Y)

# Write_gene_names
f = createFont("Arial",10)
rotate(radians(270))
translate(-190, 0) 
for name in Gene_dict.keys():
    text(name, Start_X-5, Start_Y+Row_width-4+int(Gene_dict[name])*Row_width)


endRecord()
