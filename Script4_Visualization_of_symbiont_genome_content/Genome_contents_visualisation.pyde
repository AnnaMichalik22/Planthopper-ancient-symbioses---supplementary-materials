add_library('pdf')

### defining the size of your "drawing board"
size(3800, 2400)

### stating where input data table is, and where the drawing should be output
Data_Table = loadTable("C:/Users/annam/OneDrive/Dokumenty/Ania Michalik/publications/Sul_Vid/custom scripts/A script for the visualization of symbiont genome content/Sulcia_genome_contents_table.csv")
"""
gene_tag,gene_type,funct_category,TETCRI,PL623,DB114
hisA,Protein,01_Biosynthesis_His,Functional,Functional,Functional
hisB,Protein,01_Biosynthesis_His,Functional,Functional,Absent
"""

beginRecord(PDF, "C:/Users/annam/OneDrive/Dokumenty/Ania Michalik/publications/Sul_Vid/custom scripts/A script for the visualization of symbiont genome content/Sulcia_genome_contents.pdf")



### setting background color
background(255)


### reading the first row of the table, typing column labels (genome names)
row = Data_Table.getRow(0) #0 - ID; 1 - gene_name; 2 - TETULN; 3,4 - TETUND; 5-10 - TETCHI
f = createFont("Arial",20)
textAlign(RIGHT)
textFont(f,20)
fill(0)

### Printing genome ID's --- names of columns from 4 onwards:
for k in range(3, Data_Table.getColumnCount()):
    text(row.getString(k),200,43+40*k)

### For rows from #2 onwards, for columns from #4 onwards, set fill color based on cell contents, and then draw a rectangle with that fill color
for i in range(1,Data_Table.getRowCount()):
    row = Data_Table.getRow(i) 
    if row.getString(0) != '':    
        for k in range(3,Data_Table.getColumnCount()):
            if row.getString(k) != '': 
                if row.getString(k) == 'Absent':  
                    fill(255)
                elif row.getString(k) == 'Pseudogene':  
                    fill(150,150,150)         
                elif row.getString(k) == 'Putative':  
                    fill(2,2,0)    
                elif row.getString(k) == 'Functional':  
                    fill(1,1,0)    
                ### Now drawing rectangle with previously defined fill. Row number (i) is used to define starting x position; Column no (k) to define starting y position;
                ### Rectangle width and height are fixed at 15 and 35.
                rect(200+18*i, 20+40*k,15,35)

### Now, typing the gene names (first column of the table). This requires rotating the drawing canvas ... you may need to play with the "translate" values if your labels are off!
textSize(16)
rotate(PI*.5)
translate(38, -278)
for i in range(1,Data_Table.getRowCount()):
    row = Data_Table.getRow(i) #0 - ID; 1 - gene_name; 2 - TETULN; 3,4 - TETUND; 5-10 - TETCHI
    if row.getString(0) != '':
        textAlign(RIGHT)
        fill(0)
        text(row.getString(0), 92, 75-18*i)

### Should be all done :)
endRecord()

