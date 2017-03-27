#-*- coding: utf-8 -*-

col_tab='p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}'
width=len(col_tab)/6
#print(width)


with open('promoters_dofs.fasta','r') as f1:
    list_lines=list()
    for line in f1:
        line=line.strip("\n")
        list_lines.append(line)
        #print(list_lines[0])

seq=list_lines[1]

with open('Interdistances_ARF2_DOF16.csv','r') as f1:
    list_sites=list()
    for n,site in enumerate(f1):
        site=site.strip("\n")
        site=site.split("\t")
        if (n!=0):
            list_sites.append(site)
        
orientation=list()
i=0
while (i<len(seq)):
    orientation.append(False)
    i=i+1

orientation2 = list (orientation)


for elt in list_sites:
    i=0
    while(i < 6):
        if elt[0] == 'DR':
            orientation[int(elt[2])-1 + i] =  True
            orientation[int(elt[2]) + int(elt[1]) + i -1] =  True
        if elt[0] == 'DR rev':
            orientation2[int(elt[2])-1 + i] =  True
            orientation2[int(elt[2]) + int(elt[1]) + i -1] =  True
        if elt[0] == 'ER':
            orientation[int(elt[2])-1 + i] =  True
            orientation2[int(elt[2]) + int(elt[1]) + i -1] =  True
        if elt[0] == 'IR':
            orientation2[int(elt[2])-1 + i] =  True
            orientation2[int(elt[2]) + int(elt[1]) + i -1] =  True
        i = i+1
        
couleur=list()

i=0
while (i<len(seq)):
    couleur.append(0)
    i=i+1

i=0
while (i < len(couleur)):
    if(orientation[i]):
        couleur[i] += 1
    if(orientation2[i]):
        couleur[i] += 2
    i+=1



init_tex='\\documentclass[10pt]{article}\n\usepackage[utf8]{inputenc}\n\usepackage{longtable}\n\usepackage{colortbl}\n\usepackage[margin=2cm]{geometry}\n\\begin{document}\n'
init_tab='\\begin{longtable}{'+col_tab+'}\n'
end_tab='\\end{longtable}\n'
end_tex='\\end{document}\n'
#print(init_tab)


i=0
j=0
texte = str()
texte2 = str()
while( i < len(seq)):
    texte = texte + seq[i]
    texte2 = texte2 + '.' 
    i = i+1
    j = j+1
    if(j == width):
        j=0
        texte = texte + ' \\\\\n'
        texte = texte + texte2 + ' \\\\\n'
        texte2 = str()
    else :
        texte = texte + ' & '
        texte2 = texte2 + ' & '


with open('pDOF16.tex','w') as f1:
    f1.write(init_tex)
    f1.write(init_tab)
    f1.write(texte)
    f1.write(end_tab)
    f1.write(end_tex)
