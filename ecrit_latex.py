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

init_tex='\\documentclass[10pt]{article}\n\usepackage[utf8]{inputenc}\n\usepackage{longtable}\n\usepackage[margin=2cm]{geometry}\n\\begin{document}\n'
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
