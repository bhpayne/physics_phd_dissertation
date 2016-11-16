#!/bin/bash

#texfile="dissertation_publication"
#texfile="dissertation_section"
texfile="dissertation_report"
#texfile="dissertation"
#texfile="ptotemplate"

rm -f $texfile.{aux,bbl,blg,dvi,glg,glo,gls,idx,ist,fff,nlo,nls,ilg}
rm -f $texfile.{lof,log,lot,out,toc,pdf,ps,ttt,xdy}


latex $texfile

makeindx $texfile.idx
#makeindex $texfile.glo -s $texfile.ist -t $texfile.glg -o $texfile.gls

bibtex $texfile

latex $texfile

latex $texfile

latex $texfile

#  this chokes on big dvi
#dvipdf $texfile

#  create ps from dvi. Works for large files
dvips $texfile.dvi

sleep 3

#  ps2pdf chokes on big files
ps2pdf $texfile.ps 

#pdflatex $texfile.tex


sleep 3

#evince $texfile.ps
evince $texfile.pdf
