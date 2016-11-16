
:: echo give me the name of the tex file
::set /p filename=(without the .tex): 
color 0F
: here we assume the filename
::set filename=dissertation
set filename=dissertation_section

: note that latex needs to be run twice in order for bibliography to be recompiled

latex -enable-installer -time-statistics %filename%.tex
:: another arg, not in use:
::-dont-parse-first-line

:pause

ping -n 2 127.0.0.1

color 02

latex -enable-installer -halt-on-error -time-statistics %filename%.tex

bibtex %filename%
makeindex %filename%
makeindex %filename%.glo -s %filename%.ist -t %filename%.glg -o %filename%.gls

ping -n 2 127.0.0.1

latex -enable-installer -halt-on-error -time-statistics %filename%.tex

:ping -n 2 127.0.0.1

: the following doesn't work with .eps pictures:
:latex -output-format=PDF %filename%.tex

:latex %filename%.tex

:pause

:dvips %filename%.dvi

:ps2pdf %filename%.ps

: output is a pdf file

:erase /q *.ps
:erase /q *.dvi
:erase /q *.aux
:erase /q *.log
:erase /q *.out
:erase /q *.toc
:erase /q *.brf
:erase /q *.lot
:erase /q *.lof
:: glossary
:erase /q *.glo 
:: index
:erase /q *.idx 
:: bibliography
:erase /q *.bbl 
:: bibliography
:erase /q *.blg
:erase /q *.ilg
 :: index
:erase /q *.ind
:erase /q *.ist
