%%
%% This is file `mstthesis.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% theses.dtx  (with options: `mstthesis')
%% 
%% IMPORTANT NOTICE:
%% 
%% For the copyright see the source file.
%% 
%% Any modified versions of this file must be renamed
%% with new filenames distinct from mstthesis.cls.
%% 
%% For distribution of the original source see the terms
%% for copying and modification in the file theses.dtx.
%% 
%% This generated file may be distributed as long as the
%% original source files, as listed above, are part of the
%% same distribution. (The sources need not necessarily be
%% in the same archive or directory.)

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mstthesis}[2008/01/22 MO S&T Thesis Class]
\RequirePackage{ifthen}
\newboolean{timesfree}
\newboolean{timescomm}
\newboolean{double}
\newboolean{onehalf}
\DeclareOption{10pt}{%
\PassOptionsToClass{10pt}{report}%
}

\DeclareOption{11pt}{%
\PassOptionsToClass{11pt}{report}%
}

\DeclareOption{12pt}{%
\PassOptionsToClass{12pt}{report}%
}

\DeclareOption{times}{%
  \setboolean{timesfree}{false}%
}

\DeclareOption{timesmtp}{%
  \setboolean{timescomm}{true}%
}

\DeclareOption{doublespace}{%
  \setboolean{double}{true}
}

\DeclareOption{onehalfspace}{%
  \setboolean{onehalf}{true}%
}

\newif\ifchap  % true for chap option
  \chapfalse   % false by default
\DeclareOption{chap}{\chaptrue} % option to print "Chapter" at each new chapter

\ProcessOptions
\LoadClass{report}

\ifthenelse{\boolean{timescomm}}{%
  \RequirePackage[T1]{fontenc}% http://www.ctan.org/pkg/fontenc
  \RequirePackage{textcomp}% http://www.ctan.org/pkg/textcomp
  \renewcommand{\rmdefault}{ptm}%
  \RequirePackage{mtpro2}% http://www.pctex.com/mtpro2.html
}{}

\ifthenelse{\boolean{timesfree}}{%
  \RequirePackage[T1]{fontenc}%
  \RequirePackage{textcomp}%
  \renewcommand{\rmdefault}{ptm}%
  \RequirePackage{mathptmx}  %  http://www.ctan.org/pkg/mathptmx
  \RequirePackage{pifont}% http://www.ctan.org/pkg/pifont
}{}

\RequirePackage{geometry}
% Normally the following gives 1`` right and bottom margins
%\geometry{letterpaper,total={5.985in,9.5in},includehead} % page width and length
% However, due to use of \brokenpenalty this is not always obeyed
\geometry{letterpaper,total={5.985in,9.37in},includehead} % page width and length
\geometry{top=0.515in,headheight=0.125in,headsep=0.3125in,left=1.515in} % margins

\RequirePackage{setspace}
\RequirePackage{amsmath} % ,amssymb
\RequirePackage[american]{babel}
\RequirePackage[notoday,nofancy]{rcsinfo}
\RequirePackage{indentfirst}
\RequirePackage[noprefix]{nomencl}

\RequirePackage{titlesec}
\titleformat{\chapter}{\singlespacing\center\normalsize\bfseries\uppercase}%
  {\thechapter.}{.5em}{}[\doublespacing]
\titlespacing{\chapter}{0pt}{-0.5in}{0.6875in} % -0.5in is the vertical spacing of the chapter in the body
\titleformat{\section}{\singlespacing\normalsize\bfseries\uppercase}%
  {\thesection.}{.5em}{}[\doublespacing]
\titleformat{\subsection}[runin]{\normalsize\bfseries}%
  {\thesubsection.}{0.5em}{}[.\hspace{0.5em}~]
\titlespacing{\subsection}{\parindent}{\parskip}{0pt}
\titleformat{\subsubsection}[runin]{\normalsize\bfseries}%
  {\thesubsubsection.}{0.5em}{}[.\hspace{0.5em}~]
\titlespacing{\subsubsection}{\parindent}{\parskip}{0pt}

\RequirePackage{fancyhdr}
\pagestyle{fancy}
\lhead{}\chead{}\rhead{\thepage}
\lfoot{}\cfoot{}\rfoot{}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}

\RequirePackage[authoryear,sort]{natbib}


\ifthenelse{\boolean{double}}%
{\newcommand{\paperbodysettings}{\doublespacing}}%
{\newcommand{\paperbodysettings}{}}

\ifthenelse{\boolean{onehalf}}%
{\renewcommand{\paperbodysettings}{\onehalfspacing}}%
{\renewcommand{\paperbodysettings}{}}

\def\nopagenumber{%
\thispagestyle{empty}%
}

\setlength{\parindent}{0.5in} % paragraph indent

\newcommand{\ThesisAdviser}[1]{\thadviser{#1}}
\newcommand{\CoThesisAdviser}[1][false]{%
\ifthenelse{\equal{#1}{true}}{\cothadviser{#1}}{}}

\newenvironment{ThesisTitlePage}[1]{%
  \ifthenelse{\equal{#1}{ms}}{\thesistype{1}}{\thesistype{2}}}%
{\titlepage}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \def\thesistitle#1{\gdef\@thesistitle{#1}}
 \def\author#1{\gdef\@author{#1}}
 \def\degree#1{\gdef\@degree{#1}}
 \def\department#1{\gdef\@department{#1}}
 \newcount\numcomm \numcomm=1
 \def\thesistype#1{\numcomm=#1}
 \def\thadviser#1{\gdef\@thadviser{#1}}
 \def\projadviser#1{\gdef\@projadviser{#1}} \projadviser{*}
 \def\cothadviser#1{\gdef\@cothadviser{#1}} \cothadviser{*}
 \def\coprojadviser#1{\gdef\@coprojadviser{#1}} \coprojadviser{*}
 \def\doctype{\if \@projadviser *Thesis \else Project \fi}
 \def\adviser{\if \@projadviser *\@thadviser \else \@projadviser\fi}
 \def\coadviser{\if \@coprojadviser *\@cothadviser\else \@coprojadviser\fi}
 \def\memberone#1{\gdef\@memberone{#1}}
 \def\membertwo#1{\gdef\@membertwo{#1}}
 \def\memberthree#1{\gdef\@memberthree{#1}}
 \def\memberfour#1{\gdef\@memberfour{#1}}
 \def\memberfive#1{\gdef\@memberfive{#1}}
 \def\membersix#1{\gdef\@membersix{#1}}
 \def\copyrightyear#1{\gdef\@copyrightyear{#1}} \copyrightyear{\the\year}
 \def\graddate#1{\gdef\@graddate{#1}}

\newcount\numcount

\def\titlepage{%

  \ifnum\numcomm<2 \typeout{---------------Masters Thesis----------}

     \thispagestyle{empty}
     \mbox{}\vspace*{\fill}
     \begin{singlespace}
     \begin{center}
         \parindent=0pt
         \vspace*{-0.25in} % ... koh (2011/09/06) ...
        \begin{doublespace}
        {\uppercase\expandafter{\@thesistitle}}\\ [24pt]
         by \\ [24pt]
         {\uppercase\expandafter{\@author}}\\ [24pt]
         \end{doublespace}
         A \uppercase\expandafter{\doctype} \\[12pt]
         Presented to the Faculty of the Graduate School of the\\ [18pt]
         MISSOURI UNIVERSITY OF SCIENCE AND TECHNOLOGY\\ [18pt]
         In Partial Fulfillment of the Requirements for the Degree\\[24pt]
         \uppercase\expandafter{\@degree}~IN~%
         \uppercase\expandafter{\@department} \\ [24pt]
         \@graddate\\[18pt]
         Approved by\\
     \vspace*{0.375in}
     \end{center}
     \if \coadviser *  % no coadviser
        \begin{center}
        \begin{minipage}{5in}
          \begin{center}\begin{onehalfspacing}
          \adviser, Advisor\\
          \@memberone \\
          \@membertwo\\
          \end{onehalfspacing}\end{center}
        \end{minipage}
        \end{center}
     \else
        \begin{center}
        \begin{minipage}{5in}
          \begin{center}\begin{onehalfspacing}
          \adviser, Co-advisor\\
          \@memberone, Co-advisor \\
          \@membertwo\\
          \@memberthree \\
          \end{onehalfspacing}\end{center}
        \end{minipage}
        \end{center}
     \fi
     \end{singlespace}
     \vspace*{\fill}
     \newpage
     \pagestyle{myheadings}
     \pagenumbering{roman}
     \setcounter{page}{2}

  \else  \typeout{---------------Doctoral Thesis----------}

     \thispagestyle{empty}
     \mbox{}\vspace*{\fill}
     \begin{singlespace}
     \begin{center}
         \parindent=0pt
         \vspace*{-0.25in}
        \begin{doublespace}
        {\uppercase\expandafter{\@thesistitle}}\\ [20pt]
         by \\ [20pt]
         {\uppercase\expandafter{\@author}}\\ [20pt]
         \end{doublespace}
         A DISSERTATION \\[18pt]
         Presented to the Faculty of the Graduate School of the\\ [18pt]
         MISSOURI UNIVERSITY OF SCIENCE AND TECHNOLOGY\\ [18pt]
         In Partial Fulfillment of the Requirements for the Degree\\[24pt]
         \uppercase\expandafter{\@degree}\\[18pt]
         ~in~\\[18pt]%
         \uppercase\expandafter{\@department} \\[18pt]
         \@graddate\\[18pt]
         Approved by\\
     \vspace*{0.3in}
     \end{center}
     \if \coadviser *  % no coadviser
        \begin{center}
        \begin{minipage}{5in}
          \begin{center}\begin{onehalfspacing}
          \adviser, Advisor\\
          \@memberone \\
          \@membertwo\\
          \@memberthree \\
          \@memberfour \\
          \end{onehalfspacing}\end{center}
        \end{minipage}
        \end{center}
     \else
        \begin{center}
        \begin{minipage}{5in}
          \begin{center}\begin{onehalfspacing}
          \adviser, Co-advisor\\
          \@memberone, Co-advisor \\
          \@membertwo\\
          \@memberthree \\
          \@memberfour \\
          \end{onehalfspacing}\end{center}
        \end{minipage}
        \end{center}
     \fi
     \end{singlespace}
     \vspace*{\fill}
     \newpage
     \pagestyle{myheadings}
     \pagenumbering{roman}
     \setcounter{page}{2}

  \fi
  }
%%%%
 \def\abstitlepage{%
     \if@twoside\thispagestyle{empty}\cleardoublepage\fi
     \thispagestyle{empty}
     \begin{singlespace}
     \begin{center}
         \parindent=0pt
         \vglue 24pt plus 24pt minus 12pt
         {\uppercase\expandafter{\@thesistitle}}\\ [11pt]
         By \\ [11pt]
         \@author\\ [11pt]
         An Abstract of a Thesis Submitted to the Graduate\\ [7pt]
         Faculty of Rensselaer Polytechnic Institute\\ [7pt]
         in Partial Fulfillment of the\\ [7pt]
         Requirements for the Degree of\\ [7pt]
         \uppercase\expandafter{\@degree}\\ [7pt]
         Major Subject:~~\@department\\  [7pt]
         The original of the complete thesis is on file\\
         in the Rensselaer Polytechnic Institute Library\\
       \vfill
       \sigblockabs
       \vskip 36pt plus 2pt minus 12pt
        Rensselaer Polytechnic Institute\\
        Troy, New York\\ [12pt]
         \@graddate
     \end{center}
     \end{singlespace}
     \newpage
     \pagestyle{myheadings}
     \pagenumbering{roman}
     \setcounter{page}{2}}

\newcommand{\ThesisCopyrightPage}[1]{%
  \ifthenelse{\equal{#1}{true}}%
    {\copyrightpage}%
  {\thispagestyle{empty}\vspace*{1in}\vfill\eject}%
}

\newenvironment{ThesisPublicationOption}[1]{%
  \ifthenelse{\equal{#1}{ms}}%
    {\chapter*{Publication Thesis Option}%
      \addcontentsline{toc}{chapter}{PUBLICATION THESIS OPTION}}%
    {\chapter*{Publication Dissertation Option}
      \addcontentsline{toc}{chapter}{PUBLICATION DISSERTATION OPTION}}%
}{}

\newenvironment{ThesisAbstract}{%
  \chapter*{ABSTRACT}%
  \addcontentsline{toc}{chapter}{ABSTRACT}%
}{}

\newenvironment{ThesisAcknowledgment}{%
  \chapter*{ACKNOWLEDGMENTS}%
  \addcontentsline{toc}{chapter}{ACKNOWLEDGMENTS}%
}{}

\newenvironment{ThesisFrontMatter}{}{%
  \addtocontents{toc}{\protect\vspace{4ex}}%
  \addtocontents{toc}{\hbox to \textwidth{{}\hspace*{-2mm}{\bf SECTION\ \hfill}}}%
}


\newenvironment{ThesisBackMatter}{%
  \newpage%
  %\pagenumbering{arabic}%
  %\setcounter{page}{1}
  %\addtolength{\textfloatsep}{+0.28125in}%
  \addtocontents{toc}{\protect\vspace{4ex}}%
  \addtocontents{toc}{\hbox to \textwidth{{}\hspace*{-2mm}{\bf SECTION\ \hfill}}}%
}{}

\def\copyrightpage{
    \hbox{ }
    \thispagestyle{empty}
    \vfill
    \begin{center}
     Copyright \@copyrightyear \\
    \@author \\
    All Rights Reserved \\ [12pt]
    \end{center}
    \clearpage}
\renewcommand{\tableofcontents}{%
  \renewcommand{\contentsname}{TABLE OF CONTENTS}
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \chapter*{\contentsname
      \@mkboth{\uppercase{\contentsname}}{\uppercase{\contentsname}}}
  \vspace*{-0.375in} % vertical where the contents of TOC start
  \hbox to \textwidth{{} \hspace*{-2mm}{\ \hfill Page}}
  \singlespacing\setlength{\parskip}{10\p@} % 10 is lots of space, like double. % 6 or 5 is like single spacing.
  \@starttoc{toc}%
  \if@restonecol\twocolumn\fi
  }

\renewcommand{\listoffigures}{%
  \newpage%
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \addcontentsline{toc}{chapter}{LIST OF ILLUSTRATIONS}
  \renewcommand{\listfigurename}{LIST OF ILLUSTRATIONS}
  \chapter*{\listfigurename
    \@mkboth{\uppercase{\listfigurename}}%
            {\uppercase{\listfigurename}}}%
  \vspace*{-0.375in}
  \hbox to \textwidth{{} \hspace*{-2mm}{Figure \ \hfill Page}}
  \singlespacing\setlength{\parskip}{10\p@}
  \@starttoc{lof}%
  \if@restonecol\twocolumn\fi
  \setlength{\parskip}{0pt}\doublespacing% ... added by koh ...
  }

\renewcommand{\listoftables}{%
  \newpage%
  \if@twocolumn
    \@restonecoltrue\onecolumn
  \else
    \@restonecolfalse
  \fi
  \addcontentsline{toc}{chapter}{LIST OF TABLES}
  \renewcommand{\listtablename}{LIST OF TABLES}
  \chapter*{\listtablename
    \@mkboth{\uppercase{\listtablename}}{\uppercase{\listtablename}}}%
  \vspace*{-0.375in}
  \hbox to \textwidth{{} \hspace*{-2mm}{Table \ \hfill Page}}
  \singlespacing\setlength{\parskip}{10\p@}
  \@starttoc{lot}%
  \if@restonecol\twocolumn\fi
  \setlength{\parskip}{0pt}\doublespacing% ... added by koh ...
  }

\newcommand{\listofsymbols}{%
  \newpage%
  \addcontentsline{toc}{chapter}{LIST OF SYMBOLS}%
  \renewcommand{\nomname}{LIST OF SYMBOLS}%
  \printnomenclature%
}


\newcommand{\Thesistocentry}[1]{%
\addtocontents{toc}{\protect\vspace{4ex}}
\addtocontents{toc}%
{\hbox to \textwidth{{} \hspace*{-2mm}{{#1} \ \hfill}}}
}

\newenvironment{ThesisBody}{%
  \newpage%
  \pagenumbering{arabic}%
  \setcounter{page}{1}
  \addtolength{\textfloatsep}{+0.28125in}%
}{}

\newenvironment{ThesisPublications}{\Thesistocentry{\bf PAPERS}}{}

\newcommand{\ThesisListofEntry}[2]{%
\addtocontents{#1}{\protect\vspace{2ex}}
\addtocontents{#1}%
{\hbox to \textwidth{{} \hspace*{-2mm}{{#2} \ \hfill}}}
}

\newcommand{\PaperManuscript}[2]{%
  \ifthenelse{\equal{#1}{1}}
  {\chapter*{PAPER \\ \ \\ \ \\ \ \\#1. #2}} % as of 20120222, Sharon Matson at OGS says this is required
  {\chapter*{#1. #2}}
  \addcontentsline{toc}{chapter}{#1. #2}%
  \ThesisListofEntry{lof}{PAPER #1}%
  \ifthenelse{\equal{#1}{6}} % crappy hack by BHP: I only have one table, which is in paper 6. All others are empty. 
  {\ThesisListofEntry{lot}{PAPER #1}}{}
}

\def\@chapter[#1]#2{\ifnum \c@secnumdepth >\m@ne
                         \refstepcounter{chapter}%
                         \typeout{\@chapapp\space\thechapter.}%
                         \addcontentsline{toc}{chapter}%
                                   {\protect\numberline{\thechapter}#1}% chapter in TOC
                    \else
                      \addcontentsline{toc}{chapter}{#1}%
                    \fi
                    \chaptermark{#1}%
                    \ifthenelse{\value{page}=1}{\thispagestyle{empty}}{} % first page of first chapter has no page number
                    %\addtocontents{lof}{\protect\addvspace{10\p@}}%
                    %\addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{#2}]%
                    \else
                      %\ifthenelse{\equal{#1}{CONCLUSIONS}}{ % crappy hack by BHP to have a second ``SECTION'' show up after PAPERS. Should actually be done in ThesisBackMatter
                      %  \@makechapterhead{SECTION \\ \ \\ \ \\ \ \\#2}
                      %}{
                        \@makechapterhead{#2}% this is what shows up at the top of the beginning of a chapter
                      %}
                      \@afterheading
                    \fi}

\renewcommand*{\l@chapter}[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@ % vertical skip in the TOC between chapters
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode %\bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      \uppercase{#1}\nobreak \dotfill \nobreak\hbox to\@pnumwidth{\hss #2}\par % contents of TOC
      \penalty\@highpenalty
    \endgroup
  \fi}

%% BHP messing with indentation of sections for the TOC, as requested by OGS on 20120222; I don't think these look as professional
% % % \renewcommand*\l@figure{\vskip 10pt plus 1pt minus 1pt
% % %                        \@dottedtocline{1}{0em}{2.8em}}
% % % \let\l@table\l@figure
% % % \renewcommand*\l@chapter{\pagebreak[3]\vskip 10pt plus 1pt minus 1pt
% % %                          \@dottedtocline{0}{0em}{2em}}

%%including these addresses the indentation of sections in TOC, but screws up the appendices
% 
% \renewcommand*\l@section{\vskip 6pt plus 1pt minus 1pt
%                          \@dottedtocline{1}{4em}{2em}}
% \renewcommand*\l@subsection{\vskip 4pt minus 1pt
%                          \@dottedtocline{2}{8em}{3em}}

% % % \renewcommand*\l@subsubsection{\ifnum\c@tocdepth>2 \vskip 3pt minus 1pt \fi
% % %                          \@dottedtocline{3}{10em}{2em}}
% % % \renewcommand*\l@paragraph{\ifnum\c@tocdepth>3 \vskip 2pt minus 1pt \fi
% % %                          \@dottedtocline{4}{7.2em}{3.6em}}


\newenvironment{ThesisAppendix}[1]{%
  \newpage%
  \ifthenelse{\equal{#1}{one}}%
    {\addcontentsline{toc}{chapter}{APPENDIX}%
     \addtocontents{toc}{\protect\setcounter{tocdepth}{-1}}%
     \titleformat{\chapter}{\center\normalsize\bfseries\uppercase}{}{.5em}{}%
     \titlespacing{\chapter}{0pt}{-0.5in}{0.6875in}%
     \titleformat{\section}{\normalsize\bfseries\uppercase}{}{.5em}{}%
     \titleformat{\subsection}[runin]%
        {\normalsize\bfseries}{\arabic{subsection}.}{0.5em}%
        {}[.\hspace{0.5em}~]%
     \titleformat{\subsubsection}[runin]%
       {\normalsize\bfseries}%
       {\arabic{subsection}.\arabic{subsubsection}.}%
       {0.5em}{}[.\hspace{0.5em}~]}%
    {\addtocontents{toc}{\protect\setcounter{tocdepth}{0}}%
     \titleformat{\chapter}%[block]%
       {\vspace*{\stretch{1}}\center\normalsize\bfseries\uppercase}%
       {APPENDIX \thechapter}{0em}% this goes on the page prior to the appendix
       {\\~\\\hspace{-1in}\nopagenumber}[\vspace*{\stretch{1}}\newpage] % content of appendix cover
     %\titlespacing{\chapter}{0pt}{-0.5in}{0.6875in}
     \titleformat{\section}{\normalsize\bfseries\uppercase}{}{.5em}{}%
     \titleformat{\subsection}[runin]%
       {\normalsize\bfseries}{\arabic{subsection}.}{0.5em}%
       {}[.\hspace{0.5em}~]%
     \titleformat{\subsubsection}[runin]%
       {\normalsize\bfseries}%
       {\arabic{subsection}.\arabic{subsubsection}.}{0.5em}%
       {}[.\hspace{0.5em}~]%
     \Thesistocentry{\bf APPENDICES}
     \appendix}%
}{\addtocontents{toc}{\protect\setcounter{tocdepth}{4}}}


\renewenvironment{thebibliography}[1]
     {\chapter*{\MakeUppercase{BIBLIOGRAPHY}}% ... added by koh ...
      \addcontentsline{toc}{chapter}{BIBLIOGRAPHY}
      \@mkboth{\MakeUppercase\bibname}{\MakeUppercase\bibname}%
      \singlespacing\setlength{\parskip}{2pt} % ... added by koh ...
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist
      \doublespacing % ... added by koh, 2011/10/20 ...
}

\newenvironment{Vita}{%
  \chapter*{Vita}%
  \addcontentsline{toc}{chapter}{VITA}}{}

\renewcommand{\ps@plain}{%
  \renewcommand{\@oddhead}{\hfill\textrm{\thepage}}%
  \renewcommand{\@evenhead}{\@oddhead}%
  \renewcommand{\@evenfoot}{}%
  \renewcommand{\@oddfoot}{}}

\newcommand{\psplainreset}{\renewcommand{\ps@plain}{%
  \renewcommand{\@oddhead}{\hfill\textrm{\thepage}}%
  \renewcommand{\@evenhead}{\@oddhead}%
  \renewcommand{\@evenfoot}{}%
  \renewcommand{\@oddfoot}{}}
}

\renewcommand\@dotsep{4.5}

\endinput
%%
%% End of file `mstthesis.cls'.
