% Minimal local stand-in for elsarticle, for OFFLINE preview compilation only.
\ProvidesClass{elsarticle}[2026/09/27 offline stub]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[11pt,a4paper]{article}
\RequirePackage{geometry}\geometry{margin=1in}
\RequirePackage[numbers]{natbib}
\newcommand{\journal}[1]{}
\newenvironment{frontmatter}{}{\maketitle}
\newcommand{\ead}[2][]{}
\newcommand{\address}[2][]{\par\small #2\par\normalsize}
\newcommand{\fntext}[2]{}
\newcommand{\tnotetext}[2][]{}
\newcommand{\corref}[1]{}
\newcommand{\fnref}[1]{}
\newcommand{\tnoteref}[1]{}
\newenvironment{keyword}{\par\medskip\noindent\textbf{Keywords: }}{\par}
\providecommand{\sep}{; }
\bibliographystyle{plainnat}
\newcommand{\cortext}[2][]{}
\newcommand{\nonumnote}[1]{}
