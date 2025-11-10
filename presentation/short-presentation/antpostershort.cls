% This material is declared a work of the United States Government and is not
% subject to copyright protection in the United States.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and sell copies
% of the Software, and to permit persons to whom the Software is furnished to do
% so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% ------------------------------------------------------------------------------
% Author: David Woodburn
% Description: business presentation template
% ------------------------------------------------------------------------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{antposter}[2024/10/01 ANT Center Poster]
\def\@version{20241001}

% Process options.
\newif\ifdefineul\defineultrue
\DeclareOption{noul}{\defineulfalse}

\newif\ifcompassonly\compassonlyfalse
\DeclareOption{compassonly}{\compassonlytrue}
\newif\ifafit\afitfalse
\DeclareOption{afit}{\afittrue}

\def\@colortheme{7} % default color theme
\DeclareOption{blue}{\def\@colortheme{0}}
\DeclareOption{azure}{\def\@colortheme{1}}
\DeclareOption{cyan}{\def\@colortheme{2}}
\DeclareOption{green}{\def\@colortheme{3}}
\DeclareOption{lime}{\def\@colortheme{4}}
\DeclareOption{yellow}{\def\@colortheme{5}}
\DeclareOption{orange}{\def\@colortheme{6}}
\DeclareOption{red}{\def\@colortheme{7}}
\DeclareOption{magenta}{\def\@colortheme{8}}
\DeclareOption{purple}{\def\@colortheme{9}}

\def\@layouttheme{1} % default layout theme
\DeclareOption{twopart}{\def\@layouttheme{1}}

\DeclareOption*{\PackageWarning{antposter}{Unknown '\CurrentOption'}}
\ProcessOptions\relax

% Load the base class.
\ifnum\@layouttheme=0 % center
    \LoadClass[twocolumn, 10pt, fleqn]{article}
\fi
\ifnum\@layouttheme=1 % left
    \LoadClass[onecolumn, 12pt, fleqn]{article}
\fi

% Packages
\RequirePackage[T1]{fontenc} % use 8-bit fonts
\RequirePackage{geometry} % page dimensions and margins
\RequirePackage{pgfplots} % includes tikz, xcolor, graphicx
\RequirePackage{amsmath, amsfonts, amssymb} % improved math
\RequirePackage{bm} % bold math (must be after amsmath)
\RequirePackage{cancel} % strike-through in equations
\RequirePackage{upquote} % upright single quotes in verbatim
\RequirePackage{array} % custom column alignment types
\RequirePackage[absolute]{textpos} % free placement of objects
\RequirePackage{enumitem} % custom lists and list spacing
\RequirePackage{qrcode} % QR code generation
\RequirePackage{hyperref} % hyperlinks

\pgfplotsset{compat=newest}

% ----------------------
% Added content features
% ----------------------

% Standard math macros
\ifdefineul
    \newcommand{\ul}[1] % vectors
        {{}\mkern1mu\underline{\mkern-1mu#1\mkern-1mu}\mkern1mu}
\fi
\DeclareSymbolFont{euler}{U}{eur}{m}{n}     % symbol font
\DeclareMathSymbol{\PI}{\mathord}{euler}{25} % up pi

% Convenience macros
\DeclareMathOperator{\sgn}{sgn}             % sign function
\DeclareMathOperator{\tr}{tr}               % trace of a matrix
\DeclareMathOperator{\adj}{adj}             % adjugate of a matrix
\DeclareMathOperator*{\argmax}{arg\,max}    % argmax operator
\DeclareMathOperator*{\argmin}{arg\,min}    % argmin operator
\newcommand{\abs}[1]{\left\lvert#1\right\rvert} % absolute value
\newcommand{\norm}[1]{\left\lVert#1\right\rVert} % vector or matrix norm
\newcommand{\ssm}[1]{[#1]_{\rule{0pt}{1.38ex}\times}} % skew symmetric matrix
\newcommand{\diff}[2]{\frac{\mathrm{d}{#1}}{\mathrm{d}{#2}}} % derivative
\newcommand{\piff}[2]{\frac{\partial{#1}}{\partial{#2}}} % partial derivative
\newcommand{\EE}[2]{\ensuremath{{{#1}\!\times\!10^{#2}}}} % scientific notation
\newcommand{\rom}[1]{\uppercase\expandafter{\romannumeral #1\relax}}

% Upright math letters
\newcommand{\upc}[1][]{\mathrm{c}_{#1}}     % abbreviated cosine
\newcommand{\ups}[1][]{\mathrm{s}_{#1}}     % abbreviated sine
\newcommand{\upt}[1][]{\mathrm{t}_{#1}}     % abbreviated tangent
\newcommand{\upd}{\mathrm{d}}               % differential operator
\newcommand{\upe}{\mathrm{e}}               % Euler's number
\newcommand{\upi}{\mathrm{i}}               % imaginary number
\newcommand{\upj}{\mathrm{j}}               % imaginary number
\newcommand{\upk}{\mathrm{k}}               % quaternion number

% Custom colors (xcolor package)
\definecolor{blue}{RGB}{76,76,255}          % 4c4cff
\definecolor{azure}{RGB}{0,128,255}         % 0080ff
\definecolor{cyan}{RGB}{0,209,209}          % 00d1d1
\definecolor{green}{RGB}{15,180,0}          % 0fb400
\definecolor{lime}{RGB}{180,216,0}          % b4d800
\definecolor{yellow}{RGB}{255,191,0}        % ffbf00
\definecolor{orange}{RGB}{255,128,0}        % ff8000
\definecolor{red}{RGB}{255,38,38}           % ff2626
\definecolor{magenta}{RGB}{255,0,255}       % ff00ff
\definecolor{purple}{RGB}{153,51,255}       % 9933ff

\definecolor{darkblue}{RGB}{25,25,127}
\definecolor{darkazure}{RGB}{25,76,127}
\definecolor{darkcyan}{RGB}{22,112,112}
\definecolor{darkgreen}{RGB}{27,103,21}
\definecolor{darklime}{RGB}{99,114,23}
\definecolor{darkyellow}{RGB}{127,102,25}
\definecolor{darkorange}{RGB}{127,76,25}
\definecolor{darkred}{RGB}{127,25,25}
\definecolor{darkmagenta}{RGB}{127,25,127}
\definecolor{darkpurple}{RGB}{76,25,127}

% Background color
\ifnum\@colortheme=0 % blue
    \colorlet{dark}{darkblue}
\fi
\ifnum\@colortheme=1 % azure
    \colorlet{dark}{darkazure}
\fi
\ifnum\@colortheme=2 % cyan
    \colorlet{dark}{darkcyan}
\fi
\ifnum\@colortheme=3 % green
    \colorlet{dark}{darkgreen}
\fi
\ifnum\@colortheme=4 % lime
    \colorlet{dark}{darklime}
\fi
\ifnum\@colortheme=5 % yellow
    \colorlet{dark}{darkyellow}
\fi
\ifnum\@colortheme=6 % orange
    \colorlet{dark}{darkorange}
\fi
\ifnum\@colortheme=7 % red
    \colorlet{dark}{darkred}
\fi
\ifnum\@colortheme=8 % magenta
    \colorlet{dark}{darkmagenta}
\fi
\ifnum\@colortheme=9 % purple
    \colorlet{dark}{darkpurple}
\fi

% Arbitrary object placement (textpos package)
\newcommand{\pos}[4][]{%
    \if\relax\detokenize{#1}\relax
        \sbox{\@tempboxa}{#4}
        \begin{textblock*}{\wd\@tempboxa}(#2,#3)
            \usebox{\@tempboxa}
        \end{textblock*}
    \else
        \begin{textblock*}{#1}(#2,#3)
            #4
        \end{textblock*}
    \fi}

% TikZ predefined drawings
\tikzset{
    antlogocompass/.pic={
        % (-0.866,-0.866) to (0.866,0.866) with 0.245 spikes
        \node[transform shape, scale=0.15, opacity=0.03, above]
            at (0.433,0) {DW~\@version};
        \foreach \phi in {0,90,180,270}{
            \fill[rotate=\phi]
                (-0.034,0) -- (-0.01,0.024) -- (-0.01,1.111) -- (0,1.111) --
                (0.054,0.864) arc (86.41:3.59:0.866) -- (0.687,0.093)
                arc (7.74:82.26:0.693) -- (0.198,0.212) -- (0.444,0.458) --
                (0.451,0.451) -- (0.342,0.169) -- (0.200,0.200) -- (0,0);}},
    antlogoant/.pic={
        % (-0.866,-0.334) to (0.866,0.334)
        \fill[even odd rule] % A
            (-0.866,-0.334) -- (-0.668,0.334) -- (-0.473,0.334) --
            (-0.275,-0.334) -- (-0.450,-0.334) -- (-0.476,-0.231) --
            (-0.665,-0.231) -- (-0.691,-0.334)
            (-0.571,0.164) -- (-0.506,-0.092) -- (-0.636,-0.092);
        \fill                % N
            (-0.202,-0.334) -- (-0.202,0.334) -- (-0.04,0.334) --
            (0.156,-0.009) -- (0.156,0.334) -- (0.318,0.334) --
            (0.318,-0.334) -- (0.162,-0.334) -- (-0.041,-0.02) --
            (-0.041,-0.334);
        \fill                % T
            (0.56,-0.334) -- (0.56,0.194) -- (0.414,0.194) --
            (0.414,0.334) -- (0.866,0.334) -- (0.866,0.194) --
            (0.72,0.194) -- (0.72,-0.334);},
    antlogocenter/.pic={
        % (-0.866,-0.177) to (0.866,0.177)
        \fill                % C
            (-0.655,0.085) -- (-0.608,0.106)
            .. controls (-0.625,0.146) and (-0.657,0.184) .. (-0.721,0.184)
            .. controls (-0.8,0.184) and (-0.866,0.126) .. (-0.866,0.001)
            .. controls (-0.866,-0.124) and (-0.8,-0.184) .. (-0.721,-0.184)
            .. controls (-0.657,-0.184) and (-0.625,-0.141) ..
            (-0.608,-0.101) -- (-0.655,-0.08)
            .. controls (-0.66,-0.094) and (-0.68,-0.132) .. (-0.721,-0.132)
            .. controls (-0.78,-0.132) and (-0.807,-0.064) .. (-0.807,0.001)
            .. controls (-0.807,0.066) and (-0.78,0.133) .. (-0.721,0.133)
            .. controls (-0.68,0.133) and (-0.66,0.099) .. (-0.655,0.085);
        \fill                % E
            (-0.534,-0.177) -- (-0.534,0.177) -- (-0.34,0.177) --
            (-0.34,0.125) -- (-0.475,0.125) -- (-0.475,0.029) --
            (-0.343,0.029) -- (-0.343,-0.023) -- (-0.475,-0.023) --
            (-0.475,-0.126) -- (-0.34,-0.126) -- (-0.34,-0.177);
        \fill                % N
            (-0.256,-0.177) -- (-0.256,0.177) -- (-0.197,0.177) --
            (-0.068,-0.055) -- (-0.068,0.177) -- (-0.009,0.177) --
            (-0.009,-0.177) -- (-0.066,-0.177) -- (-0.198,0.062) --
            (-0.198,-0.177);
        \fill                % T
            (0.146,-0.177) -- (0.146,0.125) -- (0.064,0.125) --
            (0.064,0.177) -- (0.285,0.177) -- (0.285,0.125) --
            (0.203,0.125) -- (0.203,-0.177);
        \fill                % E
            (0.36,-0.177) -- (0.36,0.177) -- (0.554,0.177) --
            (0.554,0.125) -- (0.419,0.125) -- (0.419,0.029) --
            (0.551,0.029) -- (0.551,-0.023) -- (0.419,-0.023) --
            (0.419,-0.126) -- (0.554,-0.126) -- (0.554,-0.177);
        \fill[even odd rule] % R
            (0.637,-0.177) -- (0.637,0.177) -- (0.77,0.177)
            .. controls (0.805,0.177) and (0.863,0.146) .. (0.863,0.062)
            .. controls (0.863,0.011) and (0.830,-0.024) .. (0.795,-0.034) --
            (0.866,-0.177) -- (0.799,-0.177) -- (0.739,-0.042) --
            (0.694,-0.042) -- (0.694,-0.177) -- (0.637,-0.177)
            (0.694,0.008) -- (0.694,0.125) -- (0.76,0.125)
            .. controls (0.78,0.125) and (0.803,0.106) .. (0.803,0.062)
            .. controls (0.803,0.036) and (0.78,0.008) .. (0.75,0.008);},
    afitlogo/.pic={
        % (-0.866,-0.259) to (0.866,0.388)
        % Atomic structure
        \fill (-0.541,0.063) circle (0.0146);
        \fill (-0.623,-0.015)
            .. controls (-0.81,-0.075) and (-0.860,-0.071) .. (-0.866,-0.056)
            .. controls (-0.876,-0.033) and (-0.8,0.035) .. (-0.548,0.140)
            .. controls (-0.27,0.244) and (0,0.29) .. (0.44,0.34)
            .. controls (0,0.273) and (-0.27,0.228) .. (-0.548,0.116)
            .. controls (-0.685,0.053) and (-0.786,-0.010) .. (-0.773,-0.032)
            .. controls (-0.760,-0.050) and (-0.700,-0.029) .. (-0.623,-0.015);
        \fill (-0.806,0.187)
            .. controls (-0.813,0.195) and (-0.815,0.204) .. (-0.813,0.211)
            .. controls (-0.803,0.238) and (-0.700,0.212) .. (-0.541,0.130)
            .. controls (-0.450,0.085) and (-0.246,-0.047) .. (-0.270,-0.088)
            .. controls (-0.298,-0.133) and (-0.610,0.014) .. (-0.745,0.126)
            .. controls (-0.485,-0.030) and (-0.341,-0.074) .. (-0.324,-0.051)
            .. controls (-0.302,0.010) and (-0.760,0.239) .. (-0.759,0.171);
        \fill (-0.591,0.197)
            .. controls (-0.595,0.285) and (-0.574,0.373) .. (-0.549,0.373)
            .. controls (-0.514,0.373) and (-0.484,0.220) .. (-0.482,0.061)
            .. controls (-0.478,-0.095) and (-0.498,-0.247) .. (-0.533,-0.247)
            .. controls (-0.567,-0.247) and (-0.590,-0.105) .. (-0.596,-0.014)
            -- (-0.578,-0.025)
            .. controls (-0.565,-0.120) and (-0.548,-0.182) .. (-0.528,-0.182)
            .. controls (-0.500,-0.176) and (-0.491,-0.050) .. (-0.497,0.061)
            .. controls (-0.503,0.195) and (-0.528,0.308) .. (-0.554,0.308)
            .. controls (-0.574,0.308) and (-0.582,0.250) .. (-0.591,0.197)
            (-0.599,0.035) -- (-0.601,0.100) -- (-0.587,0.110) --
            (-0.583,0.030);
        % Jet
        \fill (0.498,0.377) -- (0.527,0.369) -- (0.647,0.381) --
            (0.767,0.388) -- (0.702,0.377) -- (0.682,0.370) -- (0.658,0.350) --
            (0.651,0.346) -- (0.634,0.344) -- (0.563,0.358) -- (0.512,0.352);
        % AFIT
        \fill[even odd rule] (-0.285,-0.259) -- (-0.285,-0.243)
            .. controls (-0.265,-0.243) and (-0.260,-0.235) .. (-0.249,-0.203)
            -- (-0.111,0.293) -- (-0.105,0.293) -- (0.028,-0.203)
            .. controls (0.035,-0.230) and (0.040,-0.243) .. (0.065,-0.243)
            -- (0.065,-0.259)
            .. controls (0.035,-0.254) and (-0.023,-0.254) .. (-0.053,-0.259)
            -- (-0.053,-0.243)
            .. controls (-0.033,-0.243) and (-0.023,-0.239) .. (-0.030,-0.203)
            -- (-0.059,-0.081) -- (-0.19,-0.081) -- (-0.223,-0.203)
            .. controls (-0.231,-0.233) and (-0.215,-0.240) .. (-0.196,-0.243)
            -- (-0.196,-0.259)
            .. controls (-0.226,-0.254) and (-0.255,-0.254) .. (-0.285,-0.259)
            (-0.181,-0.05) -- (-0.122,0.168) -- (-0.067,-0.05);
        \fill (0.079,-0.259) -- (0.079,-0.243)
            .. controls (0.111,-0.240) and (0.117,-0.240) .. (0.120,-0.203)
            -- (0.120,0.235)
            .. controls (0.117,0.263) and (0.111,0.267) .. (0.079,0.269)
            -- (0.079,0.285)
            .. controls (0.109,0.280) and (0.303,0.280) .. (0.333,0.285)
            -- (0.338,0.166) -- (0.328,0.165)
            .. controls (0.320,0.225) and (0.300,0.252) .. (0.245,0.252)
            -- (0.182,0.252)
            .. controls (0.176,0.252) and (0.173,0.246) .. (0.173,0.240)
            -- (0.173,0.038) -- (0.253,0.039)
            .. controls (0.273,0.045) and (0.285,0.055) .. (0.289,0.112)
            -- (0.299,0.112) -- (0.300,-0.063) -- (0.290,-0.063)
            .. controls (0.284,-0.006) and (0.275,0.004) .. (0.253,0.010)
            -- (0.173,0.010) -- (0.173,-0.203)
            .. controls (0.175,-0.235) and (0.180,-0.242) .. (0.214,-0.243)
            -- (0.214,-0.259)
            .. controls (0.184,-0.254) and (0.109,-0.254) .. (0.079,-0.259);
        \fill (0.380,-0.259) -- (0.380,-0.243)
            .. controls (0.415,-0.241) and (0.418,-0.235) .. (0.420,-0.203)
            -- (0.420,0.235)
            .. controls (0.417,0.261) and (0.413,0.267) .. (0.380,0.269)
            -- (0.380,0.285)
            .. controls (0.410,0.280) and (0.485,0.280) .. (0.515,0.285)
            -- (0.515,0.269)
            .. controls (0.485,0.266) and (0.477,0.265) .. (0.474,0.235)
            -- (0.474,-0.203)
            .. controls (0.475,-0.233) and (0.481,-0.243) .. (0.515,-0.243)
            -- (0.515,-0.259)
            .. controls (0.485,-0.254) and (0.410,-0.254) .. (0.380,-0.259);
        \fill (0.641,-0.259) -- (0.641,-0.243)
            .. controls (0.674,-0.242) and (0.681,-0.235) .. (0.681,-0.203)
            -- (0.682,0.254) -- (0.652,0.254)
            .. controls (0.595,0.253) and (0.570,0.207) .. (0.567,0.149)
            -- (0.555,0.149)
            .. controls (0.559,0.179) and (0.562,0.255) .. (0.563,0.285)
            .. controls (0.593,0.280) and (0.828,0.280) .. (0.858,0.285)
            .. controls (0.859,0.255) and (0.862,0.179) .. (0.866,0.149)
            -- (0.854,0.149)
            .. controls (0.851,0.207) and (0.826,0.253) .. (0.769,0.254)
            -- (0.735,0.254) -- (0.735,-0.203)
            .. controls (0.738,-0.237) and (0.744,-0.242) .. (0.776,-0.243)
            -- (0.777,-0.259)
            .. controls (0.747,-0.254) and (0.671,-0.254) .. (0.641,-0.259);}
    }

\pgfplotsset{every axis/.append style={
    label style={font=\large\rmfamily},
    grid style={dark!90!white},
    tick label style={font=\large\rmfamily},
    tick style={dark!50!white},
}}

% ---------------------------
% General document formatting
% ---------------------------

% Default plot settings
\pgfplotsset{%
    tick label style={font=\scriptsize\sffamily},
    every axis label={font=\sffamily},
    legend style={font=\tiny\sffamily},
    label style={font=\sffamily},
    nodes near coords style=
        {font=\scriptsize,/pgf/number format/assume math mode}
}

% Page style
\pagestyle{empty}

% Paper height and width
\newlength{\@Height}\setlength{\@Height}{27.94cm}
\newlength{\@Width}\setlength{\@Width}{1.41421356\@Height}
\newlength{\@Margin}\setlength{\@Margin}{1cm}
% Abstract margin
\newlength{\@abMargin}\setlength{\@abMargin}{1.5cm}
\newlength{\@abModMargin}\setlength{\@abModMargin}{\@abMargin}
\addtolength{\@abModMargin}{-0.105cm}

% Define the layouts.
\ifnum\@layouttheme=0 % center
    % Abstract x and y position (top-left corner)
    \newlength{\@abX}\setlength{\@abX}{0.2875\@Width}
    \newlength{\@abY}\setlength{\@abY}{0cm}
    % Abstract width, textblock width, and height
    \newlength{\@abWidth}\setlength{\@abWidth}{0.425\@Width}
    \newlength{\@abMini}\setlength{\@abMini}{\@abWidth}
    \addtolength{\@abMini}{-0.2cm}
    \newlength{\@abHeight}\setlength{\@abHeight}{\@Height}
    % QR code height and x and y position (top-left corner)
    \newlength{\@qrHeight}\setlength{\@qrHeight}{2.5cm}
    \newlength{\@qrX}\setlength{\@qrX}{\@abX}
    \addtolength{\@qrX}{\@abMargin}
    \newlength{\@qrY}\setlength{\@qrY}{\@abHeight}
    \addtolength{\@qrY}{-\@abMargin}
    \addtolength{\@qrY}{-\@qrHeight}
    % Abstract inner height, inner width, and inner X position.
    \newlength{\@abInnerHeight}\setlength{\@abInnerHeight}{\@abHeight}
    \addtolength{\@abInnerHeight}{-3\@abMargin}
    \addtolength{\@abInnerHeight}{-\@qrHeight}
    \newlength{\@abInnerWidth}\setlength{\@abInnerWidth}{\@abWidth}
    \addtolength{\@abInnerWidth}{-2\@abMargin}
    % Author x position, width, and color
    \newlength{\@authorIX}\setlength{\@authorIX}{\@abX}
    \addtolength{\@authorIX}{\@abMargin}
    \newlength{\@authorIIX}\setlength{\@authorIIX}{\@authorIX}
    \newlength{\@authorWidth}\setlength{\@authorWidth}{\@abInnerWidth}
    \addtolength{\@authorWidth}{-2\@qrHeight}
    \addtolength{\@authorWidth}{-\@Margin} % margins from QR and logo
    \colorlet{authorColor}{white}
    % Author y positions
    \newlength{\@authorIY}\setlength{\@authorIY}{\@qrY}
    \addtolength{\@authorIY}{-0.18cm}
    \newlength{\@authorIVY}\setlength{\@authorIVY}{\@qrY}
    \addtolength{\@authorIVY}{\@qrHeight}
    \addtolength{\@authorIVY}{-0.525cm}
    \newlength{\@authorSpan}\setlength{\@authorSpan}{\@authorIVY}
    \addtolength{\@authorSpan}{-\@authorIY}
    \newlength{\@authorSep}\setlength{\@authorSep}{0.3333\@authorSpan}
    \newlength{\@authorIIY}\setlength{\@authorIIY}{\@authorIY}
    \addtolength{\@authorIIY}{\@authorSep}
    \newlength{\@authorIIIY}\setlength{\@authorIIIY}{\@authorIIY}
    \addtolength{\@authorIIIY}{\@authorSep}
    % Logo x and y positions
    \newlength{\@logoX}\setlength{\@logoX}{25.50392cm}
    \newlength{\@logoY}\setlength{\@logoY}{23.94cm}
    \newlength{\@afitX}\setlength{\@afitX}{23.1348cm}
    \newlength{\@afitY}\setlength{\@afitY}{24.58cm}
    \newlength{\@compassX}\setlength{\@compassX}{25.28cm}
    \newlength{\@compassY}\setlength{\@compassY}{23.59cm}
    \newlength{\@antAfitX}\setlength{\@antAfitX}{25.50392cm}
    \newlength{\@antAfitY}\setlength{\@antAfitY}{23.80917cm}
    % Overall paper geometry
    \newlength{\@columnSep}\setlength{\@columnSep}{0.5\@Width}
    \addtolength{\@columnSep}{2cm}
    \geometry{
        paperwidth=\@Width,
        paperheight=\@Height,
        margin=\@Margin,
        columnsep=17.5cm}
\fi
\ifnum\@layouttheme=1 % left
    % Abstract x and y position (top-left corner)
    \newlength{\@abX}\setlength{\@abX}{0cm}
    \newlength{\@abY}\setlength{\@abY}{0cm}
    % Abstract width, textblock width, and height
    \newlength{\@abWidth}\setlength{\@abWidth}{0.325\@Width}
    \newlength{\@abMini}\setlength{\@abMini}{\@abWidth}
    \addtolength{\@abMini}{-0.25cm}
    \newlength{\@abHeight}\setlength{\@abHeight}{\@Height}
    % QR code height and x and y position (top-left corner)
    \newlength{\@qrHeight}\setlength{\@qrHeight}{2.5cm}
    \newlength{\@qrX}\setlength{\@qrX}{\@abX}
    \addtolength{\@qrX}{\@abMargin}
    \newlength{\@qrY}\setlength{\@qrY}{\@abHeight}
    \addtolength{\@qrY}{-\@abMargin}
    \addtolength{\@qrY}{-\@qrHeight}
    % Abstract inner height, inner width, and inner X position.
    \newlength{\@abInnerHeight}\setlength{\@abInnerHeight}{\@abHeight}
    \addtolength{\@abInnerHeight}{-3\@abMargin}
    \addtolength{\@abInnerHeight}{-\@qrHeight}
    \newlength{\@abInnerWidth}\setlength{\@abInnerWidth}{\@abWidth}
    \addtolength{\@abInnerWidth}{-2\@abMargin}
    % Author x position, width, and color
    \newlength{\@authorIX}\setlength{\@authorIX}{\@abWidth}
    \addtolength{\@authorIX}{\@Margin}
    \newlength{\@authorIIX}\setlength{\@authorIIX}{27.3028cm}
    \newlength{\@authorWidth}\setlength{\@authorWidth}{10cm}
    \colorlet{authorColor}{black}
    % Author y positions
    \newlength{\@authorIY}\setlength{\@authorIY}{25.696536cm}
    \newlength{\@authorIIY}\setlength{\@authorIIY}{26.415023cm}
    \newlength{\@authorIIIY}\setlength{\@authorIIIY}{25.696536cm}
    \newlength{\@authorIVY}\setlength{\@authorIVY}{26.415023cm}
    % Logo x and y positions
    \newlength{\@logoX}\setlength{\@logoX}{10.96175cm}
    \newlength{\@logoY}\setlength{\@logoY}{23.94cm}
    \newlength{\@afitX}\setlength{\@afitX}{8.59263cm}
    \newlength{\@afitY}\setlength{\@afitY}{24.58cm}
    \newlength{\@compassX}\setlength{\@compassX}{10.73783cm}
    \newlength{\@compassY}\setlength{\@compassY}{23.59cm}
    \newlength{\@antAfitX}\setlength{\@antAfitX}{10.96175cm}% 25.50392cm
    \newlength{\@antAfitY}\setlength{\@antAfitY}{23.80917cm}% 23.80917cm
    % Overall paper geometry
    \newlength{\@leftMargin}\setlength{\@leftMargin}{\@abWidth}
    \addtolength{\@leftMargin}{\@Margin}
    \newlength{\@bottomMargin}\setlength{\@bottomMargin}{2\@Margin}
    \addtolength{\@bottomMargin}{0.5\@qrHeight}
    \geometry{
        paperwidth=\@Width,
        paperheight=\@Height,
        margin=\@Margin, lmargin=\@leftMargin, bmargin=\@bottomMargin}
\fi

% Paragraph parameters
\renewcommand{\sfdefault}{qhv} % TeX Gyre Heros (like Helvetica)
\setlength{\parindent}{0pt} % indentation of first line of paragraphs
\setlength{\parskip}{1.2ex plus 6pt} % space between paragraphs
\setlength{\emergencystretch}{\hsize} % overfull avoidance
\overfullrule=1mm % overfull black bands
\raggedbottom % page-ending empty space allowed

% Reduce list spacing (enumitem package).
\setlist{
    itemindent=0pt,
    labelsep=0.2cm,
    labelwidth=0.4cm,
    leftmargin=0.6cm,
    topsep=0pt,
    parsep=0.5\parskip,
    itemsep=0.5\parskip,
    before=\raggedright}

% Reduce spacing of other environments (enumitem package).
\setlist[trivlist]{topsep=0pt}

% Hyperlink styling
\hypersetup{
    colorlinks=false, hidelinks,
    hypertexnames=false}

% Remove pdf image warnings.
\pdfsuppresswarningpagegroup=1

% Use French spacing.
\frenchspacing

% headings
\renewcommand{\section}[1]{\vspace{0.5cm}%
    {\sffamily\large\bfseries{#1}}%
    \vspace{0.125cm}}
\renewcommand{\subsection}[1]{\vspace{0.25cm}%
    {\sffamily\large\itshape{#1}}%
    \vspace{0.125cm}}
\renewcommand{\subsubsection}[1]{\vspace{0.25cm}%
    {\sffamily\bfseries{#1}}%
    \vspace{0.125cm}}

% Preamble metadata
\def\@title{}\renewcommand{\title}[1]{\def\@title{#1}}
\def\@subtitle{}\newcommand{\subtitle}[1]{\def\@subtitle{#1}}
\def\@author{}\renewcommand{\author}[1]{\def\@author{#1}}
\newif\ifsecondauthor\secondauthorfalse
\def\@secondauthor{}
\newcommand{\secondauthor}[1]{\secondauthortrue\def\@secondauthor{#1}}
\newif\ifthirdauthor\thirdauthorfalse
\def\@thirdauthor{}
\newcommand{\thirdauthor}[1]{\thirdauthortrue\def\@thirdauthor{#1}}
\newif\iffourthauthor\fourthauthorfalse
\def\@fourthauthor{}
\newcommand{\fourthauthor}[1]{\fourthauthortrue\def\@fourthauthor{#1}}
\newif\ifuseqr\useqrfalse
\def\@qrdata{}\newcommand{\qrdata}[1]{\useqrtrue\def\@qrdata{#1}}
\def\@abstract{}\renewcommand{\abstract}[1]{\def\@abstract{#1}}
\newif\ifcustomlogo\customlogofalse
\def\@customlogo{}\newcommand{\logo}[1]{\customlogotrue\def\@customlogo{#1}}

% maketitle
\newlength{\authorpos}
\renewcommand{\maketitle}{%
    % Abstract
    \begin{textblock*}{\@abWidth}(\@abX,\@abY)% {w}(x,y)
        \colorbox{dark}{% set the color
        \begin{minipage}[t][\@abHeight][t]{\@abMini}% used to create margins
            \vspace{\@abModMargin}\hspace{\@abModMargin}% 1.5 cm margins
            \begin{minipage}[t][\@abInnerHeight][t]{\@abInnerWidth}% content
                \color{white}\sffamily\raggedright%
                \sffamily\fontsize{38pt}{57pt}\selectfont%
                \@abstract%
            \end{minipage}
        \end{minipage}}
    \end{textblock*}%
    % QR code
    \ifuseqr
        \begin{textblock*}{\@qrHeight}(\@qrX,\@qrY)
            \textcolor{dark!20!white}{\qrcode[height=\@qrHeight]{\@qrdata}}
        \end{textblock*}
        \ifnum\@layouttheme=0 % center
            \addtolength{\@authorIX}{\@qrHeight}
            \addtolength{\@authorIX}{0.5cm} % margin from QR code
            \addtolength{\@authorIIX}{\@qrHeight}
            \addtolength{\@authorIIX}{0.5cm} % margin from QR code
        \fi
    \fi
    % Authors
    \ifsecondauthor
        \begin{textblock*}{\@authorWidth}(\@authorIX,\@authorIY)
            {\sffamily\color{authorColor}%
            \parbox[b]{0.8cm}{\raggedleft\tiny FIRST\\AUTHOR}
            \,\fontsize{13.2pt}{15pt}\selectfont%
            \MakeUppercase{\@author}}
        \end{textblock*}
        \begin{textblock*}{\@authorWidth}(\@authorIX,\@authorIIY)
            {\sffamily\color{authorColor}%
            \parbox[b]{0.8cm}{\raggedleft\tiny SECOND\\AUTHOR}
            \,\fontsize{13.2pt}{15pt}\selectfont%
            \MakeUppercase{\@secondauthor}}
        \end{textblock*}
        \ifthirdauthor
            \begin{textblock*}{\@authorWidth}(\@authorIIX,\@authorIIIY)
                {\sffamily\color{authorColor}%
                \parbox[b]{0.8cm}{\raggedleft\tiny THIRD\\AUTHOR}
                \,\fontsize{13.2pt}{15pt}\selectfont%
                \MakeUppercase{\@thirdauthor}}
            \end{textblock*}
            \iffourthauthor
                \begin{textblock*}{\@authorWidth}(\@authorIIX,\@authorIVY)
                    {\sffamily\color{authorColor}%
                    \parbox[b]{0.8cm}{\raggedleft\tiny FOURTH\\AUTHOR}
                    \,\fontsize{13.2pt}{15pt}\selectfont%
                    \MakeUppercase{\@fourthauthor}}
                \end{textblock*}
            \fi
        \fi
    \else
        \begin{textblock*}{\@authorWidth}(\@authorIX,\@authorIY)
            {\sffamily\color{authorColor}%
            \parbox[b]{0.8cm}{\raggedleft\tiny AUTHOR\\\,}
            \,\fontsize{13.2pt}{15pt}\selectfont%
            \MakeUppercase{\@author}}
        \end{textblock*}
    \fi
    % Logo
    \ifcustomlogo
        \begin{textblock*}{2.5cm}(\@logoX,\@logoY)
            \hfill\includegraphics[width=2.5cm, height=2.5cm, keepaspectratio]
                {\@customlogo}
        \end{textblock*}
    \else
        \ifafit
            \begin{textblock*}{5cm}(\@afitX,\@afitY)
                \tikz{\path (0,0) pic[color=dark!20!white, scale=2.87]
                    {afitlogo};}
            \end{textblock*}
        \else
            \ifcompassonly
                \begin{textblock*}{3.3cm}(\@compassX,\@compassY)
                    \tikz{\path (0,0) pic[color=dark!20!white, scale=1.4434]
                        {antlogocompass};}
                \end{textblock*}
            \else
                \begin{textblock*}{2.7cm}(\@antAfitX,\@antAfitY)
                    \tikz{%
                        %\fill[red, opacity=0.2] (-1.25,-1.25)
                        %   rectangle (1.25,1.25);
                        \def\antscale{1.78}
                        \begin{scope}[scale=\antscale, yshift=0.305cm]
                            % (-0.702,-0.26) to (0.702,0.26)
                            \path (-0.442,0.137) pic[scale={0.3*\antscale},
                                color=dark!20!white] {antlogocompass};
                            \path (0.347,0.26) pic[scale={0.41*\antscale},
                                color=dark!20!white] {antlogoant};
                            \path (0.347,-0.05) pic[scale={0.41*\antscale},
                                color=dark!20!white] {antlogocenter};
                        \end{scope}
                        % (-0.866,-0.259) to (0.866,0.388)
                        \path (0,-0.8762) pic[scale=1.4434,
                            color=dark!20!white] {afitlogo};
                    }
                \end{textblock*}
            \fi
        \fi
    \fi
    % FIXME Remove spacing blocks
    %\ifnum\@layouttheme=0 % center
    %    \begin{textblock*}{\@abMargin}(9.8782817cm,0cm)
    %        \tikz{\fill[white, opacity=0.1]
    %            (-0.75,0) rectangle (0.75,\@abHeight);}
    %    \end{textblock*}
    %    \begin{textblock*}{\@abMargin}(28.134845cm,0cm)
    %        \tikz{\fill[white, opacity=0.1]
    %            (-0.75,0) rectangle (0.75,\@abHeight);}
    %    \end{textblock*}
    %\fi
    %\ifnum\@layouttheme=1 % left
    %    \begin{textblock*}{\@abMargin}(0cm,0cm)
    %        \tikz{\fill[white, opacity=0.1]
    %            (-0.75,0) rectangle (0.75,\@abHeight);}
    %    \end{textblock*}
    %    \begin{textblock*}{\@abMargin}(13.59267cm,0cm)
    %        \tikz{\fill[white, opacity=0.1]
    %            (-0.75,0) rectangle (0.75,\@abHeight);}
    %    \end{textblock*}
    %\fi
    % Title and subtitle
    {%
        \sffamily\raggedright%
        \itshape\LARGE\textbf{\@title}\\[0.5cm]
        \itshape\Large\@subtitle\vspace{0.25cm}}%
    }
