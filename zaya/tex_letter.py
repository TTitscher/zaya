import argparse
from pathlib import Path

r"""
\documentclass[]{scrlttr2}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[ngerman]{babel}
\usepackage{geometry}
\geometry{verbose,a4paper,tmargin=25mm,bmargin=25mm,lmargin=20mm,rmargin=20mm}
\usepackage{graphicx}

\setlength{\parindent}{0pt}

\begin{document}
	
\pagestyle{empty}

\setkomavar{fromname}{...}
\setkomavar{fromaddress}{...}
\setkomavar{signature}{...}

\setkomavar{subject}{...}


\begin{letter}{ 
	% address of recipient
...\\
...\\
...
}

\opening{Sehr geehrte Damen und Herren,}

body 1

\phantom{m}

body 2

\closing{Mit freundlichen Grüßen}

\end{letter}
\end{document}
"""

def cli():
    """
    Create and parse the command line interface.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--to", help=r"Adress of the recipient. Lines are seperated by \\", type=str)
    parser.add_argument("--body", help=r"Body of the letter. New paragraph by \\", type=str)
    parser.add_argument("--sender", help=r"Text file containing the sender adress", type=Path, default = Path(__file__).parent/"tex_letter_sender.txt")
    parser.add_argument("--date", help=r"Text file containing the sender adress", type=str, default=r"\today")
    parser.add_argument("--subject", help=r"Subject of the letter", type=str)
    parser.add_argument("--closing", help=r"Closing", type=str, default = "Mit freundlichen Grüßen")
    parser.add_argument("--tex", help="Path to .tex doc. If none, STDOUT")
    return parser.parse_args()

def to_tex_string(args):
    """
    Transform the CLI ``args`` to a LaTeX file using KOMA "scrlttr2".
    """

    # define some KOMA terms

    with open(args.sender, "r") as f:
        sender_lines = f.readlines()
    fromname = sender_lines[0]
    fromaddress = "".join(sender_lines[1:])

    recipient = "\n".join(args.to.split(r"\n"))

    print(f"{sender} \nto \n{recipient} \n  {args.subject} {args.date} \n{args.body}\n\n{args.closing}\n{name}")

def main():
    args = cli()

if __name__ == "__main__":
    main()
