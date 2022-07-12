import shutil
import subprocess
import sys
import argparse
from pathlib import Path
from string import Template
from tempfile import TemporaryDirectory


import yaml

template = Template(
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

\setkomavar{fromname}{$fromname}
\setkomavar{fromaddress}{$fromadress}
\setkomavar{signature}{$signature}

\setkomavar{subject}{$subject}


\begin{letter}{ 
$to
}

\opening{$opening}

$body

\closing{$closing}

\end{letter}
\end{document}
"""
)


def fix_defaults(letter_dict):
    letter_dict.setdefault("date", "\\today")
    letter_dict.setdefault("closing", "Mit freundlichen Grüßen")
    letter_dict.setdefault("opening", "Sehr geehrte Damen und Herren,")
    letter_dict.setdefault("signature", letter_dict["fromname"])
    letter_dict.setdefault("pdf", True)

    for key, value in letter_dict.items():
        if isinstance(value, str):
            num_newlines = value.count("\n")
            letter_dict[key] = value.replace("\n", "\\\\ \n", num_newlines - 1)


def compile_pdf(tex):
    with TemporaryDirectory() as d:
        tex_dir = Path(d).absolute()

        latex_cmd = ["latexmk", "-pdf", f"-outdir={tex_dir}", str(tex)]
        p = subprocess.run(latex_cmd, capture_output=True)

        src = tex_dir / tex.with_suffix(".pdf").name
        dst = tex.with_suffix(".pdf")
        shutil.copy2(src, dst)
        print(f"Created {dst}")


def main():
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):         
        def _fill_text(self, text, width, indent):             
            return "".join(indent + line for line in text.splitlines(keepends=True))     

    desc = r"""
    TEXIFY

    Converts an `input.yml` to a LaTeX DIN letter.

    Required YML keys
    =================
    - fromname:     name of sender
    - fromadress:   adress of sender
    - to:           name and adress of recipient
    - subject
    - body
    
    Optional YML keys
    =================
    - date:         defaults to \today
    - opening:      defaults to "Sehr geehrte Damen und Herren,"
    - closing:      defaults to "Mit freundlichen Grüßen"
    - signature:    defaults to `fromname`
    - pdf:          compile the pdf? Default: True
    """

    parser = argparse.ArgumentParser(description=desc, formatter_class=MyFormatter)
    parser.add_argument("input", help="Input YAML file.", type=Path)

    args = parser.parse_args()

    with open(args.input, "r") as f:
        letter_dict = yaml.load(f, Loader=yaml.Loader)

    fix_defaults(letter_dict)

    # write tex file in CWD in any case
    tex = Path.cwd() / args.input.with_suffix(".tex")
    with open(tex, "w") as f:
        f.write(template.substitute(letter_dict))
        print(f"Created {tex}")

    if letter_dict["pdf"]:
        compile_pdf(tex)


if __name__ == "__main__":
    main()
