import shutil
import subprocess
import sys
from pathlib import Path
from string import Template
from tempfile import TemporaryDirectory

import yaml
from loguru import logger

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


def main():
    try:
        inp = Path(sys.argv[1])
    except IndexError as e:
        raise RuntimeError("Specify a YAML file.") from e

    with open(inp, "r") as f:
        letter_dict = yaml.load(f, Loader=yaml.Loader)

    letter_dict.setdefault("date", "\\today")
    letter_dict.setdefault("closing", "Mit freundlichen Grüßen")
    letter_dict.setdefault("opening", "Sehr geehrte Damen und Herren,")
    letter_dict.setdefault("signature", letter_dict["fromname"])

    for key, value in letter_dict.items():
        num_newlines = value.count("\n")
        letter_dict[key] = value.replace("\n", "\\\\ \n", num_newlines - 1)

    with TemporaryDirectory() as d:
        outdir = Path(d).absolute()
        tex = outdir / inp.with_suffix(".tex")
        with open(tex, "w") as f:
            logger.info(f"Created {tex}")
            f.write(template.substitute(letter_dict))

        latex_cmd = ["latexmk", "-pdf", f"-outdir={outdir}", str(tex)]
        p = subprocess.run(latex_cmd, capture_output=True)
        logger.debug(p.stderr.decode("utf-8"))

        src = tex.with_suffix(".pdf")
        dst = Path.cwd() / tex.with_suffix(".pdf").name
        shutil.copy2(src, dst)
        logger.info(f"Created {dst}")


if __name__ == "__main__":
    main()
