# Jupytext

Jupyter notebooks generate files that may contain metadata, source code,
formatted text, and rich media. Unfortunately, this makes these files poor
candidates for conventional version control solutions, which works best with
plain text.

[Jupytext](https://github.com/mwouts/jupytext) creates tracable python scripts
from Jupyter notebooks.

## Installation

```bash
conda install -c conda-forge jupytext
```

## Use Jupytext

Jupytext can save Jupyter notebooks as

    - Markdown and R Markdown documents,
    - Scripts in many languages.

It can also convert these documents into Jupyter Notebooks, allowing you to
synchronize content in both directions.

The languages that are currently supported by Jupytext are: Julia, Python, R,
Bash, Scheme, Clojure, Matlab, Octave, C++, q/kdb+, IDL, TypeScript, Javascript
and Scala. In addition, jupytext users can choose between two formats for
notebooks as scripts:

- The `percent` format, compatible with several IDEs, including Spyder,
  Hydrogen, VScode and PyCharm. In that format, cells are delimited with a
  commented %%.
- The `light` format, designed for this project. Use that format to open
  standard scripts as notebooks, or to save notebooks as scripts with few cell
  markers - none when possible.

### Python script

```bash
jupytext --to py simple-nb.ipynb
```
