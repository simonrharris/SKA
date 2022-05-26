"""ska: split k-mer alignment

Usage:
  ska fasta <files>... -o <output>
  ska align -l <file-list> -o <output>
  ska map ref -l <file-list> -o <output>
  ska (-h | --help)
  ska (--version)

Options:
  -h --help     Show this help.
  --version     Show version.

  -o <output>    Output prefix.
  -l <file-list> File with a list of input files.
"""

import os, sys
import re

from .__init__ import __version__

import ska

def get_options():
    from docopt import docopt
    arguments = docopt(__doc__, version="ska v"+__version__)

    # TODO Check options here

    return arguments

def main():
    args = get_options()

    #
    # Create a database (sketch input)
    #
    if args["fasta"]:
        print("run ska fasta")
    elif args["align"]:
        print("run ska align")
        ska.align(args["file-list"])
    elif args["map"]:
        print("run ska map")
    else:
        print("Option error!")
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()