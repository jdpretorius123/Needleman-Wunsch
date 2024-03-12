"""
Run Needleman-Wunsch Algorithm.

This module allows the user to run the 
Needleman-Wunsch algorithm.
"""

import sys
import process

argv: list[str] = sys.argv
if (len(argv) != 8):
    print(
        """
        Usage: 
        main.py <infile1> <infile2> <matrixfile> <outfile> <gap> <score> <extend>
        """
        )
    sys.exit("Please enter the correct input.")

process.writeAlignment(argv)