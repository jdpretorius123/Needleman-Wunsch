"""
Write Needleman-Wunsch Output.

This module allows the user to setup and run the 
Needleman-Wunsch algorithm.

Functions
---------
writeAlignment(argv: list[str]) -> None:
    Write alignment results.
"""

import os
from sequence import Sequence
from file import MatrixFile, FastaFile
from nw import Linear, Affine

SUB_MATRIX = dict[tuple[str, str], int]

def writeAlignment(argv: list[str]) -> None:
    """Write alignment results."""
    fasta1: FastaFile = FastaFile(argv[1])
    fasta2: FastaFile = FastaFile(argv[2])

    mf: MatrixFile = MatrixFile(argv[3])
    submatrix: SUB_MATRIX = mf.generate()

    seqs1: dict[int, Sequence] = fasta1.generate()
    seqs2: dict[int, Sequence] = fasta2.generate()

    outfile: str = argv[4]
    if os.path.isfile(outfile):
            os.remove(outfile)

    gap: int = int(argv[5])
    score: int = int(argv[6])

    if not score:
        for i in range(len(seqs1)):
            l: Linear = Linear(seqs1[i], seqs2[i], submatrix, gap)
            l.execute(i + 1, 0, outfile)
    else:
        for j in range(len(seqs1)):
            a: Affine = Affine(seqs1[j], seqs2[j], submatrix, gap)
            a.extend = float(argv[7])
            a.execute(j + 1, 0, outfile)
    