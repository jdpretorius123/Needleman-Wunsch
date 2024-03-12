"""
File Class.

This module allows the user to process an input MatrixFile
or FastaFile and use the file's contents for the 
Needleman-Wunsch algorithm.

Classes
-------
File
MatrixFile
FastaFile
"""

SUB_MATRIX = dict[tuple[str, str], int]
KEYS = list[tuple[str, str]]

import regex
from typing import TextIO
from sequence import Sequence

class File:
    """A class to represent a file."""

    def __init__(self, path: str) -> None:
        """Construct all attributes for File."""
        self.path = path

    @property
    def path(self) -> str:
        """Path to file."""
        return self._path
    
    @path.setter
    def path(self, path: str) -> None:
        if isinstance(path, str):
            self._path = path
        else:
            raise ValueError('"path" must be a str')
        
    def _processFile(self, lines: list[str]) -> list[str]:
        print("Error: _processFile not defined for parent class File.")
        return [""]

    def generate(self) -> None:
        print("Error: generate not defined for parent class File.")

    def print(self) -> None:
        print("Error: print not defined for parent class File.")
    
    def remove(self, alist: list[str], element: str) -> list[str]:
        """Return alist without element."""
        count: int = alist.count(element)
        for i in range(count):
            alist.remove(element)
        return alist
    
class MatrixFile(File):
    """A class to represent a matrix file."""

    def __init__(self, path: str) -> None:
        """Construct all attributes for MatrixFile."""
        super().__init__(path)

    def _processFile(self, lines: list[str]) -> list[str]:
        """Return processed lines of substitution matrix file."""
        lines = self.remove(lines, "\n")
        lines = self.remove(lines, "\n")
        for i in range(len(lines)):
            lines[i] = regex.sub(r'[\n]', '', lines[i]) # type: ignore
            lines[i] = lines[i].split()
        return lines

    def _createKeys(self, bases: list[str]) -> KEYS:
        """Create keys for substitution matrix."""
        keys: KEYS = list()
        for i in range(len(bases)):
            for j in range(len(bases)):
                keys.append((bases[i], bases[j]))
        return keys

    def _fillSubMatrix(self, keys: KEYS, lines: list[str]) -> SUB_MATRIX:
        """Fill substitution matrix dictionary."""
        idx: int = 0
        submatrix: SUB_MATRIX = dict()
        for i in range(1, len(lines)):
            for j in range(1, len(lines[1])):
                submatrix[keys[idx]] = int(lines[i][j])
                idx += 1
        return submatrix

    def _createSubMatrix(self) -> SUB_MATRIX:
        """Create substitution matrix."""
        path: str = self.path
        file = open(path, 'r')
        lines: list[str] = file.readlines()
        file.close()

        lines = self._processFile(lines)
        keys: KEYS = self._createKeys(lines[0])
        submatrix: SUB_MATRIX = self._fillSubMatrix(keys, lines)
        return submatrix
    
    def generate(self) -> SUB_MATRIX:
        """Return substitution matrix."""
        submatrix: SUB_MATRIX = self._createSubMatrix()
        return submatrix
    
    def print(self) -> None:
        submatrix: SUB_MATRIX = self._createSubMatrix()
        keys: KEYS = list(submatrix.keys())
        for key in keys:
            print("Key: ", key)
            print("Value: ", submatrix[key])
    
class FastaFile(File):
    """A class to represent a fasta file."""

    def __init__(self, path: str) -> None:
        """Construct all attributes for FastaFile."""
        super().__init__(path)

    def _processFile(self, lines: list[str]) -> list[str]:
        """Process fasta file."""
        seqs: list[str] = list()
        seq: str = ""
        for line in lines:
            if not line.startswith(">"):
                line = regex.sub(r'[\n]', '', line)
                seq = seq + line
            else:
                seqs.append(seq)
                seq = ""
        seqs.append(seq)
        seqs = self.remove(seqs, "")
        return seqs

    def _fillSequences(self, lines: list[str]) -> list[str]:
        """Fill Sequence dictionary."""
        seqs: dict[int, Sequence] = dict()
        for idx in range(len(lines)):
            seqs[idx] = Sequence(lines[idx])
        return seqs

    def _createSequences(self) -> dict[int, Sequence]:
        """Create dictionary of Sequences."""
        path: str = self.path
        file: TextIO = open(path, 'r')
        lines: list[str] = file.readlines()
        file.close()

        lines = self._processFile(lines)
        seqs: dict[int, Sequence] = self._fillSequences(lines)
        return seqs
    
    def generate(self) -> dict[int, Sequence]:
        """Return sequences."""
        seqs: dict[int, Sequence] = self._createSequences()
        return seqs
    
    def print(self) -> None:
        seqs: dict[int, Sequence]= self._createSequences()
        keys: list[int] = list(seqs.keys())
        for key in keys:
            print("Key: ", key)
            print("Value: ")
            seqs[key].print("string")