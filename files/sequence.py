"""
Sequence Class.

This module allows the user to store one 
sequence from a fasta file in a Sequence instance.

Classes
-------
Sequence
"""

class Sequence:
    """A class to represent a sequence."""

    def __init__(self, seq: str) -> None:
        """Construct all attributes for Sequence."""
        self.seqStr = seq
        self.seqLst: list[str] = list(seq)

    @property
    def seqStr(self) -> str:
        """Sequence as string."""
        return self._seqStr
    
    @seqStr.setter
    def seqStr(self, seq: str) -> None:
        if isinstance(seq, str):
            self._seqStr = seq
        else:
            raise ValueError('"seq" must be a str')
        
    @property
    def seqLst(self) -> list[str]:
        """Sequence as list."""
        return self._seqLst
    
    @seqLst.setter
    def seqLst(self, seq: list[str]) -> None:
        self._seqLst = seq
        
    def toList(self) -> None:
        """Convert Sequence to list."""
        self.seqLst = list(self.seqStr)

    def reverse(self) -> None:
        """Reverse Sequence as list."""
        self.seqLst.reverse()

    def getBase(self, pos: int) -> str:
        """Return base pair in Sequence."""
        length: int = self.getLength()
        if (pos >= length):
            raise ValueError(f"{pos} greater than length {length}")
        base: str = self.seqStr[pos]
        return base
    
    def getLength(self) -> int:
        """Return Sequence length."""
        length: int = len(self.seqStr)
        return length
    
    def print(self, rep: str) -> None:
        """Print Sequence as string or list."""
        if rep == "string":
            print(self.seqStr)
        if rep == "list":
            print("".join(self.seqLst))
