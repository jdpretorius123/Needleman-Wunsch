"""
Matrix Class.

This module allows the user to create and fill
a matrix (score or traceback).

Classes
-------
Matrix
"""

MATRIX = list[list[float | str]]

class Matrix:
    """A class to represent a matrix."""

    def __init__(self, nrows: int, ncols: int) -> None:
        """Construct all attributes for Matrix."""
        self.nrows = nrows
        self.ncols = ncols
        self.matrix: MATRIX = [[0]]

    @property
    def nrows(self) -> int:
        """Number of rows of Matrix."""
        return self._nrows

    @nrows.setter
    def nrows(self, nrows: int) -> None:
        if isinstance(nrows, int):
            self._nrows = nrows
        else:
            raise ValueError('"nrows" must be an int')

    @property
    def ncols(self) -> int:
        """Number of columns of Matrix."""
        return self._ncols

    @ncols.setter
    def ncols(self, ncols: int) -> None:
        if isinstance(ncols, int):
            self._ncols = ncols
        else:
            raise ValueError('"ncols" must be an int')

    @property
    def matrix(self) -> MATRIX:
        """Number of rows of Matrix."""
        return self._matrix
    
    @matrix.setter
    def matrix(self, matrix: MATRIX) -> None:
        self._matrix = matrix

    def setValue(self, value: float | str, row: int, col: int) -> None:
        """Set value of row,col in Matrix."""
        self.matrix[row][col] = value

    def getValue(self, row: int, col: int) -> float | str:
        """Get value of row,col in Matrix."""
        if (row > self.nrows or col > self.ncols):
            raise ValueError(
                f"{row} or {col} is out of range: nrows={self.nrows}, ncols={self.ncols}"
                )
        value: float | str = self.matrix[row][col]
        return value

    def initialize(self, valueType: str) -> None:
        """Initialize zero Matrix."""
        matrix: MATRIX = list()
        for i in range(self.nrows):
            matrix.append([])
            for j in range(self.ncols):
                if (valueType == "string"):
                    matrix[i].append("0")
                else:
                    matrix[i].append(0)
        self.matrix = matrix

    def print(self) -> None:
        """Print Matrix."""
        for i in range(self.nrows):
            end: str = ", "
            lastIndex: int = self.ncols - 1
            for j in range(self.ncols):
                if j == lastIndex:
                    end = "\n"
                print(self.matrix[i][j], sep = "\t", end = end)
                