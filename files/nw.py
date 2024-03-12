"""
NW Class.

This module allows the user to execute the Needleman-Wunsch
algorithm with either Linear or Affine scoring.

Classes
-------
NW
Linear
Affine
"""

from typing import TextIO
from matrix import Matrix
from sequence import Sequence

SUB_MATRIX = dict[tuple[str, str], int]


class NW:
    """A class to represent the Needleman-Wunch algorithm."""

    def __init__(
        self, seq1: Sequence, seq2: Sequence, submatrix: SUB_MATRIX, gap: float
    ) -> None:
        """Construct all attributes for NW."""
        self.seq1 = seq1
        self.seq2 = seq2
        self.submatrix = submatrix
        self.gap = float(gap)
        self.extend = -0.1

    @property
    def seq1(self) -> Sequence:
        """Sequence 1 used in algorithm."""
        return self._seq1

    @seq1.setter
    def seq1(self, seq1: Sequence) -> None:
        self._seq1 = seq1

    @property
    def seq2(self) -> Sequence:
        """Sequence 2 used in algorithm."""
        return self._seq2

    @seq2.setter
    def seq2(self, seq2: Sequence) -> None:
        self._seq2 = seq2

    @property
    def submatrix(self) -> SUB_MATRIX:
        """Substition matrix used in algorithm."""
        return self._submatrix

    @submatrix.setter
    def submatrix(self, submatrix: SUB_MATRIX) -> None:
        self._submatrix = submatrix

    @property
    def gap(self) -> float:
        """Gap penalty used in algorithm."""
        return self._gap

    @gap.setter
    def gap(self, gap: float) -> None:
        self._gap = gap

    @property
    def extend(self) -> float:
        """Gap extension penalty used in affine scoring."""
        return self._extend

    @extend.setter
    def extend(self, extend: float) -> None:
        self._extend = extend

    def _createMatrix(
        self, nrows: int, ncols: int, valueType: str, matType: str
    ) -> Matrix:
        message: str = "_createMatrix not defined for parent class NW."
        raise NotImplementedError(message)

    def _reverseSeqs(self, seq1: str, seq2: str) -> tuple[str, str]:
        """Reverse annotated sequences."""
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]
        alignment: tuple[str, str] = seq1, seq2
        return alignment

    def _annotate(self, alignment: tuple[str, str]) -> str:
        """Annotate optimal alignment."""
        annotation: str = ""
        seq1: str = alignment[0]
        seq2: str = alignment[1]
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                annotation += "|"
            if seq1[i] != seq2[i]:
                if seq1[i] == "-" or seq2[i] == "-":
                    annotation += " "
                else:
                    annotation += "*"
        return annotation

    def _scoreAlignment(
        self, alignment: tuple[str, str], annotation: str
    ) -> float:
        message: str = "_scoreAlignment not defined for parent class NW."
        raise NotImplementedError(message)

    def _calcIndels(self, annotation: str) -> int:
        """Calculate the number of indels."""
        count: int = 0
        i: int = 0
        stop: int = len(annotation)
        while i < stop:
            char: str = annotation[i]
            if char == " ":
                while char == " " and i < stop:
                    char = annotation[i]
                    i += 1
                count += 1
            else:
                i += 1
        return count

    def _calcAvgIndel(self, annotation: str) -> float:
        """Calculate average indel length."""
        count: int = 0
        counts: list[int] = list()
        for i in range(len(annotation)):
            if annotation[i] == " ":
                count += 1
            else:
                if count != 0:
                    counts.append(count)
                count = 0
        avgIndel: float = sum(counts) / len(counts)
        return avgIndel

    def _calcStats(
        self, alignment: tuple[str, str], annotation: str
    ) -> list[float]:
        """Calculate alignment statistics."""
        score: float = self._scoreAlignment(alignment, annotation)
        matches: int = annotation.count("|")
        indels: int = self._calcIndels(annotation)
        avgIndel: float = self._calcAvgIndel(annotation)
        avgLength: float = (self.seq1.getLength() + self.seq2.getLength()) / 2
        percentId: int = round((matches / avgLength) * 100)
        stats: list[float] = [
            matches,
            percentId,
            indels,
            round(avgIndel, ndigits=1),
            len(annotation),
            round(score, ndigits=1),
        ]
        return stats

    def execute(self, num: int, printOutput: int, path: str) -> None:
        print("Error: execute note defined for parent class NW.")

    def _print(
        self,
        num: int,
        stats: list[float],
        alignment: tuple[str, str],
        annotation: str,
    ) -> None:
        """Print optimal alignment."""
        seq1: str = alignment[0]
        seq2: str = alignment[1]
        text: list[str] = [
            f"Alignment #{num}:\n",
            f"Sequence #1: seq{num}A",
            f"Sequence #2: seq{num}B",
            f"Matches: {stats[0]}",
            f"Percent identity: {stats[1]}%",
            f"Indels: number={stats[2]} mean length={stats[3]}",
            f"Alignment length: {stats[4]}",
            f"Score={stats[5]}\n",
        ]
        for i in text:
            print(i)

        limit: int = 60
        for j in range(0, len(seq1), limit):
            stop: int = j + limit
            print(
                seq1[j:stop],
                annotation[j:stop],
                seq2[j:stop],
                sep="\n",
            )

    def _write(
        self,
        num: int,
        stats: list[float],
        alignment: tuple[str, str],
        annotation: str,
        path: str,
    ) -> None:
        """Write optimal alignment."""
        file: TextIO = open(path, "a")
        seq1: str = alignment[0]
        seq2: str = alignment[1]
        text: list[str] = [
            f"Alignment #{num}:\n",
            f"Sequence #1: seq{num}A",
            f"Sequence #2: seq{num}B",
            f"Matches: {stats[0]}",
            f"Percent identity: {stats[1]}%",
            f"Indels: number={stats[2]} mean length={stats[3]}",
            f"Alignment length: {stats[4]}",
            f"Score={stats[5]}\n",
        ]
        for i in text:
            file.write(i + "\n")

        limit: int = 60
        for j in range(0, len(seq1), limit):
            stop: int = j + limit
            line: str = (
                f"{seq1[j:stop]}\n{annotation[j:stop]}\n{seq2[j:stop]}\n\n"
            )
            file.write(line)
        file.close()


class Linear(NW):
    """A class to represent linear scoring global alignment."""

    def __init__(
        self, seq1: Sequence, seq2: Sequence, submatrix: SUB_MATRIX, gap: float
    ) -> None:
        """Construct all attributes for Linear."""
        super().__init__(seq1, seq2, submatrix, gap)

    def _initScore(self, score: Matrix) -> Matrix:
        """Return initialized linear score matrix."""
        for i in range(1, score.ncols):  # first row
            score.setValue(i * self.gap, 0, i)
        for j in range(1, score.nrows):  # first column
            score.setValue(j * self.gap, j, 0)
        return score

    def _initTrace(self, traceback: Matrix) -> Matrix:
        """Return initialized linear traceback matrix."""
        traceback.setValue("STOP", 0, 0)
        for i in range(1, traceback.ncols):
            traceback.setValue("LEFT", 0, i)
        for j in range(1, traceback.nrows):
            traceback.setValue("UP", j, 0)
        return traceback

    def _createMatrix(
        self, nrows: int, ncols: int, valueType: str, matType: str
    ) -> Matrix:
        """Return initialized score matrix."""
        matrix: Matrix = Matrix(nrows, ncols)
        matrix.initialize(valueType)
        match matType:
            case "score":
                matrix = self._initScore(matrix)
            case "traceback":
                matrix = self._initTrace(matrix)
        return matrix

    def _traceValue(self, scores: list[float], maxScore: float) -> str:
        """Return value for traceback matrix."""
        if maxScore == scores[0]:
            return "DIAGONAL"
        if maxScore == scores[1]:
            return "LEFT"
        if maxScore == scores[2]:
            return "UP"
        return "ERROR"

    def _fillMatrices(self, score: Matrix, traceback: Matrix) -> list[Matrix]:
        """Fill score and traceback matrices."""
        for i in range(1, score.nrows):
            for j in range(1, score.ncols):
                base1: str = self.seq1.getBase(i - 1)
                base2: str = self.seq2.getBase(j - 1)
                key: tuple[str, str] = base1, base2
                scores: list[float] = [
                    score.getValue(i - 1, j - 1) + self.submatrix[key],  # type: ignore
                    score.getValue(i - 1, j) + self.gap,  # type: ignore
                    score.getValue(i, j - 1) + self.gap,  # type: ignore
                ]
                maxScore: float = max(scores)
                score.setValue(maxScore, i, j)
                traceValue: str = self._traceValue(scores, maxScore)
                traceback.setValue(traceValue, i, j)
        matrices: list[Matrix] = [score, traceback]
        return matrices

    def _reverseSeqs(self, seq1: str, seq2: str) -> tuple[str, str]:
        """Reverse annotated sequences."""
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]
        alignment: tuple[str, str] = seq1, seq2
        return alignment

    def _getTraceback(self, traceback: Matrix) -> tuple[str, str]:
        """Get traceback for optimal alignment."""
        seq1: str = ""
        seq2: str = ""
        i: int = self.seq1.getLength()
        j: int = self.seq2.getLength()
        pointer: str = traceback.getValue(i, j)  # type: ignore
        while i > 0 or j > 0:
            match pointer:
                case "DIAGONAL":
                    seq1 += self.seq1.getBase(i - 1)
                    seq2 += self.seq2.getBase(j - 1)
                    i -= 1
                    j -= 1
                case "UP":
                    seq1 += "-"
                    seq2 += self.seq2.getBase(j - 1)
                    j -= 1
                case "LEFT":
                    seq1 += self.seq1.getBase(i - 1)
                    seq2 += "-"
                    i -= 1
            pointer = traceback.getValue(i, j)  # type: ignore
        alignment: tuple[str, str] = self._reverseSeqs(seq1, seq2)
        return alignment

    def _scoreAlignment(
        self, alignment: tuple[str, str], annotation: str
    ) -> float:
        """Calculate score for global alignment."""
        seq1: str = alignment[0]
        seq2: str = alignment[1]
        score: float = 0
        for i in range(len(annotation)):
            char: str = annotation[i]
            if char == "|" or char == "*":
                key: tuple[str, str] = seq1[i], seq2[i]
                score += self.submatrix[key]
            else:
                if char == " ":
                    score += self.gap
        return score

    def execute(self, num: int, printOutput: int, path: str) -> None:
        """Run Needleman-Wunsch algorithm with linear scoring."""
        nrows: int = self.seq1.getLength() + 1
        ncols: int = self.seq2.getLength() + 1
        score: Matrix = self._createMatrix(nrows, ncols, "integer", "score")
        traceback: Matrix = self._createMatrix(
            nrows, ncols, "string", "traceback"
        )
        matrices: list[Matrix] = self._fillMatrices(score, traceback)
        alignment: tuple[str, str] = self._getTraceback(matrices[1])
        annotation: str = self._annotate(alignment)
        stats: list[float] = self._calcStats(alignment, annotation)
        if printOutput:
            self._print(num, stats, alignment, annotation)
        else:
            self._write(num, stats, alignment, annotation, path)


class Affine(NW):
    """A class to represent affine scoring global alignment."""

    def __init__(
        self, seq1: Sequence, seq2: Sequence, submatrix: SUB_MATRIX, gap: float
    ) -> None:
        """Construct all attributes for Affine."""
        super().__init__(seq1, seq2, submatrix, gap)

    def _initM(self, score: Matrix) -> Matrix:
        """Return initialized M matrix."""
        for i in range(1, score.ncols):  # first row
            score.setValue(self.gap + (i * self.extend), 0, i)
        for j in range(1, score.nrows):  # first column
            score.setValue(self.gap + (j * self.extend), j, 0)
        return score

    def _initI(self, score: Matrix) -> Matrix:
        """Return initialized I matrix."""
        for i in range(1, score.ncols):  # first row
            score.setValue(-100000, 0, i)
        for j in range(0, score.nrows):  # first column
            score.setValue(self.gap + (j * self.extend), j, 0)
        return score

    def _initD(self, score: Matrix) -> Matrix:
        """Return initialized D matrix."""
        for i in range(0, score.ncols):  # first row
            score.setValue(self.gap + (i * self.extend), 0, i)
        for j in range(0, score.nrows):  # first column
            score.setValue(-1000000, j, 0)
        return score

    def _initTM(self, trace: Matrix) -> Matrix:
        """Return initialized M Traceback matrix."""
        trace.setValue("STOP", 0, 0)
        for i in range(1, trace.ncols):  # first row
            trace.setValue("M:LEFT", 0, i)
        for j in range(1, trace.nrows):  # first column
            trace.setValue("M:UP", j, 0)
        return trace

    def _initTI(self, trace: Matrix) -> Matrix:
        """Return initialized I Traceback matrix."""
        trace.setValue("STOP", 0, 0)
        for i in range(1, trace.ncols):  # first row
            trace.setValue("I:LEFT", 0, i)
        for j in range(1, trace.nrows):
            trace.setValue("I:UP", j, 0)
        return trace

    def _initTD(self, trace: Matrix) -> Matrix:
        """Return initialized D Traceback matrix."""
        trace.setValue("STOP", 0, 0)
        for i in range(1, trace.ncols):  # first row
            trace.setValue("D:LEFT", 0, i)
        for j in range(1, trace.nrows):  # first column
            trace.setValue("D:UP", j, 0)
        return trace

    def _createMatrix(
        self, nrows: int, ncols: int, valueType: str, matType: str
    ) -> Matrix:
        """Return initialized score matrix."""
        matrix: Matrix = Matrix(nrows, ncols)
        matrix.initialize(valueType)
        match matType:
            case "M":
                matrix = self._initM(matrix)
            case "I":
                matrix = self._initI(matrix)
            case "D":
                matrix = self._initD(matrix)
            case "TM":
                matrix = self._initTM(matrix)
            case "TI":
                matrix = self._initTI(matrix)
            case "TD":
                matrix = self._initTD(matrix)
        return matrix

    def _traceM(self, scores: list[float], maxScore: float) -> str:
        """Return value for M traceback matrix."""
        if maxScore == scores[0]:
            return "M:DIAGONAL"
        if maxScore == scores[1]:
            return "I:DIAGONAL"
        if maxScore == scores[2]:
            return "D:DIAGONAL"
        return "ERROR"

    def _traceI(self, scores: list[float], maxScore: float) -> str:
        """Return value for I traceback matrix."""
        if maxScore == scores[0]:
            return "M:UP"
        if maxScore == scores[1]:
            return "I:UP"
        return "ERROR"

    def _traceD(self, scores: list[float], maxScore: float) -> str:
        """Return value for D traceback matrix."""
        if maxScore == scores[0]:
            return "M:LEFT"
        if maxScore == scores[1]:
            return "D:LEFT"
        return "ERROR"

    def _scoreLists(
        self, i: int, j: int, key: tuple[str, str], matrices: list[Matrix]
    ) -> dict[str, list[float]]:
        """Return score lists for M, I, and D matrices"""
        mismatch: Matrix = matrices[0]
        insert: Matrix = matrices[1]
        delete: Matrix = matrices[2]
        mScores: list[float] = [
            mismatch.getValue(i - 1, j - 1) + self.submatrix[key],  # type: ignore
            insert.getValue(i - 1, j - 1) + self.submatrix[key],  # type: ignore
            delete.getValue(i - 1, j - 1) + self.submatrix[key],  # type: ignore
        ]
        iScores: list[float] = [
            mismatch.getValue(i, j - 1) + self.gap,  # type: ignore
            insert.getValue(i, j - 1) + self.extend,  # type: ignore
        ]
        dScores: list[float] = [
            mismatch.getValue(i - 1, j) + self.gap,  # type: ignore
            insert.getValue(i - 1, j) + self.extend,  # type: ignore
        ]
        scoreLists: dict[str, list[float]] = {
            "M": mScores,
            "I": iScores,
            "D": dScores,
        }
        return scoreLists

    def _maxScores(self, scoreLists: dict[str, list[float]]) -> list[float]:
        """Return max scores for M, I and D matrices."""
        mismatch: float = max(scoreLists["M"])
        insert: float = max(scoreLists["I"])
        delete: float = max(scoreLists["D"])
        maxScores: list[float] = [mismatch, insert, delete]
        return maxScores

    def _updateScoreMats(
        self, i: int, j: int, mats: list[Matrix], values: list[float]
    ) -> None:
        """Update M, I, and D score matrices."""
        mats[0].setValue(values[0], i, j)
        mats[1].setValue(values[1], i, j)
        mats[2].setValue(values[2], i, j)

    def _updateTraceMats(
        self,
        i: int,
        j: int,
        matrices: list[Matrix],
        scoreLists: dict[str, list[float]],
        maxScores: list[float],
    ) -> None:
        """Update M, I, and D traceback matrices."""
        mValue: str = self._traceM(scoreLists["M"], maxScores[0])
        iValue: str = self._traceI(scoreLists["I"], maxScores[1])
        dValue: str = self._traceD(scoreLists["D"], maxScores[2])
        matrices[0].setValue(mValue, i, j)
        matrices[1].setValue(iValue, i, j)
        matrices[2].setValue(dValue, i, j)

    def _fillMatrices(
        self, scoreMats: list[Matrix], traceMats: list[Matrix]
    ) -> dict[str, list[Matrix]]:
        """Fill score and traceback matrices."""
        for i in range(1, scoreMats[0].nrows):
            for j in range(1, scoreMats[0].ncols):
                base1: str = self.seq1.getBase(i - 1)
                base2: str = self.seq2.getBase(j - 1)
                key: tuple[str, str] = base1, base2
                scoreLists: dict[str, list[float]] = self._scoreLists(
                    i, j, key, scoreMats
                )
                maxScores: list[float] = self._maxScores(scoreLists)
                self._updateScoreMats(i, j, scoreMats, maxScores)
                self._updateTraceMats(i, j, traceMats, scoreLists, maxScores)
        matrices: dict[str, list[Matrix]] = {
            "score": scoreMats,
            "traceback": traceMats,
        }
        return matrices

    def _parseMatrix(self, pointer: str) -> str:
        """Parse pointer for matrix."""
        idx: int = pointer.find(":")
        matrix: str = pointer[:idx]
        return matrix

    def _parseDirection(self, pointer: str) -> str:
        """Parse pointer for direction."""
        idx: int = pointer.find(":")
        start: int = idx + 1
        direction: str = pointer[start:]
        return direction

    def _getTraceback(self, traceMats: list[Matrix]) -> tuple[str, str]:
        """Get traceback for optimal alignment."""
        seq1: str = ""
        seq2: str = ""
        i: int = self.seq1.getLength()
        j: int = self.seq2.getLength()
        pointer: str = traceMats[0].getValue(i, j)  # type: ignore
        matrix: str = self._parseMatrix(pointer)
        direction: str = self._parseDirection(pointer)
        while i > 0 or j > 0:
            match direction:
                case "DIAGONAL":
                    seq1 += self.seq1.getBase(i - 1)
                    seq2 += self.seq2.getBase(j - 1)
                    i -= 1
                    j -= 1
                case "UP":
                    seq1 += "-"
                    seq2 += self.seq2.getBase(j - 1)
                    j -= 1
                case "LEFT":
                    seq1 += self.seq1.getBase(i - 1)
                    seq2 += "-"
                    i -= 1
            match matrix:
                case "M":
                    pointer = traceMats[0].getValue(i, j)  # type: ignore
                case "I":
                    pointer = traceMats[1].getValue(i, j)  # type: ignore
                case "D":
                    pointer = traceMats[2].getValue(i, j)  # type: ignore
            matrix = self._parseMatrix(pointer)
            direction = self._parseDirection(pointer)
        alignment: tuple[str, str] = self._reverseSeqs(seq1, seq2)
        return alignment

    def _scoreAlignment(
        self, alignment: tuple[str, str], annotation: str
    ) -> float:
        """Calculate score for global alignment."""
        seq1: str = alignment[0]
        seq2: str = alignment[1]
        score: float = 0
        for i in range(len(annotation)):
            char: str = annotation[i]
            if char == "|" or char == "*":
                key: tuple[str, str] = seq1[i], seq2[i]
                score += self.submatrix[key]
        gaps: int = annotation.count(" ")
        indels: int = self._calcIndels(annotation)
        score += self.gap * indels
        diff: int = gaps - indels
        score += self.extend * diff
        return score

    def execute(self, num: int, printOutput: int, path: str) -> None:
        """Run Needleman-Wunsch algorithm with linear scoring."""
        nrows: int = self.seq1.getLength() + 1
        ncols: int = self.seq2.getLength() + 1
        mismatch: Matrix = self._createMatrix(nrows, ncols, "integer", "M")
        insert: Matrix = self._createMatrix(nrows, ncols, "integer", "I")
        delete: Matrix = self._createMatrix(nrows, ncols, "integer", "D")
        mt: Matrix = self._createMatrix(nrows, ncols, "string", "TM")
        it: Matrix = self._createMatrix(nrows, ncols, "string", "TI")
        dt: Matrix = self._createMatrix(nrows, ncols, "string", "TD")
        scoreMats: list[Matrix] = [mismatch, insert, delete]
        traceMats: list[Matrix] = [mt, it, dt]
        matrices: dict[str, list[Matrix]] = self._fillMatrices(
            scoreMats, traceMats
        )
        alignment: tuple[str, str] = self._getTraceback(matrices["traceback"])
        annotation: str = self._annotate(alignment)
        stats: list[float] = self._calcStats(alignment, annotation)
        if printOutput:
            self._print(num, stats, alignment, annotation)
        else:
            self._write(num, stats, alignment, annotation, path)
