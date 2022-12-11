import numpy as np
from abc import ABC, abstractmethod
from numba import njit

# TODO:
# - These were all done v. quickly and need refactoring
# - Aligning on AA level would be ~10X faster, but error is compounded

# --------------------------------------------------------------------------------
# Base class for performing pairwise sequence alignment
#
# --------------------------------------------------------------------------------


class PairwiseAligner(ABC):
    def set_sequences(self, x, y):
        self.x = x
        self.y = y
        self.n = len(x)
        self.m = len(y)
    
    @abstractmethod
    def set_scoring_model(self):
        pass
    
    @abstractmethod
    def align(self):
        pass


# --------------------------------------------------------------------------------
# Concrete implementations for nucleotide based alignments
#
# --------------------------------------------------------------------------------


class NeedlemanWunsch(PairwiseAligner):
    """
    Run global pairwise sequence alignment via Needleman-Wunsch

    Note:
    - Implemented with a linear (rather than affine) gap penalty
    
    """

    def set_scoring_model(self):
        # Scoring constants
        MATCH_SCORE = 2
        MISMATCH_SCORE = -3
        GAP_SCORE = -4  # linear score

        # Define the scoring model
        self.alphabet = ["A", "T", "C", "G"]
        self.l = len(self.alphabet)

        # Create substitution matrix
        self.sub_matrix = np.full((self.l, self.l), MISMATCH_SCORE)
        np.fill_diagonal(self.sub_matrix, MATCH_SCORE)

        # Set linear gap score
        self.gap_penalty = GAP_SCORE
        
    def align(self):
        """
        Align sequences `x` and `y`
        
        """
        
        # Suffix matrix
        self.F = np.zeros((self.n + 1, self.m + 1))
        self.F[0, 1:] = np.arange(1, self.m + 1) * self.gap_penalty
        self.F[1:, 0] = np.arange(1, self.n + 1) * self.gap_penalty
        
        # Traceback
        self.tb = np.zeros((self.n + 1, self.m + 1))
        self.tb[0, 1:] = 2
        self.tb[1:, 0] = 1
        
        # Aligned sequences
        self.x_aln = ""
        self.y_aln = ""
        self.score = None
        
        # Calculate suffix matrix
        for i in range(1, self.n + 1):
            for j in range(1, self.m + 1):

                # Compute scores of subsequences
                s = [ 
                    self.F[i - 1, j - 1] + self.sub_matrix[self.alphabet.index(self.x[i - 1]),
                                                           self.alphabet.index(self.y[j - 1])],
                    self.F[i - 1, j] + self.gap_penalty,  # Insertion in sequence y (indexed by j)
                    self.F[i, j - 1] + self.gap_penalty   # Insertion in sequence x (indexed by i)
                ]

                # Assign maximum score
                self.F[i, j] = max(s)

                # Trackback
                self.tb[i, j] = np.argmax(s)
        
        # Get the final score
        self.score = self.F[i, j]
        
        # Traceback
        i = self.n
        j = self.m
        while i > 0 or j > 0:
            z = self.tb[i, j]
            if z == 0:
                self.x_aln = self.x[i-1] + self.x_aln
                self.y_aln = self.y[j-1] + self.y_aln
                i -= 1
                j -= 1
            elif z == 1:  # Insertion in y
                self.x_aln = self.x[i-1] + self.x_aln
                self.y_aln = "-" + self.y_aln
                i -= 1
            elif z == 2:  # Insertion in x
                self.x_aln = "-" + self.x_aln
                self.y_aln = self.y[j-1] + self.y_aln
                j -=1

        return self.score



class NeedlemanWunschNumba(PairwiseAligner):
    """
    Run global pairwise sequence alignment via Needleman-Wunsch,
    Implemented with JIT-compilation via Numba,
    
    """

    def set_scoring_model(self):
       # Scoring constants
        self.MATCH_SCORE = 2
        self.MISMATCH_SCORE = -3
        self.GAP_SCORE = -4  # linear score
        
    def align(self):
        """
        Wrapper for Numba implementation

        I struggled a bit passing ndarray into this,
        but almost definitely sure it is possible
        
        """
        
        score = self._align(
            n=self.n,
            m=self.m,
            x=self.x,
            y=self.y,
            match_score = self.MATCH_SCORE,
            mismatch_score = self.MISMATCH_SCORE,
            gap_penalty=self.GAP_SCORE
        )
        self.score = score
    

    @staticmethod
    @njit
    def _align(n, m, x, y, match_score, mismatch_score, gap_penalty):
        """
        Pairwise alignment with Numba

        """      
        
        # Define alphabet
        alphabet = {"A":0,
                    "T":1,
                    "C":2,
                    "G":3}
        
        # Create substitution matrix
        sub_matrix = np.zeros((len(alphabet), len(alphabet)))
        for i in range(4):
            for j in range(4):
                if i == j:
                    sub_matrix[i, j] = match_score
                else:
                    sub_matrix[i, j] = mismatch_score

        # Suffix matrix
        F = np.zeros((n + 1, m + 1))
        F[0, 1:] = np.arange(1, m + 1) * gap_penalty
        F[1:, 0] = np.arange(1, n + 1) * gap_penalty
        
        # Only compute score, remove TB for speed
        score = None
        
        # Calculate suffix matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):

                # Compute scores of subsequences
                s = [ 
                    F[i - 1, j - 1] + sub_matrix[alphabet[x[i - 1]],
                                                 alphabet[y[j - 1]]],
                    F[i - 1, j] + gap_penalty,  # Insertion in sequence y (indexed by j)
                    F[i, j - 1] + gap_penalty   # Insertion in sequence x (indexed by i)
                ]

                # Assign maximum score
                F[i, j] = max(s)
        
        # Get the final score
        score = F[i, j]
   
        return score


class NeedlemanWunschNumbaBanded(PairwiseAligner):
    """
    Run global pairwise sequence alignment via Needleman-Wunsch,
    Implemented with JIT-compilation via Numba,
    Consider only alignments inside a fixed-width band

    """

    def set_scoring_model(self):
       # Scoring constants
        self.MATCH_SCORE = 2
        self.MISMATCH_SCORE = -3
        self.GAP_SCORE = -4  # linear score
        
    def align(self, band_radius=40):
        """
        Wrapper for Numba implementation

        I struggled a bit passing ndarray into this,
        but almost definitely sure it is possible
        
        """
        
        score = self._align(
            n=self.n,
            m=self.m,
            x=self.x,
            y=self.y,
            match_score = self.MATCH_SCORE,
            mismatch_score = self.MISMATCH_SCORE,
            gap_penalty=self.GAP_SCORE,
            band_radius=band_radius
        )
        self.score = score
    
    @staticmethod
    @njit
    def _align(n, m, x, y, match_score, mismatch_score, gap_penalty, band_radius): 
        """
        Pairwise alignment with Numba

        """             

        # Adjust banding for uneequal sequence lengths
        ratio = (m + 1) / (n + 1)
        
        # Define alphabet
        alphabet = {"A":0,
                    "T":1,
                    "C":2,
                    "G":3}
        
         # Create substitution matrix
        sub_matrix = np.zeros((len(alphabet), len(alphabet)))
        for i in range(4):
            for j in range(4):
                if i == j:
                    sub_matrix[i, j] = match_score
                else:
                    sub_matrix[i, j] = mismatch_score

        # Suffix matrix
        F = np.zeros((n + 1, m + 1))
        F[0, 1:] = np.arange(1, m + 1) * gap_penalty
        F[1:, 0] = np.arange(1, n + 1) * gap_penalty
        
        # Aligned sequences
        score = None
        
        # Calculate suffix matrix
        for i in range(1, n + 1):
            
            i_adj = int(i * ratio)
            jmin = max(0, i_adj - band_radius)
            jmax = min(m + 1, i_adj + band_radius + 1)

            for j in range(jmin, jmax):

                # Compute scores of subsequences
                s = [ 
                    F[i - 1, j - 1] + sub_matrix[alphabet[x[i - 1]],
                                                 alphabet[y[j - 1]]],
                    F[i - 1, j] + gap_penalty,  # Insertion in sequence y (indexed by j)
                    F[i, j - 1] + gap_penalty   # Insertion in sequence x (indexed by i)
                ]

                # Assign maximum score
                F[i, j] = max(s)
        
        # Get the final score
        score = F[i, j]
        
        return score

# --------------------------------------------------------------------------------
# Collection
#
# --------------------------------------------------------------------------------


ALIGNER_COLLECTION = {
    "needleman": NeedlemanWunsch,
    "needleman_numba": NeedlemanWunschNumba,
    "needleman_numba_banded": NeedlemanWunschNumbaBanded
}

