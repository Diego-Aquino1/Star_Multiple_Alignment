import numpy as np
from datetime import datetime

MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_PENALTY = -2

def remove_spaces(sequence):
    return sequence.replace(" ", "")

def generate_file(matrix, alignments, multiple, star, filename, time):
    with open(filename, 'w') as file:
        file.write(f"Execution time: {time}\n")
        file.write(f"Best score: {star[1]}\n")
        file.write("\nMatrix with all the scores from the alignments:\n")
        for row in matrix:
            file.write("\t".join(map(str, row)) + "\n")
        file.write("\nAlignments with star sequence:\n")
        for i, alignment in enumerate(alignments):
            file.write(f"Alignment S{star[0] + 1} - S{i+1}\n")
            file.write(f"{alignment[0]}\n")
            file.write(f"{alignment[1]}\n")
        file.write("\nMultiple Alignment:\n")
        for seq in multiple:
            file.write(f"{seq}\n")
        print(f"Results saved to {filename}")

class NeedlemanWunsch:
    def __init__(self, seq1, seq2, match_score=MATCH_SCORE, mismatch_score=MISMATCH_SCORE, gap_penalty=GAP_PENALTY):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty
        self.m = len(seq1)
        self.n = len(seq2)
        self.score_matrix = np.zeros((self.m + 1, self.n + 1), dtype=int)
        self.traceback_matrix = [[[] for _ in range(self.n + 1)] for _ in range(self.m + 1)]
        self.alignments = []

    def align(self):
        self.initialize_matrices()
        self.fill_matrices()
        self.traceback(self.m, self.n, '', '')

    def initialize_matrices(self):
        for i in range(1, self.m + 1):
            self.score_matrix[i][0] = self.score_matrix[i - 1][0] + self.gap_penalty
            self.traceback_matrix[i][0] = [(i - 1, 0)]
        for j in range(1, self.n + 1):
            self.score_matrix[0][j] = self.score_matrix[0][j - 1] + self.gap_penalty
            self.traceback_matrix[0][j] = [(0, j - 1)]

    def fill_matrices(self):
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                match = self.score_matrix[i - 1][j - 1] + (self.match_score if self.seq1[i - 1] == self.seq2[j - 1] else self.mismatch_score)
                delete = self.score_matrix[i - 1][j] + self.gap_penalty
                insert = self.score_matrix[i][j - 1] + self.gap_penalty
                max_score = max(match, delete, insert)
                self.score_matrix[i][j] = max_score
                if match == max_score:
                    self.traceback_matrix[i][j].append((i - 1, j - 1))
                if delete == max_score:
                    self.traceback_matrix[i][j].append((i - 1, j))
                if insert == max_score:
                    self.traceback_matrix[i][j].append((i, j - 1))

    def traceback(self, i, j, alignment1, alignment2):
        if i == 0 and j == 0:
            self.alignments.append((alignment1[::-1], alignment2[::-1]))
            return
        for prev_i, prev_j in self.traceback_matrix[i][j]:
            if (prev_i, prev_j) == (i - 1, j - 1):
                self.traceback(prev_i, prev_j, alignment1 + self.seq1[i - 1], alignment2 + self.seq2[j - 1])
            elif (prev_i, prev_j) == (i - 1, j):
                self.traceback(prev_i, prev_j, alignment1 + self.seq1[i - 1], alignment2 + '-')
            elif (prev_i, prev_j) == (i, j - 1):
                self.traceback(prev_i, prev_j, alignment1 + '-', alignment2 + self.seq2[j - 1])

def find_score(sequence1, sequence2):
    nw = NeedlemanWunsch(sequence1, sequence2)
    nw.align()
    return nw.score_matrix, nw.alignments[-1]

def calculate_max_score(num_seq, sequences):
    all_scores = np.zeros((num_seq, num_seq), dtype=int)
    max_score = float('-inf')
    star = 0
    for i in range(num_seq):
        sum_scores = 0
        for j in range(i + 1, num_seq):
            score_matrix, _ = find_score(sequences[i], sequences[j])
            score = score_matrix[-1, -1]
            all_scores[i][j] = all_scores[j][i] = score
            sum_scores += score
        if sum_scores > max_score:
            max_score = sum_scores
            star = i

    row_sums = np.sum(all_scores, axis=1)

    max_row_position = np.argmax(row_sums)

    max_sum = np.max(row_sums)
    
    return max_row_position, max_score, all_scores, max_sum


# def star_alignment(num_seq, sequences):
#     star, _, all_scores = calculate_max_score(num_seq, sequences)
#     alignments = []
#     max_len = 0

#     # Find all alignments of the star sequence with the other sequences
#     for i in range(num_seq):
#         if i != star:
#             _, alignment = find_score(sequences[star], sequences[i])
#             alignments.append(alignment)
#             max_len = max(max_len, len(alignment[0]), len(alignment[1]))

#     multiple_alignment = ['-' * max_len for _ in range(num_seq)]
#     multiple_alignment[star] = sequences[star].ljust(max_len, '-')

#     for i in range(num_seq):
#         if i != star:
#             alignment1, alignment2 = alignments[i if i < star else i - 1]
#             padded_alignment1 = alignment1.ljust(max_len, '-')
#             padded_alignment2 = alignment2.ljust(max_len, '-')
#             multiple_alignment[star] = padded_alignment1

#             for j in range(max_len):
#                 if multiple_alignment[i][j] == '-' and padded_alignment2[j] != '-':
#                     multiple_alignment[i] = multiple_alignment[i][:j] + padded_alignment2[j] + multiple_alignment[i][j + 1:]
#                 elif multiple_alignment[i][j] != '-' and padded_alignment2[j] == '-':
#                     multiple_alignment[i] = multiple_alignment[i][:j] + '-' + multiple_alignment[i][j:]

#     return multiple_alignment, alignments, star, all_scores

def star_alignment(num_seq, sequences):
    star, _, all_scores, max_score = calculate_max_score(num_seq, sequences)
    alignments = []
    max_len = 0

    # Find all alignments of the star sequence with the other sequences
    for i in range(num_seq):
        if i != star:
            _, alignment = find_score(sequences[star], sequences[i])
            alignments.append(alignment)
            max_len = max(max_len, len(alignment[0]), len(alignment[1]))

    # Initialize multiple_alignment with the star sequence
    multiple_alignment = ['-' * max_len for _ in range(num_seq)]
    multiple_alignment[star] = sequences[star].ljust(max_len, '-')

    for i in range(num_seq):
        if i != star:
            alignment1, alignment2 = alignments[i if i < star else i - 1]
            padded_alignment1 = alignment1.ljust(max_len, '-')
            padded_alignment2 = alignment2.ljust(max_len, '-')

            # Create an aligned sequence with gaps
            aligned_seq = ['-'] * max_len
            star_index = 0

            for j in range(max_len):
                if padded_alignment2[j] == '-':
                    star_index += 1
                    aligned_seq[j] = '-'
                else:
                    aligned_seq[j] = padded_alignment2[star_index]
                    star_index += 1

            multiple_alignment[i] = ''.join(aligned_seq)

    return multiple_alignment, alignments, star, all_scores, max_score


def main():
    sequences = [
        "attaaaggtttataccttcc",
        "caggtaacaaaccaaccaac",
        "tttcgatctcttgtagatct",
        "gttctctaaacgaactttaa",
        "aatctgtgtggctgtcactc",
        "ggctgcatgcttagtgcact"
    ]
    filename = "alignment_results_3.txt"

    start_time = datetime.now()
    multiple_alignment, alignments, star, all_scores, max_score = star_alignment(len(sequences), sequences)
    end_time = datetime.now()

    execution_time = (end_time - start_time).total_seconds()
    generate_file(all_scores, alignments, multiple_alignment, (star, max_score), filename, execution_time)

if __name__ == "__main__":
    main()
