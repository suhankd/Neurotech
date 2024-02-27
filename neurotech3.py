import sys

def read_sequences(f1, f2):

	with open(f1) as file1, open(f2) as file2:
		return file1.read(), file2.read()

def needleman(seq1, seq2):

	match = 1
	mismatch = -1
	gap = -1

	rows = len(seq1)+1	
	cols = len(seq2)+1

	score = [[0]*cols for _ in range(rows)]

	for i in range(rows):
		score[i][0] = score[i-1][0] + gap
	for j in range(cols):
		score[0][j] = score[0][j-1] + gap

	for i in range(rows):

		for j in range(cols):

			match_mismatch = match if (seq1[i-1] == seq2[j-1]) else mismatch

			score[i][j] = max(score[i-1][j-1] + match_mismatch,
							  score[i-1][j] + gap,
							  score[i][j-1] + gap)

	align1, align2 = '', ''

	i,j = rows-1, cols-1

	while(i > 0 or j > 0):

		if i > 0 and j > 0:
		    if seq1[i-1] == seq2[j-1]:

		        current_score = score[i-1][j-1] + match

		    else:

		        current_score = score[i-1][j-1] + mismatch

		    if score[i][j] == current_score:

		        align1 = seq1[i-1] + align1
		        align2 = seq2[j-1] + align2
		        i -= 1
		        j -= 1

        elif i > 0 and score[i][j] == score[i-1][j] + gap:

            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1

        else:

            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, score[rows-1][cols-1]

if __name__ == "__main__":

	sequence1, sequence2 = read_sequences(sys.argv[1],sys.argv[2])

	alignment1, alignment2, score = needleman_wunsch(sequence1, sequence2)

    print("Alignment Score:", score)
    print("Sequence 1:", alignment1)
    print("Sequence 2:", alignment2)