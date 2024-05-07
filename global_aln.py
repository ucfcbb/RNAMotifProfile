def get_match_nucleotide_score_global_aln(i, j, seq1, seq2):
    if seq1[i] == seq2[j]:
        return 300
    return -100

def trace_back_alignment(seq1, seq2, DP_table_M, DP_table_Ix, DP_table_Iy):
    seq_len1 = len(seq1)
    seq_len2 = len(seq2)

    aln_P1 = ''
    aln_P2 = ''

    i = seq_len1
    j = seq_len2
    # print(i, j)
    # sys.exit()
    while (i>0 and j>0):
        match_score = get_match_nucleotide_score_global_aln(i-1, j-1, seq1, seq2)
        if (DP_table_M[i][j] == DP_table_M[i-1][j-1] + match_score):
            aln_P1 += seq1[i-1]
            aln_P2 += seq2[j-1]
            i -= 1
            j -= 1

        elif (DP_table_M[i][j] == DP_table_Ix[i-1][j-1] + match_score):
            # aln_P1 += '-'
            # aln_P2 += seq2[j-1]
            # j -= 1
            aln_P1 += seq1[i-1]
            aln_P2 += '-'
            i -= 1

        elif (DP_table_M[i][j] == DP_table_Iy[i-1][j-1] + match_score):
            # aln_P1 += seq1[i-1]
            # aln_P2 += '-'
            # i -= 1
            aln_P1 += '-'
            aln_P2 += seq2[j-1]
            j -= 1

    while (i>0):
        aln_P1 += seq1[i-1]
        aln_P2 += '-'
        i -= 1
    while (j>0):
        aln_P1 += '-'
        aln_P2 += seq2[j-1]
        j -= 1

    # print(aln_P1)
    # print(aln_P2)
    # aln_P1.reverse()
    # aln_P2.reverse()
    # print('reversed')
    # print(aln_P1)
    # print(aln_P2)
    return aln_P1[::-1], aln_P2[::-1]

def global_align_with_gap(seq1, seq2):
    Infinity = float('inf')

    gap_opening_penalty = -500
    gap_extension_penalty = -300

    seq_len1 = len(seq1)
    seq_len2 = len(seq2)

    rows, cols = (seq_len1 + 1, seq_len2 + 1)
    DP_table_M = [[0 for i in range(cols)] for j in range(rows)]
    DP_table_Ix = [[0 for i in range(cols)] for j in range(rows)]
    DP_table_Iy = [[0 for i in range(cols)] for j in range(rows)]

    DP_table_M[0][0] = 0
    DP_table_Ix[0][0] = gap_opening_penalty
    DP_table_Iy[0][0] = gap_opening_penalty

    for i in range(1, rows):
        DP_table_M[i][0] = -Infinity
        DP_table_Ix[i][0] = gap_opening_penalty + i * gap_extension_penalty
        DP_table_Iy[i][0] = -Infinity

    for j in range(1, cols):
        DP_table_M[0][j] = -Infinity
        DP_table_Ix[0][j] = -Infinity
        DP_table_Iy[0][j] = gap_opening_penalty + j * gap_extension_penalty

    for i in range(1, rows):
        for j in range(1, cols):

            DP_table_M[i][j] = get_match_nucleotide_score_global_aln(i-1, j-1, seq1, seq2) + max(
                    DP_table_M[i-1][j-1],
                    DP_table_Ix[i-1][j-1],
                    DP_table_Iy[i-1][j-1]
            )

            DP_table_Ix[i][j] = max(
                    gap_opening_penalty + gap_extension_penalty + DP_table_M[i-1][j],
                    # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Ix[i-1][j],
                    gap_extension_penalty + DP_table_Ix[i-1][j]
            )

            DP_table_Iy[i][j] = max(
                    gap_opening_penalty + gap_extension_penalty + DP_table_M[i][j-1],
                    gap_extension_penalty + DP_table_Iy[i][j-1]
                    # scores.penalties.gap_opening_penalty + scores.penalties.gap_extension_penalty + DP_table_Iy[i][j-1]
            )

    opt = max(DP_table_M[seq_len1][seq_len2], DP_table_Ix[seq_len1][seq_len2], DP_table_Iy[seq_len1][seq_len2])

    print(opt)
    print('tracing')
    aln_P1, aln_P2 = trace_back_alignment(seq1, seq2, DP_table_M, DP_table_Ix, DP_table_Iy)
    print(aln_P1)
    print(aln_P2)

def main():
    seq1 = 'CAAGCGAC'
    seq2 = 'UAAGCAA'

    global_align_with_gap(seq1, seq2)


if __name__ == "__main__":
    main()