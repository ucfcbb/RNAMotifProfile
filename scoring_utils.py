import copy
import sys

def add_a_region(aligned_regions, start, end):
    aligned_regions.append((start, end))

def extend_a_region(aligned_regions, start, end):
    # print(aligned_regions)
    if len(aligned_regions) == 0:
        add_a_region(aligned_regions, start, end)
        return
    last_ind = aligned_regions[-1][1]
    if last_ind + 1 != start:
        print("Error")
    s, e = aligned_regions[-1]
    aligned_regions[-1] = (s, end)
    # aligned_regions[-1][1] = end

def is_indices_across_break_points(break_points, a, b):
    for k in range(len(break_points)):
        if (a < break_points[k] and b >= break_points[k]) or (b < break_points[k] and a >= break_points[k]):
            return True
    return False

def match_isosteric_basepair(i, j, basepair_index, scores, isosteric_substitution_matrix):
    if i >= 0 and j >= 0 and basepair_index[i] > 0 and basepair_index[j] > 0:

        i1 = (i - 1) / 16
        j1 = (j - 1) / 16

        iori = i1 % 2
        jori = j1 % 2

        i2 = i1 / 2
        j2 = j1 / 2

        iedga = i2 % 3
        iedgb = i2 / 3
        jedga = j2 % 3
        jedgb = j2 / 3

        flip_bonus = 0
        if iori == jori and iedga != jedga and iedgb != jedgb and iedga == jedgb and iedgb == jedga:
            flip_bonus = (isosteric_substitution_matrix[basepair_index[i] - 1][basepair_index[j] - 1] * 0.5)

        return isosteric_substitution_matrix[basepair_index[i] - 1][basepair_index[j] - 1] + flip_bonus;

    elif i >= 0 and j >= 0 and (basepair_index[i] == 0 or basepair_index[j] == 0):
        return scores.hbond_match_base

    elif i >= 0 and j >= 0 and (basepair_index[i] == 0 and basepair_index[j] == 0):
        return scores.hbond_match_base + scores.bonuses.hbond_match_bonus

    return 0

# def match_isosteric_basepair_profile(bp_dict1, bp_dict2, scores, isosteric_substitution_matrix):
#     bp_cnt1 = len(bp_dict1)
#     bp_cnt2 = len(bp_dict2)

#     print(bp_dict1)
#     print(bp_dict2)

#     for h_bp1 in bp_dict1:
#         weighted_score = 0.0
#         for h_bp2 in bp_dict2:
#             weighted_score += bp_dict2[h_bp2] * match_isosteric_basepair(h_bp1, h_bp2, scores, isosteric_substitution_matrix) / bp_cnt2
#         weighted_score_list.append(weighted_score)

#     return sum(weighted_score_list) / len(weighted_score_list)


def match_isosteric_basepair_profile(hbp_list1, hbp_list2, basepair_index, scores, isosteric_substitution_matrix):
    bp_cnt1 = len(hbp_list1)
    bp_cnt2 = len(hbp_list2)

    # print(hbp_list1)
    # print(hbp_list2)

    weighted_score_list = []
    for h_bp1 in hbp_list1:
        weighted_score = 0.0
        for h_bp2 in hbp_list2:
            weighted_score += hbp_list2.count(h_bp2) * match_isosteric_basepair(h_bp1, h_bp2, basepair_index, scores, isosteric_substitution_matrix) / bp_cnt2
        weighted_score_list.append(weighted_score)

    return sum(weighted_score_list) / len(weighted_score_list)





def match_stacking(i, j, scores, stacking_substitution_matrix):
    #    return the corresponding score of matching two basestackings

    if i < 0 and j < 0:
        i = abs(i) - 1
        j = abs(j) - 1
        return stacking_substitution_matrix[i][j]

    elif (i == 0 and j < 0) or (i < 0 and j == 0):
        #    note that we do not consider an hydrogen bond match with a base stacking
        return 0
    
    elif i == 0 and j == 0:
        # cout << "ERROR!! H-bond match call in match_stacking_basepair(...)!!" << endl;
        # exit(0);
        return scores.hbond_match_base + scores.bonuses.hbond_match_bonus
    
    return 0;

# def match_stacking_profile(stk_dict1, stk_dict2, scores, stacking_substitution_matrix):
#     stk_cnt1 = len(stk_dict1)
#     stk_cnt2 = len(stk_dict2)

#     weighted_score_list = []
#     for h_stk1 in stk_dict1:
#         weighted_score = 0.0
#         for h_stk2 in stk_dict2:
#             weighted_score += stk_dict2[h_stk2] * match_stacking(h_stk1, h_stk2, scores, stacking_substitution_matrix) / stk_cnt2
#         weighted_score_list.append(weighted_score)

#     return sum(weighted_score_list) / len(weighted_score_list)

def match_stacking_profile(hstk_list1, hstk_list2, scores, stacking_substitution_matrix):
    stk_cnt1 = len(hstk_list1)
    stk_cnt2 = len(hstk_list2)

    weighted_score_list = []
    for h_stk1 in hstk_list1:
        weighted_score = 0.0
        for h_stk2 in hstk_list2:
            weighted_score += hstk_list2.count(h_stk2) * match_stacking(h_stk1, h_stk2, scores, stacking_substitution_matrix) / stk_cnt2
        weighted_score_list.append(weighted_score)

    return sum(weighted_score_list) / len(weighted_score_list)

def get_adjacent_stacking_score(b1, e1, b2, e2, stksP1, stksP2, scores, stacking_substitution_matrix):
    if (b1, e1) in stksP1 and (b2, e2) in stksP2:
        return scores.weight_adjacent_stacking * match_stacking_profile(stksP1[(b1, e1)], stksP2[(b2, e2)], scores, stacking_substitution_matrix)
    return 0

    # if(b1 < 0 || e1 >= (int)seq1.length() || b2 < 0 || e2 >= (int)seq2.length())
    #     return 0;
    # else
    # {
    #     # Condition added to stop calling match_stacking_basepair(0, 0) [Edited by Shahidul]
    #     # It used to return score for h-bond match, while (0, 0) actually represent the initialization value (it has nothing to do with h-bond)
    #     if(!INCLUDE_HBOND_IN_ALL_PAIR_STACK)
    #       if(M1_effective_stacking[b1][e1] == 0 && M2_effective_stacking[b2][e2] == 0)
    #         return 0;
    #     return scoring_function.weight_adjacent_stacking
    #     * scoring_function.match_stacking_profile(M1_effective_stacking[b1][e1], M2_effective_stacking[b2][e2]);
    # }



def match_nucleotide(a, b, nucleotide_substitution_matrix):
    #    return the corresponding score of matching two nucleotides
    if a == -1 or b == -1:
        return 0
    return nucleotide_substitution_matrix[a][b]

def match_nucleotide_profile(nucl_dict1, nucl_dict2, nucleotide_substitution_matrix):
    nucl_cnt1 = len(nucl_dict1)
    nucl_cnt2 = len(nucl_dict2)

    weighted_score_list = []
    for h_nucl1 in nucl_dict1:
        weighted_score = 0.0
        for h_nucl2 in nucl_dict2:
            weighted_score += nucl_dict2[h_nucl2] * match_nucleotide(h_nucl1, h_nucl2, nucleotide_substitution_matrix) / nucl_cnt2
        weighted_score_list.append(weighted_score)

    return sum(weighted_score_list) / len(weighted_score_list)

def compute_nucleotide_deletion_penalties(profile, penalties):
    l = len(profile.sequenceP)
    rows, cols = (l, l)
    nucleotide_deletion_penalties = [[0 for i in range(cols)] for j in range(rows)]

    # for i in range(l):
    #     for j in range(i+2, l):
    #         k = i + 2
    #         while k - i < l:
    #             j = k % l
    #             print(i, j)
    #             nucleotide_deletion_penalties[i][j] += penalties.gap_opening_penalty + penalties.gap_extension_penalty * (k - i - 1)

    for i in range(l-2):
        for j in range(i+2, l):
            if is_indices_across_break_points(profile.break_points, i, j) == True:
                continue

            nucleotide_deletion_penalties[i][j] += (penalties.gap_opening_penalty + penalties.gap_extension_penalty * (j - i - 1))
            nucleotide_deletion_penalties[j][i] = nucleotide_deletion_penalties[i][j]
            # k += 1

    # print(nucleotide_deletion_penalties)
    return nucleotide_deletion_penalties

def get_nucleotide_deletion_penalty(b, e, scores, nucleotide_deletion_penalties, seqP):

    if b < 0:
        b = 0
    if e >= len(seqP):
        e = len(seqP) - 1
    if e >= b:
        return scores.weight_sequence * nucleotide_deletion_penalties[b][e]
    else:
        return 0

def  get_match_nucleotide_score(i, j, seqP1, seqP2, scores, nucleotide_substitution_matrix):
    if i < 0 or i >= len(seqP1) or j < 0 or j >= len(seqP2):
        return 0
    else:
        return scores.weight_sequence * match_nucleotide_profile(seqP1[i], seqP2[j], nucleotide_substitution_matrix)


def is_isosteric(interaction_ind_pair1, interaction_ind_pair2, basepair_index):
    if interaction_ind_pair1[2] == 'b' and interaction_ind_pair2[2] == 'b' and \
    len(interaction_ind_pair1[3]) == 1 and len(interaction_ind_pair2[3]) == 1 and \
    interaction_ind_pair1[3].keys()[0] < 99999 and interaction_ind_pair2[3].keys()[0] < 99999 and \
    basepair_index[interaction_ind_pair1[3].keys()[0]] == basepair_index[interaction_ind_pair2[3].keys()[0]]:
        # if both are canonical base pair, then they does not counted as isosteric matching
        if basepair_index[interaction_ind_pair1[3].keys()[0]] >= 2:
            return True
    return False

def asymmetric_penalty(matched_pairs, interaction_ind_pair_list1, interaction_ind_pair_list2, scores):
    penalty_score = 0.0

    # sorting matched_pairs
    matched_pairs_size = len(matched_pairs)
    for i in range(matched_pairs_size - 1):
        for j in range(i+1, matched_pairs_size):
            # Sort based on the info of the base pairs from the align entry 1 (sequence or partial profile) (???)
            if (interaction_ind_pair_list1[matched_pairs[i][0]][0] > interaction_ind_pair_list1[matched_pairs[j][0]][0]) or \
            (interaction_ind_pair_list1[matched_pairs[i][0]][0] == interaction_ind_pair_list1[matched_pairs[j][0]][0] and \
            interaction_ind_pair_list1[matched_pairs[i][0]][1] > interaction_ind_pair_list1[matched_pairs[j][0]][1]):
                temp = matched_pairs[i]
                matched_pairs[i] = matched_pairs[j]
                matched_pairs[j] = temp
    

    for bp_count in range(1, matched_pairs_size):
        p1 = matched_pairs[bp_count - 1][0] # Index of the bp from align entry 1 (sequence or partial profile)
        p2 = matched_pairs[bp_count - 1][1] # Index of the bp from align entry 2 (sequence or partial profile)

        #The stack condition may be unneccesary. Change to allow stack (by Shahidul)
        if interaction_ind_pair_list1[p1] != 'b':
            continue;

        #Record where previous interation is in the adjacent enclosed one (???)
        #This code block might work for descending order sorted list (Shahidul)
        record_i = 0
        for i in range(bp_count, matched_pairs_size):
            #Compare BP indices
            #If i-th BP is enclosed by bp_count-th BP
            if interaction_ind_pair_list1[matched_pairs[i][0]][0] > interaction_ind_pair_list1[p1][0] and \
            interaction_ind_pair_list1[matched_pairs[i][0]][1] < interaction_ind_pair_list1[p1][1]:
                #Recoding the first BP that is enclosed inside
                record_i = i
                break

        
        adjacent_pairs = []
        adjacent_pairs.append(record_i)

        max_index = interaction_ind_pair1[matched_pairs[record_i][0]][1]

        #Decision: this code block is fine. It's the BP order in matched that messed things up (Shahidul)
        for i in range(record_i+1, matched_pairs_size):
            if interaction_ind_pair_list1[matched_pairs[i][0]][1] > max_index and interaction_ind_pair_list1[matched_pairs[i][0]][1] < interaction_ind_pair_list1[p1][1]:
                max_index = interaction_ind_pair1[matched_pairs[i][0]][1]
                adjacent_pairs.append(i)
        

    

        #  compute asymmetricity between the loops
        adjacent_pairs_size = len(adjacent_pairs)
        for i in range(adjacent_pairs_size):
            k = adjacent_pairs[i]
        
            sp1_l = interaction_ind_pair1[matched_pairs[k][0]][0] - interaction_ind_pair_list1[p1][0]
            sp1_r = interaction_ind_pair1[p1][1] - interaction_ind_pair_list1[matched_pairs[k][0]][1]
            
            sp2_l = (interaction_ind_pair_list2[matched_pairs[k][1]][0] - interaction_ind_pair_list2[p2][0])
            sp2_r = (interaction_ind_pair_list2[p2][1] - interaction_ind_pair_list2[matched_pairs[k][1]][1])

            a1 = sp1_l - sp1_r
            a2 = sp2_l - sp2_r
            asym_size = abs(a1 - a2)

            penalty_score += scores.weight_isosteric * asym_size * scores.asym_nuc

            if a1 * a2 < 0:
                penalty_score += scores.weight_isosteric * scores.asym_loop

            if sp1_l == 1 and sp1_r == 1:
                #  in case of a stacking case
                if sp2_l == 1 and sp2_r == 1 and \
                is_isosteric(interaction_ind_pair1[matched_pairs[k][0]], interaction_ind_pair2[matched_pairs[k][1]]) and \
                is_isosteric(interaction_ind_pair1[p1], interaction_ind_pair_list2[p2]):
                    penalty_score += scores.weight_isosteric * scores.cons_stacking

                elif (sp2_l > 1 or sp2_r > 1) and sp2_l != sp2_r:
                    #  asymmetric internal loop or bulge
                    if sp2_l == 1 or sp2_r == 1:
                        #  a bulge
                        penalty_score += scores.weight_isosteric * scores.stack_to_bulge
                    else:
                        #  an aymmetric internal loop
                        penalty_score += scores.weight_isosteric * scores.stack_to_internal_asym
                    
                elif (sp2_l > 1 or sp2_r > 1) and sp2_l == sp2_r:
                    #  symmetric internal loop
                    penalty_score += scores.weight_isosteric * scores.stack_to_internal_sym
                
    return penalty_score

def count_deleted_interactions(aligned_pairs, interaction_ind_pair_list, aligned_regions, count_single_end=False):
    # initialization
    num_bp_deleted = 0
    num_stk_deleted = 0

    print('aligned_pairs')
    print(aligned_pairs)
    # for each base pair that is not aligned
    for i in range(len(aligned_pairs)):
        # Test if the annotation is unmatched
        if aligned_pairs[i] == 0:
            interaction_end_in_the_region = 0
            print('aligned_regions')
            print(aligned_regions)
            for (s, e) in aligned_regions:
                if interaction_ind_pair_list[i][0] >= s and interaction_ind_pair_list[i][0] <= e:
                    interaction_end_in_the_region += 1
                if interaction_ind_pair_list[i][1] >= s and interaction_ind_pair_list[i][1] <= e:
                    interaction_end_in_the_region += 1

                if interaction_end_in_the_region == 2:
                    break
            # Check if at least on end in the region
            # if(interaction_end_in_the_region == 1)  {
            # Check if both ends in the region
            print('interaction_end_in_the_region')
            print(interaction_end_in_the_region)
            if (count_single_end == False and interaction_end_in_the_region == 2) or \
            (count_single_end == True and interaction_end_in_the_region == 1):
                # Check if it is base pair (not h_bond)
                # print(interaction_ind_pair_list)
                if interaction_ind_pair_list[i][2] == 'b' and len(interaction_ind_pair_list[i][3]) == 1 and interaction_ind_pair_list[i][3][0] < 99999:
                    num_bp_deleted += 1
                # Check if it is stack
                elif interaction_ind_pair_list[i][2] == 's':
                    num_stk_deleted += 1
                # What if it is h-bond?
                # Is it ignored as a potential mistake to identify a weak bond?

    return num_bp_deleted, num_stk_deleted

def get_missing_bp_and_stack_penalty(interaction_ind_pair_list1, interaction_ind_pair_list2, aligned_regions1, aligned_regions2, pairs_record1, pairs_record2, scores):
    penalty_score = 0.0

    num_bp_deleted1, num_stk_deleted1 = count_deleted_interactions(pairs_record1, interaction_ind_pair_list1, aligned_regions1)
    print('num_bp_deleted1')
    print(num_bp_deleted1)
    print('num_stk_deleted1')
    print(num_stk_deleted1)
    num_bp_deleted2, num_stk_deleted2 = count_deleted_interactions(pairs_record2, interaction_ind_pair_list2, aligned_regions2)
    print('num_bp_deleted2')
    print(num_bp_deleted2)
    print('num_stk_deleted2')
    print(num_stk_deleted2)

    penalty_score += scores.weight_isosteric * (scores.penalties.missing_bp_penalty * num_bp_deleted1) \
    + scores.weight_nonadjacent_stacking * (scores.penalties.missing_stk_penalty * num_stk_deleted1)
    penalty_score += scores.weight_isosteric * (scores.penalties.missing_bp_penalty * num_bp_deleted2) \
    + scores.weight_nonadjacent_stacking * (scores.penalties.missing_stk_penalty * num_stk_deleted2)

    return penalty_score

def bp_and_stack_penalty_inclusion_exclusion(break_points1, break_points2, interaction_ind_pair_list1, interaction_ind_pair_list2, k, i, l, j, offset1, offset2, current_region_i, current_region_j, align_region_i, align_region_j, scores):
    # vector<int> M1_pairs_record(M1_interaction_union.size(), 0)
    # vector<int> M2_pairs_record(M2_interaction_union.size(), 0)

    pairs_record1 = [0 for i in range(len(interaction_ind_pair_list1))]
    pairs_record2 = [0 for i in range(len(interaction_ind_pair_list2))]

    current_region_i = copy.deepcopy(align_region_i[k][l])
    current_region_j = copy.deepcopy(align_region_j[k][l])

    if is_indices_across_break_points(break_points1, offset1 + k, offset1 + i) or is_indices_across_break_points(break_points2, offset2 + l, offset2 + j):
        add_a_region(current_region_i, offset1 + i, offset1 + i)
        add_a_region(current_region_j, offset2 + j, offset2 + j)
    else:
        extend_a_region(current_region_i, offset1 + k + 1, offset1 + i)
        extend_a_region(current_region_j, offset2 + l + 1, offset2 + j)

    # # if (is_indices_across_break_points(M1_break_points[0], offset1 + k, offset1 + i) ||
    # # is_indices_across_break_points(M2_break_points[0], offset2 + l, offset2 + j)):
    # #     add_a_region(current_region_i, offset1 + i, offset1 + i);
    # #     add_a_region(current_region_j, offset2 + j, offset2 + j);
    # # else:
    # # print('from scoring_utils')
    # extend_a_region(current_region_i, offset1 + k + 1, offset1 + i)
    # # print('from scoring_utils')
    # extend_a_region(current_region_j, offset2 + l + 1, offset2 + j)

    missing_interaction_penalty = get_missing_bp_and_stack_penalty(interaction_ind_pair_list1, interaction_ind_pair_list2, current_region_i, current_region_j, pairs_record1, pairs_record2, scores)
    print('missing_interaction_penalty')
    print(missing_interaction_penalty)
    missing_interaction_penalty -= get_missing_bp_and_stack_penalty(interaction_ind_pair_list1, interaction_ind_pair_list2, align_region_i[k][l], align_region_j[k][l], pairs_record1, pairs_record2, scores)
    print('missing_interaction_penalty')
    print(missing_interaction_penalty)
    # sys.exit()
    return missing_interaction_penalty