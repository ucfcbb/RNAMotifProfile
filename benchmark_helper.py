import numpy as np

# from profile import *
from search import *

def initiate_stat_file():
	fp = open('stat_data.csv', 'w')
	fp.write('Profile file name,Mean,STD\n')
	# fp.write(family_id + ',' + str(max_score) + ',' + str(max_s) + ' (' + str(round(max_s * 100 / max_score, 2)) + '%),' + str(min_s) + ',' + str(round(mean, 2)) + ',' + str(round(std, 2)) + '\n')
	fp.close()

def generate_stats_for_aligning_a_profile_with_source_loops(family_id, profile1, profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG, profile_fname, less_freq_nucl_and_intera_filter_threshold, frequent_gap_filtering_threshold):

	profile1 = filter_less_significant_nucl_and_interactions_from_a_profile(profile1, less_freq_nucl_and_intera_filter_threshold)
	profile1 = revise_profile_based_on_gap(profile1, frequent_gap_filtering_threshold)

	interaction_ind_pair_list1 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile1))
	profile1_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile1, scores.penalties)

	max_s = -1.0
	min_s = 10000000.0
	score_list = []
	for profile2 in profiles:	
		interaction_ind_pair_list2 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile2))
		interaction_compatibility_matrix = compute_interaction_compatibility(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, include_stackings)
		# print(interaction_compatibility_matrix)
		profile2_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile2, scores.penalties)
		# print('calling from 192')
		# max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, rank_dict, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
		max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
		max_score = get_scaled_score(max_score)

		max_s = max(max_s, max_score)
		min_s = min(min_s, max_score)

		score_list.append(max_score)


	interaction_ind_pair_list1 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile1))
	interaction_ind_pair_list2 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile1))
	interaction_compatibility_matrix = compute_interaction_compatibility(profile1, profile1, interaction_ind_pair_list1, interaction_ind_pair_list2, include_stackings)
	profile1_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile1, scores.penalties)
	profile2_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile1, scores.penalties)

	max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile1, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
	max_score = get_scaled_score(max_score)

	print(family_id)
	print('self score: ' + str(max_score))
	print('max with participating motifs: ' + str(max_s))
	print('min with participating motifs: ' + str(min_s))

	mean = np.mean(score_list)
	std = np.std(score_list)
	fp = open('stat_data.csv', 'a')
	fp.write(profile_fname + ',' + str(mean) + ',' + str(std) + '\n')
	# fp.write(family_id + ',' + str(max_score) + ',' + str(max_s) + ' (' + str(round(max_s * 100 / max_score, 2)) + '%),' + str(min_s) + ',' + str(round(mean, 2)) + ',' + str(round(std, 2)) + '\n')
	fp.close()

def load_stat_data(profile_fname):
	fp = open('stat_data.csv')
	lines = fp.readlines()
	fp.close()
	has_title_row = True


	mean_std_data_dict = {}
	if has_title_row:
		lines = lines[1:]
		
	for line in lines:
		pieces = line.strip().split(',')
		mean_std_data_dict[pieces[0].strip()] = (float(pieces[1].strip()), float(pieces[2].strip()))
	
	if profile_fname in mean_std_data_dict:
		return mean_std_data_dict[profile_fname]

	return None, None