import sys
import os
import time
import argparse
import logging
import multiprocessing as mp

from logger import *
from utils import *
from input_utils import *
from interaction_utils import *
from scoring_utils import *

def generate_single_profile_multi_process(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, mp_number_of_process, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG):
	total = len(profiles)

	if total == 1:
		return profiles[0]

	item_per_piece = items_per_chunk
	if total > item_per_piece:
		# print(total)
		# total_pieces = total // item_per_piece if total % item_per_piece == 0 else total // item_per_piece + 1
		# print(total_pieces)
		logger.info('- - ' + str(total) + ' items to align, distributing them into small chunks of ' + str(items_per_chunk) + ' items.')

		manager = mp.Manager()
		temp_profiles = manager.list()

		# temp_profiles = []
		parameter_list = []
		k = 0
		for ii in range(0, total, item_per_piece):
			temp_profiles.append(None)

		for ii in range(0, total, item_per_piece):
			st = ii
			fi = min(ii + item_per_piece - 1, total)
			# print('generating profiles for:')
			# print(st, fi)
			# p = generate_single_profile(profiles[st:fi+1], scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index)
			# temp_profiles.append(p)
			# temp_profiles.append(None)
			parameter_list.append((profiles[st:fi+1], scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, temp_profiles, k, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG))
			# ii += 5
			k += 1

		pool = mp.Pool(mp_number_of_process)
		pool.map(_generate_single_profile_worker, parameter_list)


		# print(temp_profiles)
		# sys.exit()
		while check_if_all_profiles_are_generated(temp_profiles) == False:
			print('Waiting for previous stage to be completed')
			time.sleep(10)

		if check_if_all_profiles_are_generated(temp_profiles):
			return generate_single_profile_multi_process(temp_profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, mp_number_of_process, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)

	return generate_single_profile(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)

def check_if_all_profiles_are_generated(temp_profiles):
	for p in temp_profiles:
		if p == None:
			return False
	return True

def _generate_single_profile_worker(p):
	generate_single_profile_mp(*p)

def generate_single_profile_mp(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, temp_profiles, ind, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG):
	p = generate_single_profile(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)
	temp_profiles[ind] = p
	# print('assigning ' + p + ' to ind: ' + str(ind))
	# return

def generate_single_profile(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG):
	# return "dummy"
	total = len(profiles)
	if dBG:
		print('received ' + str(total) + ' profiles')
		print('loop counts are following')
		for profile in profiles:
			print(profile.loopCount)

	if total == 1:
		if dBG:
			print('returning')
		return profiles[0]

	item_per_piece = items_per_chunk
	if total > item_per_piece:
		logger.info('- - ' + str(total) + ' items to align, distributing them into small chunks of ' + str(items_per_chunk) + ' items.')
		# if dBG:
		# 	# print(total)
		# 	logger.debug('Total ' + str(total) + ' items found. Distributing them into multiple chunks of ' + str(items_per_chunk) + ' items.')
		# total_pieces = total // item_per_piece if total % item_per_piece == 0 else total // item_per_piece + 1
		# print(total_pieces)

		temp_profiles = []
		for ii in range(0, total, item_per_piece):
			st = ii
			fi = min(ii + item_per_piece - 1, total)
			if dBG:
				# print('generating profiles for:')
				logger.debug('generating profiles for the items in the following indices range: ' + str((st, fi)))
				# print(st, fi)
				# print('sending')
				# print(fi-st+1)
			p = generate_single_profile(profiles[st:fi+1], scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)
			temp_profiles.append(p)
			# print('found')
			# print(p.loopCount)
			# break
			# ii += 5
		# if len(temp_profiles) <= 3:
		# 	print('switching into debug mode')
		# 	dBG = True
		return generate_single_profile(temp_profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)
		# sys.exit()
		# for ii in range(total/5)

	# logger.info(str(total) + ' items to align.')
	max_score, maximal_clique, ind1, ind2, profile1, profile2 = find_best_aligned_profile_pair(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
	
	if dBG:
		print('Best profile pair index:')
		print(ind1, ind2)
		print('loop counts: ')
		print(profile1.loopCount, profile2.loopCount)

	if dBG:
		logger.debug('max_score: ' + str(max_score))
		# print('max_score')
		# print(max_score)
		logger.debug('maximal_clique: ' + str(maximal_clique))
		# print('maximal_clique')
		# print(maximal_clique)
		logger.debug('Best profile pair index: ' + str((ind1, ind2)))
		# print('Best profile pair index:')
		# print(ind1, ind2)
		logger.debug(str((profiles[ind1].loopName, profiles[ind2].loopName)))
		# print(profiles[ind1].loopName, profiles[ind2].loopName)
		print_formatted_profile(profiles[ind1])
		print_formatted_profile(profiles[ind2])

	# profile1 = profiles[ind1]
	# profile2 = profiles[ind2]

	interaction_ind_pair_list1 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile1))
	interaction_ind_pair_list2 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile2))


	# interaction_compatibility_matrix = compute_interaction_compatibility(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2)
	# profile1_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile1, scores.penalties)
	# profile2_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile2, scores.penalties)
	# # print('calling from 136')
	# max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, dBG)
	# if dBG:
	# 	logger.debug('recheck: ' + str((max_score, maximal_clique)))
	# 	# print('recheck')
	# 	# print(max_score, maximal_clique)

	s2t, t2s, ind1_list, ind2_list = generate_index_mapping(profile1, profile2, maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, nucl_scoring_data, dBG)

	merged_profile = generate_merged_profile_based_on_alignment(profile1, profile2, ind1_list, ind2_list, maximal_clique, interaction_ind_pair_list1, interaction_ind_pair_list2, scores, nucl_scoring_data, dBG)
	# print('got merged_profile')

	merged_profile = revise_profile_based_on_gap(merged_profile, frequent_gap_filtering_threshold)
	# print('got revised merged_profile')

	new_profiles = profiles[:ind1] + profiles[ind1+1:ind2] + profiles[ind2+1:]
	new_profiles.append(merged_profile)
	# print('new loopcounts')
	# for profile in new_profiles:
	# 	print(profile.loopCount)

	return generate_single_profile(new_profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)

def find_best_aligned_profile_pair(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG):
	# print_formatted_profile(profiles[0])
	total = len(profiles)
	maximal_clique = []
	max_score = -Infinity
	max_ind1 = -1
	max_ind2 = -1

	max_profile1 = None
	max_profile2_o = None

	# return 5.0, [], 0, 1

	for i in range(total-1):
		profile1 = profiles[i]
		if dBG:
			logger.info('Profile1:')
			print_formatted_profile(profile1)
		for j in range(i+1, total):
			# print('checking index ' + str(i) + ' and ' + str(j))
			profile2 = profiles[j]
			if dBG:
				logger.info('Profile1:')
				print_formatted_profile(profile1)

			score, clique, profile2_o = find_best_aligned_data(profile1, generate_all_profile_combinations(profile2, dBG), scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)

			if dBG:
				logger.info('Best orientation:')
				print_formatted_profile(profile2_o)

			if score > max_score:
				max_score = score
				maximal_clique = copy.deepcopy(clique)
				max_ind1 = i
				max_ind2 = j
				max_profile1 = profile1
				max_profile2_o = profile2_o

	return max_score, maximal_clique, max_ind1, max_ind2, max_profile1, max_profile2_o

def find_best_aligned_data(profile1, profile2s, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG):
	max_max_score = -Infinity
	max_maximal_clique = []
	max_profile2_orientation = None

	for profile2 in profile2s:
		interaction_ind_pair_list1 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile1))
		interaction_ind_pair_list2 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile2))

		# interaction_compatibility_matrix = compute_interaction_compatibility(interaction_ind_pair_list1, interaction_ind_pair_list2)
		# interaction_compatibility_matrix, rank_dict = compute_interaction_compatibility(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, include_stackings)
		interaction_compatibility_matrix = compute_interaction_compatibility(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, include_stackings)
		# print(interaction_compatibility_matrix)

		profile1_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile1, scores.penalties)
		profile2_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile2, scores.penalties)

		# print('calling from 192')
		# max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, rank_dict, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
		max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
		if dBG:
			print('max_score', 'maximal_clique')
			print(max_score, maximal_clique)
			print_formatted_profile(profile2)
		if max_score > max_max_score:
			max_max_score = max_score
			max_maximal_clique = copy.deepcopy(maximal_clique)
			max_profile2_orientation = profile2

	return max_max_score, max_maximal_clique, max_profile2_orientation

def get_scaled_score(score):
	return round(score / 10000.0, 2)