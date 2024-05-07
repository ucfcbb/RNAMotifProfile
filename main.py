import sys
import os
import time
import argparse
import logging
import collections
import multiprocessing as mp

from classes import *
from config import *
from logger import *
from prepare import *
from converter import *
from profile import *
from utils import *
from input_utils import *
from interaction_utils import *
from scoring_utils import *
from search import *

from benchmark_helper import *

def main():

	process_start_time = time.time()

	parser = argparse.ArgumentParser(description='Prepare input for RNAMotifComp')

	# Input/Output
	parser.add_argument('-i', nargs='?', default='input/sample1.in', const='input/sample1.in', help='Input file containing motifs.')
	parser.add_argument('-o', nargs='?', default='', const='', help='Subdirectory inside Output directory to save the results.')

	parser.add_argument('-f1', nargs='?', default='0.0', const='0.1', help='Filter less occurring interactions and nucleotides (value range 0.0-1.0).')
	parser.add_argument('-f2', nargs='?', default='100.0', const='75.0', help='Filter more occuring gaps that are not associated with any interactions (value range 100.0-0.0).')
	parser.add_argument('-bnb', nargs='?', default=False, const=True, help='Use branch and bound to get results faster.')	

	# Searching mode
	parser.add_argument('-m', nargs='?', default='', const='', help='Mode selector. Provide "search" to search motifs.')
	parser.add_argument('-p', nargs='?', default='', const='', help='Profile (.pfl) file selector. Required, if "search" mode is selected.')
	parser.add_argument('-c', nargs='?', default='', const='', help='Chain selector (<PDB_ID>_<Chain_ID>). Required, if "search" mode is selected.')
	parser.add_argument('-z', nargs='?', default=False, const=True, help='Calculate z-score for each query motif using mean and std of the used profile with the motifs used to generate this profile.')
	# parser.add_argument('-x', nargs='?', default=False, const=True, help='Include ScanX alignment score (if particular structure file provided) in search mode.')
	# parser.add_argument('-xs', nargs='?', default='', const='', help='Include ScanX structure (.struct) to generate alignment. Required, if -x is selected.')


	# Debug
	parser.add_argument('-d', nargs='?', default=False, const=True, help='Debug mode selector.')

	# Bonus and penalties
	parser.add_argument('-g', nargs='?', default='-5.0', const='-5.0', help='Gap opening penalty.')
	parser.add_argument('-e', nargs='?', default='-3.0', const='-3.0', help='Gap extension penalty.')
	parser.add_argument('-b', nargs='?', default='-3.0', const='-3.0', help='Missing BP penalty.')
	parser.add_argument('-s', nargs='?', default='-2.0', const='-2.0', help='Missing stack penalty.')
	parser.add_argument('-y', nargs='?', default='2.0', const='2.0', help='H-bond match bonus.')
	parser.add_argument('-t', nargs='?', default='3.0', const='3.0', help='Triple-interaction bonus.')

	try:
		args = parser.parse_args()
	except Exception as e:
		parser.print_help()
		sys.exit()

	validate_params(args)

	dBG = args.d

	gap_opening_penalty = float(args.g)
	gap_extension_penalty = float(args.e)
	missing_bp_penalty = float(args.b)
	missing_stk_penalty = float(args.s)

	hbond_match_bonus = float(args.y)
	triple_interaction_bonus = float(args.t)

	input_fname = args.i
	output_subdir_name = args.o

	mode = str(args.m)
	target_profile_fname = str(args.p)
	calculate_zscore = args.z
	# include_scanx_score = args.x
	# scanx_structure_file = str(args.xs)
	include_scanx_score = False
	scanx_structure_file = None
	pdb_chains_to_search = args.c

	less_freq_nucl_and_intera_filter_threshold = float(args.f1)
	frequent_gap_filtering_threshold = float(args.f2)
	use_branch_and_bound = args.bnb

	# dBG = False	

	# gap_opening_penalty = -5.0
	# gap_extension_penalty = -3.0
	# missing_bp_penalty = -3.0
	# missing_stk_penalty = -2.0
	penalties = Penalties(gap_opening_penalty, gap_extension_penalty, missing_bp_penalty, missing_stk_penalty)

	# hbond_match_bonus = 2.0
	# triple_interaction_bonus = 3.0
	bonuses = Bonuses(hbond_match_bonus, triple_interaction_bonus)

	# hbond_match_base = 2.0
	# weight_isosteric = 2.0
	# weight_nonadjacent_stacking = 1.0
	# weight_adjacent_stacking = 0.2
	# weight_sequence = 0.1
	
	# asym_nuc = -0.2
	# asym_loop = -2.0
	# cons_stacking = 0.5
	# stack_to_bulge = -1.5
	# stack_to_internal_asym = -1.0
	# stack_to_internal_sym = -0.5

	# scores = Scores(penalties, bonuses, hbond_match_base, weight_isosteric, weight_nonadjacent_stacking, weight_adjacent_stacking, weight_sequence, asym_nuc, asym_loop, cons_stacking, stack_to_bulge, stack_to_internal_asym, stack_to_internal_sym)
	scores = get_scores(penalties, bonuses)

	# print(scores.weight_sequence)
	# sys.exit()

	# isosteric_file_name = 'data/mat/iso.mat'
	# stacking_file_name = 'data/mat/stk.mat'
	# nucleotide_file_name = 'data/mat/nuc.mat'
	isosteric_file_name, stacking_file_name, nucleotide_file_name = get_mat_file_names()

	initialized_basepair_index = get_initialized_basepair_index()
	bp_scoring_data = read_isosteric_scoring_data(isosteric_file_name)
	stk_scoring_data = read_stacking_scoring_data(stacking_file_name)
	nucl_scoring_data = read_nucleotide_scoring_data(nucleotide_file_name)
	# compute_edge_substitution()


	# loop_dir = 'data/loops/'
	# output_dir = 'output/'

	directories = get_base_dir_names()
	create_required_directories(directories)
	input_index_type, annotation_source, items_per_chunk, mp_number_of_process, content_download_attempts = get_misc_params()

	include_stackings = get_include_stackings_flag(mode)

	# input_fname = 'example_cluster_KT1.csv'
	# input_fname = 'KT.csv'
	# input_fname = 'KT-sub1.csv'
	# input_fname = 'KT-sub2.csv'
	# input_fname = 'KT-sub3.csv'
	# input_fname = 'KT-sub4.csv'

	# if load_search_input_from_chain == False:
	fp_input = open(input_fname)
	loop_list = csv_to_list(fp_input.readlines())
	fp_input.close()

	families = {}
	loop_count = 0
	for item in loop_list:
		if len(item) > 1:
			# families[item[0]] = list(map(lambda x: str(strToNode(x)), item[1:])) # item[1:]
			families[item[0]] = item[1:]

	print('\n')
	if mode != 'search':
		# profile generation

		# if load_search_input_from_chain == True:
		# 	fp_input = open(input_fname)
		# 	loop_list = csv_to_list(fp_input.readlines())
		# 	fp_input.close()

		# 	families = {}
		# 	loop_count = 0
		# 	for item in loop_list:
		# 		if len(item) > 1:
		# 			# families[item[0]] = list(map(lambda x: str(strToNode(x)), item[1:])) # item[1:]
		# 			families[item[0]] = item[1:]

		t = len(families)
		if t == 0:
			logger.info('No families found to generate profile.')
		elif t == 1:
			logger.info(str(t) + ' family found to generate profile.')
		else:
			logger.info(str(t) + ' different families found to generate profiles.')

		print('')

		prepare_data(families, directories, annotation_source, content_download_attempts, mp_number_of_process)
		if input_index_type == 'pdb':
			families = convert_a_cluster_from_PDB_to_FASTA(families, directories)
			for family_id in families:
				families[family_id] = list(map(lambda x: str(strToNode(x)), families[family_id]))

		# else:
		# 	families = convert_a_cluster_from_FASTA_to_PDB(families, directories)
		# 	# for family_id in families:
		# 	# 	families[family_id] = list(map(lambda x: str(strToNode(x)), families[family_id]))

		# 	bs_name = os.path.basename(input_fname)
		# 	pcs = bs_name.strip().split('.')
		# 	fn = '.'.join(pcs[:-1]) + '_PDB.' + pcs[-1]
		# 	fp = open(fn, 'w')
		# 	for fam_id in families:
		# 		fp.write(fam_id+',')
		# 		fp.write(','.join(families[fam_id]))
		# 		fp.write('\n')
		# 	fp.close()
		# 	sys.exit()

		# print_input_stat(families)
		# print_loop_intera_stats(families, directories)
		# sys.exit()


		loop_count = 0
		loop_node_list_str = []
		for family_id in families:
			loops = families[family_id]
			loop_count += len(loops)
			for loop in loops:
				loop = str(strToNode(loop))
				loop_node_list_str.append(loop)
		
		duplicates = [item for item, count in collections.Counter(loop_node_list_str).items() if count > 1]
		if len(duplicates) > 0:
			print('duplicates:')
			print(duplicates)
			# sys.exit()

		loop_node_list_str = sorted(list(set(loop_node_list_str)))
		# loop_node_list_str = sorted(loop_node_list_str)

		logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in ' + str(len(families)) + ' famil' + ('ies' if len(families) > 1 else 'y') + '.')
		print('')

		prepare_loop_files(loop_node_list_str, directories, annotation_source, mp_number_of_process, get_env())	#chkd
		# print('Done till here')
		# sys.exit()

		output_dir = directories.output_dir
		if len(output_subdir_name) > 0:
			output_dir = os.path.join(output_dir, output_subdir_name)
		create_directory(output_dir)

		if not os.path.exists('stat_data.csv'):
			initiate_stat_file()
		for family_id in families:
			sub_process_start_time = time.time()
			logger.info('Processing ' + str(family_id) + ' to generate profile.')
			logger.info('- Gathering loop data')

			loops = families[family_id]

			way_counts = get_way_counts(loops)
			if len(way_counts) > 1:
				logger.error('Input contains motifs with varying segments. Skipping this family/cluster.')
				continue

			motifs = []
			not_found_list = []
			for loop in loops:
				file_name = os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf')
				if not os.path.exists(file_name):
					not_found_list.append(loop)
					continue
				loop_data = load_loop_data(os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf'), False)
				# print(loop_data[1])
				# sys.exit()
				motif = Motif(loop_data[0], loop_data[1], loop_data[3], loop_data[5], loop)
				# print(motif)
				motifs.append(motif)
			
			if len(not_found_list) > 0:
				# print(str(len(not_found_list)) + ' loop files not found, skipping these loops in the profile generation.')
				logger.warning(str(len(not_found_list)) + ' loop files not found, skipping these loops in the profile generation.')
			# sys.exit()
			profiles = list(map(lambda x: MotifToProfile(x), motifs))

			# print(profiles[0].sequenceP)
			# print_formatted_profile(profiles[0])
			# scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index

			# process_start_time = time.time()
			logger.info('- Aligning loops to generate profile')
			profile = generate_single_profile(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)
			# profile = generate_single_profile_multi_process(profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, items_per_chunk, mp_number_of_process, include_stackings, frequent_gap_filtering_threshold, use_branch_and_bound, dBG)
			sub_process_end_time = time.time()
			logger.info('Profile generation complete for ' + str(family_id) + '.')
			print('Time taken: ' + str(round((sub_process_end_time - sub_process_start_time), 3)) + ' seconds.')

			if dBG:
				logger.debug('Single profile for family_id: ' + family_id)
				# print('Single profile for family_id: ' + family_id)
				print_formatted_profile(profile)

			time_str = time.strftime('%Y%m%d-%H%M%S')
			output_fname = family_id + '_' + time_str + '.pfl'
			
			write_formatted_profile_to_file(profile, output_dir, output_fname)

			# if way_counts[0] > 4:
			# 	logger.error('Skipping stat generation for profile with loops haveing 5 or more way junctions.')
			# 	continue
			if family_id == 'ML31_1':
				logger.error('Skipping stat generation for ML31_1 profile as an exception.')
				continue
			generate_stats_for_aligning_a_profile_with_source_loops(family_id, profile, profiles, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG, output_fname, less_freq_nucl_and_intera_filter_threshold, frequent_gap_filtering_threshold)
			
		process_end_time = time.time()
		print('Total time taken: ' + str(round((process_end_time - process_start_time), 3)) + ' seconds.\n')

		logger.info('Please check the \'' + output_dir[len(os.path.dirname(directories.base_dir))+1:] + '\' directory for generated profiles.\n')

	else:
		#profile based search

		profile_fname = target_profile_fname

		if load_search_input_from_chain == True:
			pdb_chains_list_to_search = pdb_chains_to_search.strip().split(',')
			output_fname = os.path.basename(profile_fname).strip().split('_')[0] + '_vs_' + '-'.join(pdb_chains_list_to_search) + '.out'
		else:
			query_loops_fname = os.path.basename(input_fname)
			output_fname = os.path.basename(profile_fname).strip().split('_')[0] + '_vs_' + query_loops_fname.strip().split('.')[0] + '.out'

		output_dir = directories.output_dir
		if len(output_subdir_name) > 0:
			output_dir = os.path.join(output_dir, output_subdir_name)
		create_directory(output_dir)
		
		output_fname = os.path.join(output_dir, output_fname)

		fq = open(output_fname, 'w')

		logger.info('Loading profile data.')
		target_profile = load_profile_data(profile_fname, os.path.basename(profile_fname).strip().split('_')[0])

		target_profile = filter_less_significant_nucl_and_interactions_from_a_profile(target_profile, less_freq_nucl_and_intera_filter_threshold)
		target_profile = revise_profile_based_on_gap(target_profile, frequent_gap_filtering_threshold)

		loop_node_list_str = []
		if load_search_input_from_chain == True:
			junction_count = len(target_profile.break_points)

			pdb_chains_list_to_search = pdb_chains_to_search.strip().split(',')
			loop_count = 0
			loops_list = []
			for pdb_chain_to_search in pdb_chains_list_to_search:
				prepare_data_for_a_chain(pdb_chain_to_search, directories, annotation_source, content_download_attempts, mp_number_of_process)
				loops = get_loops_from_a_chain(pdb_chain_to_search, directories, annotation_source, junction_count)
				loop_count += len(loops)
				loops_list.extend(loops)
			
			loop_node_list_str = list(map(lambda x: str(strToNode(x)), loops_list))
			loop_node_list_str = sorted(list(set(loop_node_list_str)))
			# print(loop_node_list_str)
			# print(len(loop_node_list_str))
			# sys.exit()
			logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in the provided chain(s) for searching.')
			print('')

			prepare_loop_files(loop_node_list_str, directories, annotation_source, mp_number_of_process, get_env())	#chkd
		else:
			prepare_data(families, directories, annotation_source, content_download_attempts, mp_number_of_process)
			if input_index_type == 'pdb':
				families = convert_a_cluster_from_PDB_to_FASTA(families, directories)
				for family_id in families:
					families[family_id] = list(map(lambda x: str(strToNode(x)), families[family_id]))


			loop_count = 0
			loop_node_list_str = []
			for family_id in families:
				loops = families[family_id]
				loop_count += len(loops)
				for loop in loops:
					loop = str(strToNode(loop))
					loop_node_list_str.append(loop)
			
			duplicates = [item for item, count in collections.Counter(loop_node_list_str).items() if count > 1]
			if len(duplicates) > 0:
				print('duplicates:')
				print(duplicates)

			loop_node_list_str = sorted(list(set(loop_node_list_str)))

			logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in ' + str(len(families)) + ' famil' + ('ies' if len(families) > 1 else 'y') + ' for searching using provided profile.')
			print('')

			prepare_loop_files(loop_node_list_str, directories, annotation_source, mp_number_of_process, get_env())	#chkd


		# sys.exit()

		query_loop_list = copy.deepcopy(loop_node_list_str)
		# query_loop_list = loop_list[0][1:]
		# if input_index_type == 'pdb':
		# 	query_loop_list = list(map(lambda x: convert_a_loop_from_PDB_to_FASTA(x, directories), query_loop_list))
		# convert_a_loop_from_PDB_to_FASTA(loop, directories)

		# fq.write(os.path.basename(profile_fname).strip().split('_')[0] + ' profile vs ' + query_loops_fname.strip().split('.')[0] + ' loops:\n')
		scanx_column_titles = ''
		if include_scanx_score:
			# scanx_column_titles = 'Scanx aln score' + '\t' + 'Scanx matching bp count' + '\t' + 'Scanx rank' + '\t' + 'Profile rank'
			scanx_column_titles = 'Scanx aln score' + '\t' + 'Scanx matching bp count'

		if dBG:
			fq.write('Motif' + '\t' + 'Profile aln score' + '\t' + 'Matching interactions' + '\t' + 'Matching bp count')
			if calculate_zscore:
				fq.write('\t' + 'Z-score')
			fq.write('\t' + 'Family name' + '\t' + 'Subfamily name' + '\t' + scanx_column_titles + '\n')
		else:
			fq.write('Motif' + '\t' + 'Profile aln score' + '\t' + 'Matching bp count')
			if calculate_zscore:
				fq.write('\t' + 'Z-score')
			fq.write('\t' + scanx_column_titles + '\n')

		fq.close()
		# print(target_profile.loopName)
		# print(query_loop_list)
		# sys.exit()

		search_for_motifs_using_a_profile(os.path.basename(profile_fname), target_profile, query_loop_list, directories, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG, output_fname, calculate_zscore, include_scanx_score, scanx_structure_file)
		logger.info('\nSearch complete.\nA profile-to-profile alignment score is generated for each input motif.\nMotifs with higher score are more likely to be of the same family of the provided profile.\n')

		logger.info('Please check the \'' + output_dir[len(os.path.dirname(directories.base_dir))+1:] + '\' directory for the search results.\n')

if __name__ == "__main__":
	main()
