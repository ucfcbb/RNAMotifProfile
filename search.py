import os
import sys
import time
import collections
import copy
import platform

from classes import *
from config import *
from logger import *
from utils import *
from prepare import *
from input_utils import *
from interaction_utils import *
from scoring_utils import *
from converter import *
from profile import *
from benchmark_helper import *

def parse_scanx_alignment_block_raw(lines, line_index):
    # r1 = lines[line_index].split('::')[1].split(' and ')[0].strip().strip(':')
    # r2 = lines[line_index].split('::')[1].split(' and ')[1].strip().strip(':')
    r1 = ''
    r2 = ''

    score_text = lines[line_index+1].split(':')[1].strip()
    if score_text == '':
        # score = -50.
        logger.error('ERROR: No alignment score found for: ' + r1 + ' and ' + r2)
        sys.exit()
    else:
        score = float(score_text)

    cr1 = lines[line_index+3].split(':')[1].strip()
    cr2 = lines[line_index+4].split(':')[1].strip()

    aln1 = lines[line_index+6].strip()
    aln2 = lines[line_index+7].strip()

    matching_bp_info = []
    matching_stk_info = []
    i = line_index+10
    while not lines[i].startswith('#  Matched base-stacking interactions: '):
        matching_bp_info.append(list(map(lambda x: x.strip(), lines[i].strip().split('MATCHES'))))
        i += 1

    i += 1
    # while not lines[i].startswith('Total Elapsed Time :'):
    while i < len(lines):
        matching_stk_info.append(list(map(lambda x: x.strip(), lines[i].strip().split('MATCHES'))))
        i += 1

    is_copied = False
    # if lines[i].strip().endswith('(copied)'):
    #     is_copied = True
    # elapsed_time = lines[i].strip().split(':')[1].strip().split(' ')[0].strip()
    elapsed_time = 0

    return (r1, r2, cr1, cr2, aln1, aln2, score, matching_bp_info, matching_stk_info, elapsed_time, is_copied, i)

def rotate(l, x):
    return l[-x:] + l[:-x]

def get_all_loop_combination(loop):
    loop_combinations = []
    pdb_chain, regions = loop.strip().split(':')
    regions = regions.strip().split('_')
    loop = []
    for region in regions:
        s, e = region.strip().split('-')
        loop.append((s, e))
    
    for i in range(len(loop)):
        loop_i = rotate(loop, i)
        loop_combinations.append(pdb_chain + ':' + '_'.join(list(map(lambda x: '-'.join(x), loop_i))))

    # print(loop_combinations)
    return loop_combinations

def generate_scanx_data(directories, loopName, scanx_structure_file):
	ScanX_dir = get_scanx_dir()
	scanx_aln_executable = os.path.join(ScanX_dir, 'bin/align')
	if platform.system() == 'Darwin':
		scanx_aln_executable = os.path.join(ScanX_dir, 'bin/align.mac')
	alignment_dir = directories.temp_dir
	output_fn = os.path.join(alignment_dir, 'alignment_result.aln')
	os.environ['RNAMOTIFSCANX_PATH'] = ScanX_dir
	os.environ['RNAVIEW'] = os.path.join(ScanX_dir, 'thirdparty/RNAVIEW')

	# file1 = '/home/mahfuz/Downloads/RNAMotifScanX-release_v0.0.5_x86-64_rhel/RNAMotifScanX-release/models/k-turn_consensus.struct'
	file1 = scanx_structure_file
	# file1 = '/home/mahfuz/Downloads/RNAMotifScanX-release_v0.0.5_x86-64_rhel/RNAMotifScanX-release/models/sarcin-ricin_consensus.struct'

	if os.path.exists(output_fn):
		os.remove(output_fn)

	max_alignment_data = None
	for loop in get_all_loop_combination(loopName):
		# print('generating alignment with ' + loop)
		logger.info('Generating alignment with ' + loop)
		file2 = os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf')
		# print('%s %s %s >> %s' % (scanx_aln_executable, file1, file2, output_fn))
		os.system('%s %s %s >> %s' % (scanx_aln_executable, file1, file2, output_fn))
		time.sleep(3)
		# print('parsing')
		fp = open(output_fn)
		lines = fp.readlines()
		fp.close()
		line_index = 0
		alignment_data = None
		while line_index < len(lines):
			if lines[line_index].startswith('#  Aligning'):
				alignment_data = parse_scanx_alignment_block_raw(lines, line_index)
				break
			line_index += 1
		# if alignment_data == None:
		# 	print('why?')
		# 	sys.exit()
		if max_alignment_data == None:
			max_alignment_data = alignment_data
		elif max_alignment_data[6] < alignment_data[6]:
			max_alignment_data = alignment_data
		os.remove(output_fn)
		time.sleep(0.5)
		# print('removed')

	r1, r2, cr1, cr2, aln1, aln2, score, matching_bp_info, matching_stk_info, elapsed_time, is_copied, line_index = max_alignment_data
	return score, len(matching_bp_info)		

def search_for_motifs_using_a_profile(profile_fname, target_profile, loop_list, directories, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG, output_fname, calculate_zscore, include_scanx_score, scanx_structure_file):
	# profile1 = target_profile
	# profile2 = target_profile
	# interaction_ind_pair_list1 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile1))
	# interaction_ind_pair_list2 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(profile2))
	# interaction_compatibility_matrix = compute_interaction_compatibility(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, include_stackings)
	# # print(interaction_compatibility_matrix)

	# profile1_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile1, scores.penalties)
	# profile2_nucleotide_deletion_penalties = compute_nucleotide_deletion_penalties(profile2, scores.penalties)

	# max_score, maximal_clique = find_maximal_optimal_clique(profile1, profile2, interaction_ind_pair_list1, interaction_ind_pair_list2, interaction_compatibility_matrix, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, profile1_nucleotide_deletion_penalties, profile2_nucleotide_deletion_penalties, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
	# # max_score = round(max_score / 10000.0, 2)
	# max_score = get_scaled_score(max_score)
	# print('self alignment score')
	# print(max_score)

	# calculate_zscore = True

	mean_std_data = load_stat_data(profile_fname)
	mean, std = mean_std_data
	if mean == None and std == None:
		logger.warning('Stat data not found for the target profile. No prediction will be provided.')
		calculate_zscore = False

	logger.info('Loading motif information.')
	
	motifs = []
	for loop in loop_list:
		loop_data = load_loop_data(os.path.join(directories.loop_dir, loop.replace(':', '_') + '.smf'), False)
		# print(loop_data[1])
		# sys.exit()
		motif = Motif(loop_data[0], loop_data[1], loop_data[3], loop_data[5], loop)
		# print(motif)
		motifs.append(motif)

	logger.info('Converting all motifs into profiles.')
	query_profiles = list(map(lambda x: MotifToProfile(x), motifs))

	interaction_ind_pair_list1 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(target_profile))
	interaction_count1 = len(interaction_ind_pair_list1)

	result_list = []

	logger.info('Performing profile-to-profile alignment of the query motifs with target profile.')
	skipped_count = 0
	for q_profile in query_profiles:
		matched_bp_cnt = 0
		matched_stk_cnt = 0

		# print('here')
		try:
			score, clique, q_profile_o = find_best_aligned_data(target_profile, generate_all_profile_combinations(q_profile, dBG), scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
			# score = round(score / 10000.0, 2)
			score = get_scaled_score(score)
		except Exception as e:
			# raise e
			skipped_count += 1
			continue
		# score, clique, q_profile_o = find_best_aligned_data(target_profile, generate_all_profile_combinations(q_profile, dBG), scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG)
		# print('here1')
		# print(score)
		# print(clique)

		
		interaction_ind_pair_list2 = sort_interaction_ind_pairs(get_interaction_ind_pair_list(q_profile_o))


		
		interaction_count2 = len(interaction_ind_pair_list2)
		clique_size = len(clique)

		for i in range(clique_size):
			#  the two base interactions idx_i and idx_j match
			ind_i = clique[i] // interaction_count2
			ind_j = clique[i] % interaction_count2

			interaction_ind_pair_list1[ind_i][2]
			if interaction_ind_pair_list2[ind_j][2] == 'b':
				matched_bp_cnt += 1
			elif interaction_ind_pair_list2[ind_j][2] == 's':
				matched_stk_cnt += 1

		# print('matched_bp_cnt: ' + str(matched_bp_cnt))
		# print('matched_stk_cnt: ' + str(matched_stk_cnt))

		if calculate_zscore:
			zscore = get_zscore(score, mean, std)
			result_list.append((q_profile_o.loopName, score, clique, matched_bp_cnt, matched_stk_cnt, zscore))
		else:
			result_list.append((q_profile_o.loopName, score, clique, matched_bp_cnt, matched_stk_cnt))

	# result_list = sorted(result_list, key = lambda x: (x[3], x[1]), reverse = True)
	result_list = sorted(result_list, key = lambda x: (x[1], x[3]), reverse = True)

	# print('Skipped: ' + str(skipped_count))
	if skipped_count > 0:
		logger.warning('Skipped ' + str(skipped_count) + ' motifs while searching due to some technical issues in their structure(s).')

	# output_fname = os.path.basename(profile_fname).strip().split('_')[0] + '_vs_' + query_loops_fname.strip().split('.')[0] + '.out'
	
	# fq = open(output_fname, 'a')

	sub_families = {}
	if dBG:
		# fr = open('subfamily_cluster_cluster.csv')
		rnamotifcontrast_data = 'subfamily_cluster_contrast.csv'
		rnamotifcontrast_data = 'subfamily_cluster_cluster.csv'
		fr = open(rnamotifcontrast_data)
		loop_list = csv_to_list(fr.readlines())
		fr.close()

		for item in loop_list:
			if len(item) > 1:
				# families[item[0]] = list(map(lambda x: str(strToNode(x)), item[1:])) # item[1:]
				sub_families[item[0]] = item[1:]

		for sf in sub_families:
			sub_families[sf] = map(lambda x: str(strToNode(x)), sub_families[sf])

		# print(sub_families)

	for result in result_list:
		if len(result) == 5:
			loopName, score, clique, matched_bp_cnt, matched_stk_cnt = result
		else:
			loopName, score, clique, matched_bp_cnt, matched_stk_cnt, zscore = result

		fq = open(output_fname, 'a')
		loopName_PDB = convert_a_loop_from_FASTA_to_PDB(loopName, directories)
		############################## TEMPORARY #############################
		# loopName_PDB = loopName
		############################## TEMPORARY #############################
		if dBG:
			# fq.write(str(loopName_PDB) + '\t' + str(score) + '\t' + str(clique) + '\t' + str(matched_bp_cnt) + '\t' + str(matched_stk_cnt) + '\t')
			fq.write(str(loopName_PDB) + '\t' + str(score) + '\t' + str(clique) + '\t' + str(matched_bp_cnt))
			if calculate_zscore:
				fq.write('\t' + str(round(zscore, 2)))

			family_name = ''
			subfamily_name = ''
			for sf in sub_families:
				if loopName_PDB in sub_families[sf]:
					family_name = '-'.join(sf.strip().split('-')[:-1])
					subfamily_name = sf.strip().split('-')[-1]
					break
			fq.write('\t' + str(family_name) + '\t' + str(subfamily_name))

		else:
			fq.write(str(loopName_PDB) + '\t' + str(score) + '\t' + str(matched_bp_cnt))

			if calculate_zscore:
				fq.write('\t' + str(round(zscore, 2)))

		if include_scanx_score:
			logger.info('Generating ScanX alignment data for ' + str(loopName))
			scanx_aln_score, num_bp_match = generate_scanx_data(directories, loopName, scanx_structure_file)
			fq.write('\t' + str(scanx_aln_score) + '\t' + str(num_bp_match))
		
		fq.write('\n')

		fq.close()

# ### TEMPORARY CODE TO TEST ###
# def main():

# 	dBG = False
# 	use_branch_and_bound = True

# 	gap_opening_penalty = -5.0
# 	gap_extension_penalty = -3.0
# 	missing_bp_penalty = -3.0
# 	missing_stk_penalty = -2.0
# 	penalties = Penalties(gap_opening_penalty, gap_extension_penalty, missing_bp_penalty, missing_stk_penalty)

# 	hbond_match_bonus = 2.0
# 	triple_interaction_bonus = 3.0
# 	bonuses = Bonuses(hbond_match_bonus, triple_interaction_bonus)

# 	# hbond_match_base = 2.0
# 	# weight_isosteric = 2.0
# 	# weight_nonadjacent_stacking = 1.0
# 	# weight_adjacent_stacking = 0.2
# 	# weight_sequence = 0.1
	
# 	# asym_nuc = -0.2
# 	# asym_loop = -2.0
# 	# cons_stacking = 0.5
# 	# stack_to_bulge = -1.5
# 	# stack_to_internal_asym = -1.0
# 	# stack_to_internal_sym = -0.5

# 	# scores = Scores(penalties, bonuses, hbond_match_base, weight_isosteric, weight_nonadjacent_stacking, weight_adjacent_stacking, weight_sequence, asym_nuc, asym_loop, cons_stacking, stack_to_bulge, stack_to_internal_asym, stack_to_internal_sym)
# 	scores = get_scores(penalties, bonuses)

# 	# print(scores.weight_sequence)
# 	# sys.exit()

# 	# isosteric_file_name = 'data/mat/iso.mat'
# 	# stacking_file_name = 'data/mat/stk.mat'
# 	# nucleotide_file_name = 'data/mat/nuc.mat'
# 	isosteric_file_name, stacking_file_name, nucleotide_file_name = get_mat_file_names()

# 	include_stackings = get_include_stackings_flag('search')

# 	initialized_basepair_index = get_initialized_basepair_index()
# 	bp_scoring_data = read_isosteric_scoring_data(isosteric_file_name)
# 	stk_scoring_data = read_stacking_scoring_data(stacking_file_name)
# 	nucl_scoring_data = read_nucleotide_scoring_data(nucleotide_file_name)
# 	# compute_edge_substitution()


# 	# loop_dir = 'data/loops/'
# 	# output_dir = 'output/'

# 	directories = get_base_dir_names()
# 	create_required_directories(directories)
# 	input_index_type, annotation_source, items_per_chunk, mp_number_of_process, content_download_attempts = get_misc_params()


# 	# profile_fname = 'output/Kink-turn-Sub1_20230808-123909.pfl'
# 	profile_fname = 'output/Stable2/Kink-turn_20230808-114533.pfl'
# 	# profile_fname = 'output/Kink-turn-Sub4_20230727-125146.pfl'
# 	# profile_fname = 'output/Stable2/Sarcin-ricin_20230808-123344.pfl'
# 	# profile_fname = 'output/C-loop_20230808-113357.pfl'
# 	# profile_fname = 'output/E-loop_20230808-113413.pfl'
# 	# profile_fname = 'output/reverse-Kink-turn_20230808-114536.pfl'
# 	# profile_fname = 'output/Hook-turn_20230808-113414.pfl'
# 	# query_loops_fname = 'KT-sub1.csv'
# 	# query_loops_fname = 'KT-sub4.csv'
# 	query_loops_fname = '4V9F_loops.csv'	#pdb
# 	# query_loops_fname = '5J7L_loops.csv'
# 	# query_loops_fname = '4V9F_loops_from_IL_XL_partial.csv'
# 	# query_loops_fname = '4V9F_loops_from_IL_XL_fasta.csv'
# 	# query_loops_fname = '5J7L_loops_from_IL_XL_fasta.csv'
# 	# query_loops_fname = '7KGB_A_loops.csv'
# 	# query_loops_fname = 'IL_cluster_input_exclusive_fasta.csv'
# 	# query_loops_fname = '6ERI_AA_excl_loops_fasta.csv'
# 	# query_loops_fname = '5O60_A_excl_loops_fasta.csv'
# 	# query_loops_fname = '5V7Q_A_excl_loops_fasta.csv'
# 	# query_loops_fname = '6HA1_A_excl_loops_fasta.csv'
# 	# query_loops_fname = '6HA1_A_a_excl_loops_fasta.csv'
# 	# query_loops_fname = '6D9J_5_excl_loops_fasta.csv'
# 	# query_loops_fname = 'EL-KT-SR_from_cluster_result_exclusive_fasta.csv'
# 	# query_loops_fname = 'partial_excl_loops_fasta.csv'
# 	# query_loops_fname = '6ERI_excl_loops_fasta.csv'
	

# 	fp_input = open(query_loops_fname)

# 	output_fname = os.path.basename(profile_fname).strip().split('_')[0] + '_vs_' + query_loops_fname.strip().split('.')[0] + '.out'
# 	fq = open(output_fname, 'w')

# 	target_profile = load_profile_data(profile_fname, os.path.basename(profile_fname).strip().split('_')[0])

# 	target_profile = filter_less_significant_nucl_and_interactions_from_a_profile(target_profile)
# 	target_profile = revise_profile_based_on_gap(target_profile, 80.0)
# 	# sys.exit()
# 	loop_list = csv_to_list(fp_input.readlines())
	
# 	families = {}
# 	families[loop_list[0][0]] = loop_list[0][1:]
# 	# print(families)
# 	# input_index_type, annotation_source, items_per_chunk, mp_number_of_process, content_download_attempts = get_misc_params()
# 	# input_index_type = 'fasta'
# 	prepare_data(families, directories, annotation_source, content_download_attempts, mp_number_of_process)
# 	if input_index_type == 'pdb':
# 		families = convert_a_cluster_from_PDB_to_FASTA(families, directories)
# 		for family_id in families:
# 			families[family_id] = map(lambda x: str(strToNode(x)), families[family_id])


# 	loop_count = 0
# 	loop_node_list_str = []
# 	for family_id in families:
# 		loops = families[family_id]
# 		loop_count += len(loops)
# 		for loop in loops:
# 			loop = str(strToNode(loop))
# 			loop_node_list_str.append(loop)
	
# 	duplicates = [item for item, count in collections.Counter(loop_node_list_str).items() if count > 1]
# 	if len(duplicates) > 0:
# 		print('duplicates:')
# 		print(duplicates)

# 	loop_node_list_str = sorted(list(set(loop_node_list_str)))

# 	logger.info(str(loop_count) + ' loops (' + str(len(loop_node_list_str)) + ' unique) found in ' + str(len(families)) + ' famil' + ('ies' if len(families) > 1 else 'y') + '.')
# 	print('')

# 	prepare_loop_files(loop_node_list_str, directories, annotation_source, mp_number_of_process, get_env())	#chkd


# 	# sys.exit()

# 	query_loop_list = copy.deepcopy(loop_node_list_str)
# 	# query_loop_list = loop_list[0][1:]
# 	# if input_index_type == 'pdb':
# 	# 	query_loop_list = list(map(lambda x: convert_a_loop_from_PDB_to_FASTA(x, directories), query_loop_list))
# 	# convert_a_loop_from_PDB_to_FASTA(loop, directories)

# 	fq.write(os.path.basename(profile_fname).strip().split('_')[0] + ' profile vs ' + query_loops_fname.strip().split('.')[0] + ' loops:\n')
# 	fq.close()
# 	# print(target_profile.loopName)
# 	# print(query_loop_list)
# 	# sys.exit()

# 	search_for_motifs_using_a_profile(target_profile, query_loop_list, directories, scores, bp_scoring_data, stk_scoring_data, nucl_scoring_data, initialized_basepair_index, include_stackings, use_branch_and_bound, dBG, output_fname)

# if __name__ == "__main__":
# 	main()
