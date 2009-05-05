#!/usr/bin/python

#################################################################
#                                                               #
#                      MAD MAPPING PROGRAM                      #
#                     ( ENTROPY  ORDERING )                     #
#                                                               #
#         COPYRIGHT  2005  2006  2007  2008  2009               #
#                        Alexander Kozik                        #
#                                                               #
#                      http://www.atgc.org/                     #
#                                                               #
#             UCD Genome Center. R.Michelmore group             #
#                                                               #
#################################################################

from __future__ import generators

def Define_Shuffle_Steps():

	### SHUFFLE IS EXPENSIVE ###
	### NUMBER OF SHUFFLE RUNS LIMITED TO EIGHT ROUNDS ###
	### EVERYTIME WHEN CURRENT MAP LENGTH IS EQUAL TO DEFINED VALUES THEN SHUFFLE TAKES PLACE ###

	global sh_r_1
	global sh_r_2
	global sh_r_3
	global sh_r_4
	global sh_r_5
	global sh_r_6
	global sh_r_7
	global sh_r_8

	sh_r_1 = 12	; # FIRST ROUND
	sh_r_2 = 16	; # + 4
	sh_r_3 = 24	; # + 8
	sh_r_4 = 36	; # +12
	sh_r_5 = 48	; # +12
	sh_r_6 = 64	; # +16
	sh_r_7 = 80	; # +16
	sh_r_8 = 96	; # +16

def Read_Data_File(in_name1, in_name2, in_name3, out_name, lk_gr_id, check_map):

	global max_frame
	global fixed_frame
	global default_diff
	global dummy_debug

	global shuffle_map
	global shuffle_final
	global shuffle_blk
	global shuffle_stp
	global round_scale

	global pwd_matrix_column
	global item_list_sorting
	global max_len_shuffle

	global sh_r_1
	global sh_r_2
	global sh_r_3
	global sh_r_4
	global sh_r_5
	global sh_r_6
	global sh_r_7
	global sh_r_8

	global out_file1
	global out_file2
	global out_file3
	global out_file4

	print "                                             "
	print "============================================="
	print " MATRIX (ALL PAIRS) :  " + in_name1
	print " MATRIX COLUMN      :  " + str(pwd_matrix_column)
	print " MARKERS  TO  MAP   :  " + in_name2
	print " WHAT TODO WITH LIST:  " + item_list_sorting
	print " FRAME MARKERS LIST :  " + in_name3
	print " OUTPUT  MAP  FILE  :  " + out_name
	print " MAX FRAME LENGTH   :  " + `max_frame`
	print " FIXED FRAME ORDER  :  " + fixed_frame
	print " LINKAGE GROUP ID   :  " + lk_gr_id
	print " DUMMY DEBUG        :  " + dummy_debug
	print " SHUFFLE MAP        :  " + shuffle_map
	print " SHUFFLE FINAL      :  " + shuffle_final
	print " SHUFFLE BLOCK      :  " + str(shuffle_blk)
	print " SHUFFLE STEP       :  " + str(shuffle_stp)
	print " MAX SHUFFLE LENGTH :  " + str(max_len_shuffle)
	print "============================================="
	print "                                             "

	time.sleep(2)

	in_file1   = open(in_name1,  "rb")
	in_file2   = open(in_name2,  "rb")
	in_file3   = open(in_name3,  "rb")
	out_file1  = open(out_name + '.mad_map_log', "wb")
	out_file2  = open(out_name + '.mad_map_final', "wb")
	out_file2M = open(out_name + '.mad_map_final.map', "wb")
	out_file3  = open(out_name + '.mad_map_temp', "wb")
	out_file4  = open(out_name + '.mad_map_xjump', "wb")
	out_file5M = open(out_name + '.matrix_2D.tab', "wb")
	out_file5A = open(out_name + '.matrix_2D.ASCII', "wb")

	out_file1.write("=============================================" + '\n')
	out_file1.write(" MATRIX (ALL PAIRS) :  " + in_name1 + '\n')
	out_file1.write(" MATRIX COLUMN      :  " + str(pwd_matrix_column) + '\n')
	out_file1.write(" MARKERS  TO  MAP   :  " + in_name2 + '\n')
	out_file1.write(" WHAT TODO WITH LIST:  " + item_list_sorting + '\n')
	out_file1.write(" FRAME MARKERS LIST :  " + in_name3 + '\n')
	out_file1.write(" OUTPUT  MAP  FILE  :  " + out_name + '\n')
	out_file1.write(" MAX FRAME LENGTH   :  " + `max_frame` + '\n')
	out_file1.write(" FIXED FRAME ORDER  :  " + fixed_frame + '\n')
	out_file1.write(" LINKAGE GROUP ID   :  " + lk_gr_id + '\n')
	out_file1.write(" DUMMY DEBUG        :  " + dummy_debug + '\n')
	out_file1.write(" SHUFFLE MAP        :  " + shuffle_map + '\n')
	out_file1.write(" SHUFFLE FINAL      :  " + shuffle_final + '\n')
	out_file1.write(" SHUFFLE BLOCK      :  " + str(shuffle_blk) + '\n')
	out_file1.write(" SHUFFLE STEP       :  " + str(shuffle_stp) + '\n')
	out_file1.write(" MAX SHUFFLE LENGTH :  " + str(max_len_shuffle) + '\n')
	out_file1.write("=============================================" + '\n')

	global best_map
	global id_frame
	global id_list
	global id_array
	global matrix_array
	global map_order
	global delta_array
	global map_array
	global nrl_array

	global depth_counter	; # RECURSIVE LOOP COUNTER

	id_frame  = []		; # LIST OF FRAME MARKERS
	id_list   = []		; # LIST OF ALL NON-REDUNDANT IDs
	id_array  = {}		; # ARRAY OF IDs ( -1 IS NOT MAPPED   +1 IS MAPPED )
	matrix_array = {}	; # DISTANCE PAIRWISE DATA FOR MARKERS
	map_order    = []	; # CURRENT MAP BEST ORDER
	delta_array  = {}	; # ARRAY OF CURRENT DELTA-DIFF
	map_array    = {}	; # ARRAY OF CURRENT ALL ORDERS
	nrl_array    = {}	; # ARRAY OF NON-REDUNDANT MAPS (REVERSE CASE IS CONSIDERED)

	out_file4.write("MARKER" + '\t' + "DELTA_ALL" + '\t' + "DELTA_PM" + '\n')

	#####################################################
	#####################################################
	print "READ MARKERS TO MAP LIST"
	time.sleep(2)
	n = 0
	while 1:
		t = in_file2.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		# tl = t.split('\t')
		####

		marker = t

		if marker in id_list:
			print ""
			print " MARKER ID DUPLICATION:   " + marker
			print ""
			sys.exit()

		if marker not in id_list:
			id_list.append(marker)
			n = n + 1
			print `n` + '\t' + marker

	print ""
	print "============================================="
	print `n` + " IDs FOUND IN MARKER LIST "
	print "============================================="
	print ""

	time.sleep(2)

	print ""
	print "============================================="
	print id_list
	print "============================================="
	print "LENGTH:  " + `len(id_list)`
	print "============================================="
	print ""

	time.sleep(2)

	#####################################################
	#####################################################

	print "READ FRAME MARKER LIST"
	time.sleep(2)
	n = 0
	while 1:
		t = in_file3.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		# tl = t.split('\t')
		####

		marker = t

		if marker not in id_list:
			print ""
			print " MARKER IS NOT IN ID LIST:   " + marker
			print ""
			sys.exit()

		if marker in id_frame:
			print ""
			print " MARKER ID DUPLICATION:   " + marker
			print ""
			sys.exit()

		if marker not in id_frame and marker in id_list:
			id_frame.append(marker)
			n = n + 1
			print `n` + '\t' + marker

	print ""
	print "============================================="
	print `n` + " IDs FOUND IN FRAME LIST "
	print "============================================="
	print ""

	time.sleep(2)

	print ""
	print "============================================="
	print id_frame
	print "============================================="
	print "LENGTH:  " + `len(id_frame)`
	print "============================================="
	print ""

	time.sleep(2)

	if len(id_frame) == 0:

		print ""
		print " NO MARKERS IN FRAME LIST:   " + `len(id_frame)`
		print ""
		sys.exit()

	if len(id_frame) > max_frame and fixed_frame == "FALSE":

		print ""
		print " TOO MANY MARKERS IN FRAME LIST:   " + `len(id_frame)`
		print ""
		sys.exit()

	if fixed_frame == "TRUE":

		print ""
		print " WORKING WITH FIXED ORDER "
		print ""
		time.sleep(2)

	if fixed_frame == "FALSE":

		print ""
		print " BEST FRAME ORDER MUST BE ESTIMATED "
		print ""
		time.sleep(2)

	#####################################################
	#####################################################

	print "  READ MATRIX FILE  "
	time.sleep(2)
	k = 0
	while 1:
		t = in_file1.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		tl = t.split('\t')
		####
		value0    = tl[0]		; # ID1
		value1    = tl[1]		; # ID2
		value2    = float(tl[pwd_matrix_column])	; # LINKAGE DATA
		# value2    = float(tl[2])	; # LINKAGE DATA
		if value0 in id_list and value1 in id_list:
			matrix_array[value0,value1] = value2
			matrix_array[value1,value0] = value2	; # SYMMETRICAL MATRIX
			# print value0 + '\t' + value1 + '\t' + str(value2)
			sys.stdout.write(`k` + " ")
			k = k + 1

	print ""
	print "============================================="

	if k > 0:
		print ""
		print `k` + " PAIRS FOUND IN MATRIX"
		print ""
		time.sleep(2)
	else:
		print ""
		print "NO DATA IN MATRIX FILE"
		print "NOTHING TO DO ....  EXIT"
		print "TERMINATED"
		print "============================================="
		print ""
		sys.exit()

	print "============================================="
	print "        PREPARE TO FRY CPU ...               "
	print "    STARTING BEST FIT ESTIMATION             "
	print "============================================="

	time.sleep(2)

	if fixed_frame == "TRUE":
		best_map = id_frame[:]

	### FRAME MAP ###
	### CHECK ALL POSSIBLE COMBINATIONS ###
 	if fixed_frame == "FALSE":

		number_of_maps = len(id_frame)

		print ""
		print " GENERATE LIST OF ALL POSSIBLE MAPS "
		print ""

		time.sleep(2)

		id_frame_sorted = id_frame[:]
		id_frame_sorted.sort()

		print " SORTED FRAME LIST: "
		print id_frame_sorted

		k = 1
		n = 1
		x = 1
		for p in xpermutations(id_frame_sorted):

			count_10000 = math.fmod(n,10000)

			# print '\t'.join(p) + '\t' + `k`
			string_map_frw = '\t'.join(p)
			p_r = p[:]
			p_r.reverse()
			string_map_rev = '\t'.join(p_r)
			try:
				query1 = nrl_array[string_map_frw]
			except:
				query1 = 0
			try:
				query2 = nrl_array[string_map_rev]
			except:
				query2 = 0
			if query1 == 0 and query2 == 0:
				map_array[k] = p
				if number_of_maps <= 8:
					print '\t'.join(map_array[k]) + '\t' + `k`
				if number_of_maps > 8 and count_10000 == 0:
					print " == " + `k` + " == " + " MAPS WERE GENERATED "
				nrl_array[string_map_frw] = 1
				nrl_array[string_map_rev] = 1
				# out_file1.write(string_map_frw + '\n')
				x = k
				k = k + 1
			n = n + 1

		print ""
		print map(''.join, list(xpermutations('ABCDE')))
		print ""

		time.sleep(2)

		print ""
		print map(''.join, list(xpermutations('DONE')))
		print ""

		print `k-1` + "  UNIQUE PERMUTATIONS GENERATED OUT OF  " +  `n-1`
		print ""

		time.sleep(2)

		print ""
		print " CHECK ALL DELTA-DIFF "
		print ""

		time.sleep(2)

		y = 1
		out_file1.write("=======" + '\n')
		while y <= x:
			print " WORKING WITH MAP:  " + `y`
			# print map_array[y]

			current_map = map_array[y]

			### FREE MEMORY ###
			del map_array[y]

			current_map_str = '\t'.join(current_map)

			### FREE MEMORY ###
			del nrl_array[current_map_str]
			current_map_rev = current_map[:]
			current_map_rev.reverse()
			current_map_rev = '\t'.join(current_map_rev)
			del nrl_array[current_map_rev]

			delta_array[current_map_str] = 0

			z1 = 0
			### ROW ###
			while z1 < len(current_map):
				item1 = current_map[z1]
				### COLUMN ###
				z2 = 0
				while z2 < len(current_map)-1:
					item2 = current_map[z2]
					item3 = current_map[z2 + 1]
					try:		
						dist1 = matrix_array[item1,item2]
						dist2 = matrix_array[item1,item3]
						# print item1 + "   " + item2 + "   " + `dist1`
						# print item1 + "   " + item3 + "   " + `dist2`
					except:
						print " NO VALUE FOR PAIR:   " + item1 + "  -=-  " + item2
						dist1 = default_diff
						dist2 = default_diff
					current_delta = math.fabs(dist1 - dist2)
					delta_array[current_map_str] = delta_array[current_map_str] + current_delta
					z2 = z2 + 1
				z1 = z1 + 1		

			if y == 1:
				smallest_delta = delta_array[current_map_str]
				best_map = current_map[:]
			if y > 1:
				if delta_array[current_map_str] < smallest_delta:
					smallest_delta = delta_array[current_map_str]
					best_map = current_map[:]
				### FREE MEMORY ###
				# if delta_array[current_map_str] > smallest_delta:
				#	del map_array[y]

			delta_per_marker = str(round(delta_array[current_map_str]/(len(current_map)),4))

			if number_of_maps <= 8:
				print current_map_str + '\t' + "***" + '\t' + str(round(delta_array[current_map_str],4)) + '\t' + "***" + '\t' + delta_per_marker
			if number_of_maps > 8:
				print "  " + str(round(delta_array[current_map_str],4)) + '\t' + "***"  + '\t' + delta_per_marker
			if dummy_debug == "TRUE":
				out_file1.write(current_map_str + '\t' + "***" + '\t' + str(round(delta_array[current_map_str],4)) \
					+ '\t' + "***" + '\t' + delta_per_marker + '\t' + "***" + '\t' + `y` + '\n')

			### FREE MEMORY ###
			del delta_array[current_map_str]

			y = y + 1

		out_file1.write("=======" + '\n')
		m = 0
		out_file3.write("=LG=" + '\t' + "=MARKER=" + '\t' + "=POS=" + '\t' +  '\n')
		for marker in best_map:
			out_file3.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\n')
			m = m + 1
		out_file3.write("=XX=" + '\t' + "=XXXXXX=" + '\t' + "=XXX=" + '\t' +  '\n')


	### MAP CONSTRUCTION ###

	time.sleep(2)

	if check_map == "TRUE":
		best_map = id_list

	print ""
	print ""
	print "ORIGINAL LIST: "
	print id_list
	print ""
	time.sleep(2)


	### WHAT TO DO WITH ITEM LIST - SORT - SHUFFLE - AS_IS ###
	if check_map != "TRUE" and item_list_sorting == "ABC_ORDER":
		id_list.sort()
	if check_map != "TRUE" and item_list_sorting == "DONT_TOUCH":
		mooba = "booba"
	if check_map != "TRUE" and item_list_sorting == "SHUFFLE_ME":
		random.seed()
		random.shuffle(id_list)

	print item_list_sorting
	print ""

	print "WORKING LIST: "
	print id_list
	print ""
	time.sleep(2)

	m = 1
	for item in id_list:

		### NEXT ITEM ###
		if item not in best_map:
			print " MAPPING MARKER :  "  + item + "  |  " +  `m`
			s = 1
			while s <= len(best_map)+1:
				current_map = best_map[:]
				current_map.insert(s-1,item)
				map_array[s] = current_map
				# print map_array[s]
				# print " MAP: " + `m` + "    VERSION: " + `s`
				s = s + 1
			############################################
			###            FIND BEST MAP             ###
			print " ==================================== "
			print " WORKING WITH MAP SERIES:   " + `m`
			print " ==================================== "
			time.sleep(1)
			Find_Best_Map(s, lk_gr_id, item)
			############################################
			### SHUFFLE ###
			map_len = len(best_map)
			if map_len == sh_r_1 or map_len == sh_r_2 or map_len == sh_r_3 or map_len == sh_r_4 or map_len == sh_r_5 or map_len == sh_r_6 or map_len == sh_r_7 or map_len == sh_r_8:
				if shuffle_map == "TRUE":
					print "============================="
					print "  IT IS TIME TO SHUFFLE !!!  "
					print "        " + str(m) + "       "
					print " MAP LENGTH:  " + str(map_len)
					print "============================="
					Shuffle_Best_Map()
					#################################################
					# m = 0
					# out_file3.write("=LG=" + '\t' + "=MARKER=" + '\t' + "=POS=" + '\t' + "AFTER_SHUFFLING" + '\n')
					# for marker in best_map:
					#	out_file3.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\n')
					#	m = m + 1
					# out_file3.write("=XX=" + '\t' + "=XXXXXX=" + '\t' + "=XXX=" + '\t' + "===============" + '\n')
					#################################################

				if shuffle_map != "TRUE":
					print "==================================="
					print "  ... I DO NOT WANT TO SHUFFLE ... "
					print " MAP LENGTH:  " + str(map_len)
					print "==================================="
			############################################
			m = m + 1

	### LAST ROUND OF SHUFFLING ###
	if shuffle_final == "TRUE":
		Shuffle_Best_Map()
		print "               "
		print " FINAL SHUFFLE "
		print "               "
	if shuffle_final != "TRUE":
		print "==================================="
		print "  ... I DO NOT WANT TO SHUFFLE ... "
		print "==================================="
	### BEST MAP ###
	m = 0
	map_position = 0
	out_file2.write("LG" + '\t' + " MARKER " + '\t' + " POS " + '\t' + "#1#" + '\t' + " DST1 " \
				+ '\t' + "#2#" + '\t' + " DST2 " + '\t' + "#3#" + '\t' + " DST3 " + '\t' \
				+ "#S#" + '\t' + " SUMM " + '\t' + "#D#" + '\t' + " DIFF " + '\t' + "STATUS" + '\t' + "CLASS" + '\n')
	for marker in best_map:
		if m == 0:
			### FIRST ROW ###
			out_file2M.write(lk_gr_id + '\t' + marker + '\t' + str(round(map_position,2)) + '\n')
			out_file2.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\t' + "#1#" + '\t' + `m` \
				+ '\t' + "#2#" + '\t' + "NNNNNN" + '\t' + "#3#" + '\t' + "NNNNNN" + '\t' \
				+ "#S#" + '\t' + "NNNNNN" + '\t' + "#D#" + '\t' + "NNNNNN" + '\t' + "NNNNNN" + '\t' + "NNNNN" + '\n')
		if m > 0 and m < len(best_map)-1:
			### EVERYTHING ELSE ###
			dst1 = "YYYYYY"
			dst2 = "YYYYYY"
			dst3 = "YYYYYY"
			delta_status = "GOOD"
			delta_class  = "YYYYY"
			### MARKER B - MIDDLE MARKER ###
			###        A - ABOVE         ###
			###        C - BELOW         ###
			try:
				### MARKER ABOVE ###
				mrkA = best_map[m-1]
				mrkB = best_map[m]
				dst1 = matrix_array[mrkA,mrkB]
				dst1s = str(round(dst1,4))
				map_position = map_position + dst1*100
			except:
				dst1s = "XX11XX"
				delta_status = "#BAD"
				map_position = map_position + 50
			try:
				### MARKER BELOW ###
				mrkB = best_map[m]
				mrkC = best_map[m+1]
				dst2 = matrix_array[mrkB,mrkC]
				dst2s = str(round(dst2,4))
			except:
				dst2s = "XX22XX"
				delta_status = "#BAD"
			try:
				### FLANKING MARKERS ###
				mrkA = best_map[m-1]
				mrkC = best_map[m+1]
				dst3 = matrix_array[mrkA,mrkC]
				dst3s = str(round(dst3,4))
			except:
				dst3s = "XX33XX"
				delta_status = "#BAD"
			if delta_status == "GOOD":
				summ12 = dst1 + dst2
				diff123 = summ12 - dst3
				summ12s = str(round(summ12,4))
				diff123s = str(round(diff123,4))
				if math.fabs(diff123) > 0.50:
					delta_class = "=HUGE="
				if math.fabs(diff123) > 0.40 and math.fabs(diff123) <= 0.50:
					delta_class = "LARGE"
				if math.fabs(diff123) > 0.30 and math.fabs(diff123) <= 0.40:
					delta_class = "_XXX_"
				if math.fabs(diff123) > 0.20 and math.fabs(diff123) <= 0.30:
					delta_class = "__X__"
				if math.fabs(diff123) > 0.10 and math.fabs(diff123) <= 0.20:
					delta_class = "__1__"
				if math.fabs(diff123) <= 0.10:
					delta_class = "__0__"
			if delta_status == "#BAD":
				summ12s = "XXSSXX"
				diff123s = "XXDDXX"
				delta_class  = "XXXXX"
			out_file2M.write(lk_gr_id + '\t' + marker + '\t' + str(round(map_position,2)) + '\n')
			out_file2.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\t' + "#1#" + '\t' + dst1s \
				+ '\t' + "#2#" + '\t' + dst2s + '\t' + "#3#" + '\t' + dst3s + '\t' \
				+ "#S#" + '\t' + summ12s + '\t' + "#D#" + '\t' + diff123s + '\t' + delta_status + '\t' + delta_class + '\n')
		if m == len(best_map)-1:
			### LAST  ROW ###
			dst1 = "YYYYYY"
			try:
				mrkA = best_map[m-1]
				mrkB = best_map[m]
				dst1 = matrix_array[mrkA,mrkB]
				dst1s = str(round(dst1,4))
			except:
				dst1s = "XX11XX"
			out_file2M.write(lk_gr_id + '\t' + marker + '\t' + str(round(map_position,2)) + '\n')
			out_file2.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\t' + "#1#" + '\t' + dst1s \
				+ '\t' + "#2#" + '\t' + "NNNNNN" + '\t' + "#3#" + '\t' + "NNNNNN" + '\t' \
				+ "#S#" + '\t' + "NNNNNN" + '\t' + "#D#" + '\t' + "NNNNNN" + '\t' + "NNNNNN" + '\t' + "NNNNN" + '\n')
		m = m + 1
	# out_file2.write("=XX=" + '\t' + "=XXXXXX=" + '\t' + "=XXX=" + '\t' +  '\n')

	### TWO DIMENSIONAL MATRIX ###

	out_file5M.write(";" + '\t')

	mega_length = len(best_map)

	### FIRST ROW ###
	counter = 0
	for item in best_map:
		current_id = item
		out_file5M.write(current_id)
		counter = counter + 1
		if counter < mega_length:
			out_file5M.write('\t')
		if counter == mega_length:
			out_file5M.write('\n')

	### PAIRWISE DATA ROWS ###
	for item in best_map:
		current_id = item
		out_file5M.write(current_id + '\t')
		out_file5A.write(current_id + '\t')
		counter = 0
		for other_item in best_map:
			counter = counter + 1
			other_id = other_item
			try:
				cur_diff = float(matrix_array[current_id,other_id])
				cur_diff_string = str(round(cur_diff,round_scale))
				cur_diff_fr = cur_diff
				cur_value = str(round(cur_diff_fr,round_scale))
				cur_value_ascii = (cur_diff)*10
				cur_value_ascii = round(cur_value_ascii,0)
				cur_value_ascii = int(cur_value_ascii)
				if cur_value_ascii >= 10:
					cur_value_ascii = 9
				cur_value_ascii = str(cur_value_ascii)
			except:
				cur_value = "-"
				cur_value_ascii = "-"
			out_file5M.write(cur_value)
			out_file5A.write(cur_value_ascii)
			if counter < mega_length:
				out_file5M.write('\t')
			if counter == mega_length:
				out_file5M.write('\n')
				out_file5A.write('\n')


	print "                         "
	print "     WELL DONE ...       "
	print "    ENJOY NEW ORDER      "
	print "                         "

	### THE  END ###

	in_file1.close()
	in_file2.close()
	in_file3.close()
	out_file1.close()
	out_file2.close()
	out_file2M.close()
	out_file3.close()
	out_file4.close()
	out_file5M.close()
	out_file5A.close()

	################

def Shuffle_Best_Map():

	global out_file1
	global out_file3

	global matrix_array

	global shuffle_map
	global shuffle_blk
	global shuffle_stp

	global best_map

	global shfl_array
	global dummy_debug

	shfl_array = {}

	print ""
	print "======================================"
	print "  ... SHUFFLE - SHUFFLE - SHUFFLE ... "
	print "======================================"
	print ""

	if dummy_debug == "TRUE":
		out_file1.write("  ... SHUFFLE - SHUFFLE - SHUFFLE ... " + '\n')

	time.sleep(2)

	current_map_length = len(best_map)

	st = 0
	ed = shuffle_blk

	while ed <= current_map_length:

		current_block = best_map[st:ed]
		current_head  = best_map[:st]
		current_tail  = best_map[ed:]

		print "---------------------------------"
		print current_head
		print current_block
		print current_tail
		print "---------------------------------"
		print ""

		time.sleep(2)

		#########################################
		#########################################
		k = 1
		n = 1
		x = 1
		for p in xpermutations(current_block):
			shuffled_map = current_head + p + current_tail
			string_map = '\t'.join(shuffled_map)
			shfl_array[k] = shuffled_map
			# print string_map + '\t' + "***" + '\t' + `k`
			print `k` + '\t' + "PERMUTATED MAP"

			x = k
			k = k + 1
			n = n + 1

		print ""
		print map(''.join, list(xpermutations('DONE')))
		print ""

		print `k-1` + "  UNIQUE PERMUTATIONS GENERATED OUT OF  " +  `n-1`
		print ""

		time.sleep(2)

		print ""
		print " CHECK ALL DELTA-DIFF "
		print ""

		time.sleep(2)

		y = 1
		# out_file1.write("=======" + '\n')
		while y <= x:
			print " WORKING WITH MAP:  " + `y`

			current_map = shfl_array[y]

			### FREE MEMORY ###
			del shfl_array[y]

			current_map_str = '\t'.join(current_map)

			delta_array[current_map_str] = 0

			z1 = 0
			### ROW ###
			while z1 < len(current_map):
				item1 = current_map[z1]
				### COLUMN ###
				z2 = 0
				while z2 < len(current_map)-1:
					item2 = current_map[z2]
					item3 = current_map[z2 + 1]
					try:		
						dist1 = matrix_array[item1,item2]
						dist2 = matrix_array[item1,item3]
						# print item1 + "   " + item2 + "   " + `dist1`
						# print item1 + "   " + item3 + "   " + `dist2`
					except:
						sys.stdout.write(".")
						# print " NO VALUE FOR PAIR:   " + item1 + "  -=-  " + item2
						dist1 = default_diff
						dist2 = default_diff
					current_delta = math.fabs(dist1 - dist2)
					delta_array[current_map_str] = delta_array[current_map_str] + current_delta
					z2 = z2 + 1
				z1 = z1 + 1		

			if y == 1:
				smallest_delta = delta_array[current_map_str]
				best_map = current_map[:]
			if y > 1:
				if delta_array[current_map_str] < smallest_delta:
					smallest_delta = delta_array[current_map_str]
					best_map = current_map[:]

			delta_per_marker = str(round(delta_array[current_map_str]/(len(current_map)),4))

			print "  " + str(round(delta_array[current_map_str],4)) + '\t' + "***"  + '\t' + delta_per_marker

			### FREE MEMORY ###
			del delta_array[current_map_str]

			y = y + 1

		### LOG FILE RECORD AFTER SHUFFLING ###

		current_best_str = '\t'.join(best_map)
		delta_per_marker_best = str(round(smallest_delta/(len(best_map)),4))

		if dummy_debug == "TRUE":
			out_file1.write(current_best_str + '\t' + "***" + '\t' + str(round(smallest_delta,4)) \
				+ '\t' + "***" + '\t' + delta_per_marker_best + '\t' + "***" + '\t' + `y` + '\n')


		#################################################
		m = 0
		out_file3.write("=LG=" + '\t' + "=MARKER=" + '\t' + "=POS=" + '\t' + "AFTER_SHUFFLING" + '\n')
		for marker in best_map:
			if m == 0:
				out_file3.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\t' + "=->" + '\t' \
					+ str(round(smallest_delta,4)) + '\t' + delta_per_marker_best + '\n')
			if m > 0:
				out_file3.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\n')
			m = m + 1
		out_file3.write("=XX=" + '\t' + "=XXXXXX=" + '\t' + "=XXX=" + '\t' + "===============" + '\n')
		#################################################

		#########################################
		#########################################		

		ed = ed + shuffle_stp
		st = st + shuffle_stp


def Find_Best_Map(s, lk_gr_id, itemx):

	global matrix_array
	global map_array
	global best_map
	global default_diff

	global out_file1
	global out_file2
	global out_file3
	global out_file4

	x = s - 1

	y = 1
	if dummy_debug == "TRUE":
		out_file1.write("=======" + '\n')
	while y <= x:
		print ""
		print " WORKING WITH TEMP MAP:  " + `y`
		# print map_array[y]

		current_map = map_array[y]

		### FREE MEMORY ###
		del map_array[y]

		current_map_str = '\t'.join(current_map)

		delta_array[current_map_str] = 0

		z1 = 0
		### ROW ###
		while z1 < len(current_map):
			item1 = current_map[z1]
			### COLUMN ###
			z2 = 0
			while z2 < len(current_map)-1:
				item2 = current_map[z2]
				item3 = current_map[z2 + 1]
				try:		
					dist1 = matrix_array[item1,item2]
					dist2 = matrix_array[item1,item3]
					# print item1 + "   " + item2 + "   " + `dist1`
					# print item1 + "   " + item3 + "   " + `dist2`
				except:
					# print " NO VALUE FOR PAIR:   " + item1 + "  -=-  " + item2
					sys.stdout.write(".")
					dist1 = default_diff
					dist2 = default_diff
				current_delta = math.fabs(dist1 - dist2)
				delta_array[current_map_str] = delta_array[current_map_str] + current_delta
				z2 = z2 + 1
			z1 = z1 + 1		

		delta_per_marker = str(round(delta_array[current_map_str]/(len(current_map)),4))

		if y == 1:
			smallest_delta = delta_array[current_map_str]
			best_map = current_map[:]
			best_delta_per_marker = str(round(delta_array[current_map_str]/(len(current_map)),4))
		if y > 1:
			if delta_array[current_map_str] < smallest_delta:
				smallest_delta = delta_array[current_map_str]
				best_map = current_map[:]
				best_delta_per_marker = str(round(delta_array[current_map_str]/(len(current_map)),4))

		# delta_per_marker = str(round(delta_array[current_map_str]/(len(current_map)),4))

		# print current_map_str + '\t' + "***" + '\t' + str(round(delta_array[current_map_str],4)) + '\t' + "***" + '\t' + delta_per_marker
		if dummy_debug == "TRUE":
			out_file1.write(current_map_str + '\t' + "***" + '\t' + str(round(delta_array[current_map_str],4)) \
				+ '\t' + "***" + '\t' + delta_per_marker + '\t' + "***" + '\t' + `y` + '\n')

		### FREE MEMORY ###
		del delta_array[current_map_str]

		y = y + 1

	# out_file1.write("=======" + '\n')
	m = 0
	print "                      "
	print " ==================== "
	print " BEST CURRENT MAP :   "
	print "                      "
	print best_map
	out_file4.write(itemx + '\t' + str(round(smallest_delta,4)) + '\t' + best_delta_per_marker + '\n')
	out_file3.write("=LG=" + '\t' + "=MARKER=" + '\t' + "=POS=" + '\t' +  '\n')
	time.sleep(1)
	for marker in best_map:
		marker_label = "---"
		sm_delta_label = " - "
		pm_delta_label = " - "
		if marker == itemx:
			marker_label = "<-="
			sm_delta_label = str(round(smallest_delta,4))
			pm_delta_label = best_delta_per_marker
		out_file3.write(lk_gr_id + '\t' + marker + '\t' + `m` + '\t' + marker_label \
			+ '\t' + sm_delta_label + '\t' + pm_delta_label + '\n')
		m = m + 1
	out_file3.write("=XX=" + '\t' + "=XXXXXX=" + '\t' + "=XXX=" + '\t' +  '\n')

	print ""
	print "  CONTINUE WITH ... "
	print ""

###    PERMUTATIONS     ###
def xcombinations(items, n):
	if n==0: yield []
	else:
		for i in xrange(len(items)):
			for cc in xcombinations(items[:i]+items[i+1:],n-1):
				yield [items[i]]+cc

def xpermutations(items):
	return xcombinations(items, len(items))

### END OF PERMUTATIONS ###

### MAIN BODY ###

import math
import re
import sys
import string
import time
import random

# from __future__ import generators

global max_frame
global fixed_frame
global default_diff
global dummy_debug
global shuffle_map
global shuffle_final
global shuffle_blk
global shuffle_stp
global round_scale
global max_len_shuffle

### COLUMN WITH PW VALUES ###
### COUNT FROM 0 : 2 == 3 ###
global pwd_matrix_column
pwd_matrix_column = 2

### WHAT TO DO WITH INITIAL ITEM LIST ###
###  1. RANDOM SHUFFLE - DEFAULT      ### 
###  2. NOTHING - ORDERED AS IS       ###
###  3. SORT BY ITEM ID IN ABC ORDER  ###
global item_list_sorting
item_list_sorting = "SHUFFLE_ME"
# item_list_sorting = "DONT_TOUCH"
# item_list_sorting = "ABC_ORDER"

### MAX LENGTH OF LIST TO SHUFFLE ###
### IF LIST LARGER THAN MAX THEN SHUFLING DOES NOT TAKE PLACE ###
# max_len_shuffle = 80
max_len_shuffle = 120
# max_len_shuffle = 160

### MAX LENGTH OF FRAME LIST ###
max_frame = 10
default_diff = 0

### TO CHECK OR NOT TO CHECK FRAME ORDER ###
fixed_frame = "FALSE"
# fixed_frame = "TRUE"

# shuffle_final = "TRUE"
shuffle_final = "FALSE"

dummy_debug = "FALSE"
# dummy_debug = "TRUE"

shuffle_blk = 6
shuffle_stp = 3
round_scale = 6

check_map = "FALSE"

if __name__ == "__main__":
	if len(sys.argv) <= 9 or len(sys.argv) > 10:
		print "                                                                                "
		print " Program usage:                                                                 "
		print " [matrix(all_pairs)] [items_to_order] [frame_list] [output_file] [link_gr_id] [frame_fixed/frame_flex] [shuffle/noshuffle] [shuffle_block] [shuffle_step]"
		print " frame_fixed == \"FIXED\"    frame_flexible == \"FLEX\"      "
		print " if frame is \"FIXED\" then it will be used with fixed order                    "
		print " if frame is \"FLEX\" (flexible) then all possible combinations will be checked "
		print " shuffle/noshuffle option is \"SHUFFLE\" or \"NOSHUFFLE\"                       "
		print " shuffle_block is 6, 7 or 8                                                     "
		print " shuffle_step is 3, 4, 5 or 6                                                   "
		print " default recommended shuffle options are: [ SHUFFLE 6 3 ]                       "
		print "                                                                                "
		exit
	if len(sys.argv) == 10:
		in_name1  = sys.argv[1]
		in_name2  = sys.argv[2]
		in_name3  = sys.argv[3]
		out_name  = sys.argv[4]
		lk_gr_id  = sys.argv[5]
		fix_frm   = sys.argv[6]

		shuffle_opt1 = sys.argv[7]
		shuffle_opt2 = sys.argv[8]
		shuffle_opt3 = sys.argv[9]

		if shuffle_opt1 == "SHUFFLE":
			shuffle_map = "TRUE"

		if shuffle_opt1 == "NOSHUFFLE":
			shuffle_map = "FALSE"

		if shuffle_opt1 != "SHUFFLE" and shuffle_opt1 != "NOSHUFFLE" and shuffle_opt1 != "SHUFFLE_FINAL":
			print ""
			print " argument must be \"SHUFFLE\" or \"NOSHUFFLE\" "
			print ""
			sys.exit()

		if int(shuffle_opt2) >= 6 and int(shuffle_opt2) <= 8:
			shuffle_blk = int(shuffle_opt2)

		if int(shuffle_opt3) >= 3 and int(shuffle_opt3) <= 6:
			shuffle_stp = int(shuffle_opt3)

		if fix_frm == "FIXED":
			fixed_frame = "TRUE"

		if fix_frm == "CHECK_MAP":
			fixed_frame = "TRUE"
			check_map = "TRUE"

		if fix_frm == "FLEX":
			fixed_frame = "FALSE"


		if fix_frm != "FIXED" and fix_frm != "FLEX" and fix_frm != "CHECK_MAP":
			print ""
			print "[frame_fixed/frame_flex] argument must be \"FIXED\" or \"FLEX\" or \"CHECK_MAP\" "
			print ""
			sys.exit()

		sys.setrecursionlimit(1000000000)

		Define_Shuffle_Steps()
		Read_Data_File(in_name1, in_name2, in_name3, out_name, lk_gr_id, check_map)

#### THE END ####
