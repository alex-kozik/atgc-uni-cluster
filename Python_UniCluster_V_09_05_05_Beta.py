#!/usr/bin/python

#################################################################
#                                                               #
#                          MAD MAPPER                           #
#                                                               #
#                      MAD MAPPING PROGRAM                      #
#                                                               #
#                      DIGITAL  CLUSTERING                      #
#                                                               #
#         COPYRIGHT  2004 2005 2006 2007 2008 2009              #
#                        Alexander Kozik                        #
#                                                               #
#                      http://www.atgc.org/                     #
#                                                               #
#             UCD Genome Center. R.Michelmore group             #
#                                                               #
#################################################################


##################################################################################
#1   5   10   15   20   25   30   35   40   45   50   55   60   65   70   75   80#
#-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5-=3+5#
##################################################################################


##################################################################################
#                           --+============+--                                   #
#                              README_FIRST:                                     #
#                            PYTHON MAD MAPPER                                   #
#                        HOW IT WORKS AND WHAT IT DOES:                          #
#                  RULES AND ALGORITHMS OF MADMAPPER APPROACH                    #
#             ---++==========================================++---               #
#                                                                                #
#                                                                                #
#  = HOW MADMAPPER PAIRWISE DISTANCES ARE CALCULATED =                           #
#                                                                                #
#  Example dataset:                                                              #
#                                                                                #
#  ;     1    2    3    4    5    6    7    8    9    10                         #
#                                                                                #
#  MA    1    1    1    1    1    1    1    1    1    1                          #
#  MB    0    0    0    0    0    0    0    0    1    1                          #
#  MC    0    0    0    0    0    0    1    1    1    1                          #
#  MD    0    0    0    0    0    0    -    -    1    1                          #
#                                                                                #
#  ABS_DIFF  - ABSOLUTE DIFFERENCE                                               #
#  DIFF_PP   - DIFFERENCE PER POINT                                              #
#  NORM_DIFF - NORMALIZE DIFFERENCE                                              #
#  (largest difference considered as maximum with assigned value '1')            #
#                                                                                #
#  PAIR    ABS_DIFF     DIFF_PP           NORM_DIFF                              #
#                                                                                #
#  A - B      8       (8/10)=0.80      (0.80/0.80)=1.00 - MAX_VALUE              #
#  A - C      6       (6/10)=0.60      (0.60/0.80)=0.75                          #
#  A - D      6       (6/ 8)=0.75      (0.75/0.80)=0.94                          #
#  B - C      2       (2/10)=0.20      (0.20/0.80)=0.25                          #
#  B - D      0       (0/ 8)=0.00      (0.00/0.80)=0.00                          #
#  C - D      0       (0/ 8)=0.00      (0.00/0.80)=0.00                          #
#                                                                                #
#  TWO DIMENSIONAL MATRIX FILE                                                   #
#                                                                                #
#         MA      MB     MC     MD                                               #
#    MA  0.00    1.00   0.75   0.94                                              #
#    MB  1.00    0.00   0.25   0.00                                              #
#    MC  0.75    0.25   0.00   0.00                                              #
#    MD  0.94    0.00   0.00   0.00                                              #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#                                                                                #
#  DESCRIPTION OF DISTANCE MATRIX FILES                                          #
#                                                                                #
#  *.pairs_all                                                                   #
#                                                                                #
#  [ 1 ]   [ 2 ]   [ 3 ]        [ 4 ]          [ 5 ]           [ 6 ]   [7]   [8] #
# M_A_01  M_B_55  116.5  *[DFR>* 1.0    *[VPP>* 2.33      *[PCC>* 50    0    50  #
# M_A_01  M_B_56  110.2  *[DFR>* 0.96   *[VPP>* 2.295833  *[PCC>* 48    2    50  #
# M_A_01  M_B_57  116.5  *[DFR>* 1.0    *[VPP>* 2.33      *[PCC>* 50    0    50  #
# M_A_01  M_B_58  113.1  *[DFR>* 0.98   *[VPP>* 2.308163  *[PCC>* 49    1    50  #
#                                                                                #
# [ 1 ] - item 1                                                                 #
# [ 2 ] - item 2                                                                 #
# [ 3 ] - absolute difference between item 1 and 2                               #
# [ 4 ] - data fraction (valid datapoints / max number of possible datapoints)   #
#         ( abbreviation: *[DFR>* - data fraction )                              #
# [ 5 ] - difference per point                                                   #
#         ( abbreviation: *[VPP>* - value per point )                            #
# [ 6 ] - number of valid datapoints (paired comparisons)                        #
#         ( abbreviation: *[PCC>* - pairwise comparison count )                  #
# [ 7 ] - number of missing comparisons                                          #
# [ 8 ] - max number of possible datapoints                                      #
#                                                                                #
#                                                                                #
#  *.pairs_all.2D_Matrix -                                                       #
#   - 2D Matrix with pairwise distances normalized between 0 and 1               #
#                                                                                #
#  *.pairs_all.2D_ASCII -                                                        #
#   - compact form of 2D Matrix with values normalized between 1 and 9           #
#                                                                                #
#                                                                                #
#  *.pairs_all.PW_0_1.tab -                                                      #
#   - distance matrix file with values normalized between 0 and 1                #
#                                                                                #
#  *.pairs_all.PW_Val.tab -                                                      #
#   - distance matrix file with real (non-normalized) values per point           #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#                                                                                #
#  DESCRIPTION OF 24 GROUPING/CLUSTERING ITERATIONS                              #
#  AND HOW TO CHOOSE/DEFINE CUTOFF VALUES                                        #
#                                                                                #
#  1. read into memory input file                                                #
#  2. pairwise comparisons                                                       #
#  3. creation of pairwise distance matrix                                       #
#  4. twenty four iterations (rounds) of grouping (clustering)                   #
#  5. dendro-sorting of all iterations                                           #
#  6. writing of output files                                                    #
#                                                                                #
#  cutoff values for iterations are defined under functions:                     #
#                                "Set_CutOff_Values_Type_#"                      #
#                                                                                #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#                                                                                #
#  = CLUSTERING/GROUPING OUTPUT FILES =                                          #
#     Information about grouping is stored in three files per iteration:         #
#     '*.matrix'     - pairwise distances for a given group                      #
#     '*.adj_list'   - adjacency list                                            #
#     '*.group_info' - group info                                                #
#                                                                                #
#     STRUCTURE OF '*.matrix' FILE:                                              #
#    * first column:   item ID "A"                                               #
#    * second column:  item ID "B"                                               #
#    * third column:   distance value between item "A" and "B"                   #
#    * fourth column:  number of datapoints (comparisons)                        #
#    * fifth column:   normalized (0-1) distance value between item "A" and "B"  #
#                                                                                #
#     STRUCTURE OF '*.group_info' FILE:                                          #
#    * first column:   item ID                                                   #
#    * second column:  length of an adjacency list for given item or how many    #
#                      other items are linked to given marker directly           #
#    * third column:   size of the given group (how many items in this           #
#                      particular linked group)                                  #
#    * fourth column:  arbitrary group number                                    #
#    * fifth column:   visual mark ("*****" separates different group)           #
#    * sixth column:   information about framework items                         #
#    * seventh column: type of graph [SINGLETON/LINKED/COMPLETE].                #
#      If node is a singleton (is not connected to any other node) then it is    #
#      labeled as 'SINGLE____NODE'                                               #
#      If nodes form complete graph (all nodes linked to each other directly)    #
#      then such group is labeled as 'COMPLETE_GRAPH'                            #
#      If group is not a complete graph (some nodes do not have direct links     #
#      or connections to other nodes then such group is labeled as               #
#      'LINKED___GROUP'.                                                         #
#    * eighth column:  type of node (SATURATED or DILUTED). If a node has all    #
#      possible connections to all other nodes in a group [in other words: node  #
#      is connected directly to all other nodes in a group] then such node is    #
#      labeled as 'SATURATED_NODE'. 'DILUTED___NODE' is an indication that a     #
#      group does not form complete graph. Group can be considered as a bin if   #
#      all nodes in a given group are SATURATED and graph is COMPLETE.           #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#                                                                                #
#  = CLUSTERING/GROUPING SUMMARY FOR ALL 24 ITERATIONS =                         #
#     [ DENDRO-CLUSTERING ]                                                      #
#     STRUCTURE OF *.x_tree_clust FILE:                                          #
#    * 1-st column: group ID for clustering with cutoff val_1                    #
#    * 2-nd column: group ID for clustering with cutoff val_2                    #
#    * 3-d  column: group ID for clustering with cutoff val_3                    #
#    * 4-th column: group ID for clustering with cutoff val_4                    #
#    * 5-th column: group ID for clustering with cutoff val_5                    #
#    * 6-th column: group ID for clustering with cutoff val_6                    #
#    * 7-th column: group ID for clustering with cutoff val_7                    #
#    * 8-th column: group ID for clustering with cutoff val_8                    #
#    * 9-th column: group ID for clustering with cutoff val_9                    #
#    *10-th column: group ID for clustering with cutoff val_A                    #
#    *11-th column: group ID for clustering with cutoff val_B                    #
#    *12-th column: group ID for clustering with cutoff val_C                    #
#    *13-th column: group ID for clustering with cutoff val_D                    #
#    *14-th column: group ID for clustering with cutoff val_E                    #
#    *15-th column: group ID for clustering with cutoff val_F                    #
#    *16-th column: group ID for clustering with cutoff val_G                    #
#    *17-th column: group ID for clustering with cutoff val_H                    #
#    *18-th column: group ID for clustering with cutoff val_I                    #
#    *19-th column: group ID for clustering with cutoff val_J                    #
#    *20-th column: group ID for clustering with cutoff val_K                    #
#    *21-th column: group ID for clustering with cutoff val_L                    #
#    *22-th column: group ID for clustering with cutoff val_M                    #
#    *23-th column: group ID for clustering with cutoff val_N                    #
#    *24-th column: group ID for clustering with cutoff val_O                    #
#    ************************************************************************    #
#                                                                                #
#    *25-th column: type of graph (COMPLETE or LINKED) for the last iteration    #
#    *26-th column: type of node (SATURATED or DILUTED) for the last iteration   #
#    *27-th column: '***' - visual mark                                          #
#    *28-th column: group class from frame work item list                        #
#    *29-th column: coordinates of frame item (if available)                     #
#    *30-th column: '***' - visual mark                                          #
#    *31-st column: "ABC" (alphabetical) order of items prior clustering         #
#    *32-nd column: '***' - visual mark                                          #
#    *33-d  column: "LG" - reserved field for manipulation in excel like editor  #
#    *34-th column: item ID                                                      #
#    *35-th column: '*' - visual mark                                            #
#    *36-th column: number of scores below '0' (negative)                        #
#    *37-th column: number of scores equal to '0' (neutral)                      #
#    *38-th column: number of scores above '0' (positive)                        #
#    *39-th column: '*' visual mark                                              #
#    *40-th column: total number of possible scores                              #
#    *41-st column: missing datapoints                                           #
#    *42-nd column: number of fixed values                                       #
#    *43-d  column: '*' - visual mark                                            #
#    *44-th column: order of markers according to dendro-clustering              #
#                                                                                #
#  ____________________________________________________________________________  #
#                                                                                #
#                                                                                #
#                                                                                #
##################################################################################

def Set_CutOff_Values_Type_1():

	### GLOBAL ###
	global cutoff_value_1
	global cutoff_value_2
	global cutoff_value_3
	global cutoff_value_4
	global cutoff_value_5
	global cutoff_value_6

	global cutoff_value_7
	global cutoff_value_8
	global cutoff_value_9
	global cutoff_value_A
	global cutoff_value_B
	global cutoff_value_C

	global cutoff_value_D
	global cutoff_value_E
	global cutoff_value_F
	global cutoff_value_G
	global cutoff_value_H
	global cutoff_value_I

	global cutoff_value_J
	global cutoff_value_K
	global cutoff_value_L
	global cutoff_value_M
	global cutoff_value_N
	global cutoff_value_O

	### VALUES ###
	cutoff_value_1 = 0.800
	cutoff_value_2 = 0.600
	cutoff_value_3 = 0.500
	cutoff_value_4 = 0.400
	cutoff_value_5 = 0.300
	cutoff_value_6 = 0.250

	cutoff_value_7 = 0.200
	cutoff_value_8 = 0.180
	cutoff_value_9 = 0.160
	cutoff_value_A = 0.140
	cutoff_value_B = 0.120
	cutoff_value_C = 0.100

	cutoff_value_D = 0.090
	cutoff_value_E = 0.080
	cutoff_value_F = 0.070
	cutoff_value_G = 0.060
	cutoff_value_H = 0.050
	cutoff_value_I = 0.040

	cutoff_value_J = 0.030
	cutoff_value_K = 0.020
	cutoff_value_L = 0.015
	cutoff_value_M = 0.010
	cutoff_value_N = 0.005
	cutoff_value_O = 0.000

def Set_CutOff_Values_Type_2():

	### GLOBAL ###
	global cutoff_value_1
	global cutoff_value_2
	global cutoff_value_3
	global cutoff_value_4
	global cutoff_value_5
	global cutoff_value_6

	global cutoff_value_7
	global cutoff_value_8
	global cutoff_value_9
	global cutoff_value_A
	global cutoff_value_B
	global cutoff_value_C

	global cutoff_value_D
	global cutoff_value_E
	global cutoff_value_F
	global cutoff_value_G
	global cutoff_value_H
	global cutoff_value_I

	global cutoff_value_J
	global cutoff_value_K
	global cutoff_value_L
	global cutoff_value_M
	global cutoff_value_N
	global cutoff_value_O

	### VALUES ###
	cutoff_value_1 = 0.600
	cutoff_value_2 = 0.500
	cutoff_value_3 = 0.400
	cutoff_value_4 = 0.350
	cutoff_value_5 = 0.300
	cutoff_value_6 = 0.250

	cutoff_value_7 = 0.225
	cutoff_value_8 = 0.200
	cutoff_value_9 = 0.175
	cutoff_value_A = 0.150
	cutoff_value_B = 0.125
	cutoff_value_C = 0.100

	cutoff_value_D = 0.090
	cutoff_value_E = 0.080
	cutoff_value_F = 0.070
	cutoff_value_G = 0.060
	cutoff_value_H = 0.050
	cutoff_value_I = 0.040

	cutoff_value_J = 0.030
	cutoff_value_K = 0.020
	cutoff_value_L = 0.015
	cutoff_value_M = 0.010
	cutoff_value_N = 0.005
	cutoff_value_O = 0.000

def Set_CutOff_Values_Type_3():

	### GLOBAL ###
	global cutoff_value_1
	global cutoff_value_2
	global cutoff_value_3
	global cutoff_value_4
	global cutoff_value_5
	global cutoff_value_6

	global cutoff_value_7
	global cutoff_value_8
	global cutoff_value_9
	global cutoff_value_A
	global cutoff_value_B
	global cutoff_value_C

	global cutoff_value_D
	global cutoff_value_E
	global cutoff_value_F
	global cutoff_value_G
	global cutoff_value_H
	global cutoff_value_I

	global cutoff_value_J
	global cutoff_value_K
	global cutoff_value_L
	global cutoff_value_M
	global cutoff_value_N
	global cutoff_value_O

	### VALUES ###
	cutoff_value_1 = 0.600
	cutoff_value_2 = 0.550
	cutoff_value_3 = 0.500
	cutoff_value_4 = 0.450
	cutoff_value_5 = 0.400
	cutoff_value_6 = 0.350

	cutoff_value_7 = 0.300
	cutoff_value_8 = 0.275
	cutoff_value_9 = 0.250
	cutoff_value_A = 0.225
	cutoff_value_B = 0.200
	cutoff_value_C = 0.175

	cutoff_value_D = 0.150
	cutoff_value_E = 0.140
	cutoff_value_F = 0.130
	cutoff_value_G = 0.120
	cutoff_value_H = 0.110
	cutoff_value_I = 0.100

	cutoff_value_J = 0.075
	cutoff_value_K = 0.050
	cutoff_value_L = 0.025
	cutoff_value_M = 0.020
	cutoff_value_N = 0.010
	cutoff_value_O = 0.000

def Set_CutOff_Values_Type_4():

	### GLOBAL ###
	global cutoff_value_1
	global cutoff_value_2
	global cutoff_value_3
	global cutoff_value_4
	global cutoff_value_5
	global cutoff_value_6

	global cutoff_value_7
	global cutoff_value_8
	global cutoff_value_9
	global cutoff_value_A
	global cutoff_value_B
	global cutoff_value_C

	global cutoff_value_D
	global cutoff_value_E
	global cutoff_value_F
	global cutoff_value_G
	global cutoff_value_H
	global cutoff_value_I

	global cutoff_value_J
	global cutoff_value_K
	global cutoff_value_L
	global cutoff_value_M
	global cutoff_value_N
	global cutoff_value_O

	### VALUES ###
	cutoff_value_1 = 0.300
	cutoff_value_2 = 0.260
	cutoff_value_3 = 0.240
	cutoff_value_4 = 0.220
	cutoff_value_5 = 0.200
	cutoff_value_6 = 0.180

	cutoff_value_7 = 0.160
	cutoff_value_8 = 0.150
	cutoff_value_9 = 0.140
	cutoff_value_A = 0.130
	cutoff_value_B = 0.120
	cutoff_value_C = 0.110

	cutoff_value_D = 0.100
	cutoff_value_E = 0.090
	cutoff_value_F = 0.080
	cutoff_value_G = 0.070
	cutoff_value_H = 0.060
	cutoff_value_I = 0.050

	cutoff_value_J = 0.040
	cutoff_value_K = 0.030
	cutoff_value_L = 0.020
	cutoff_value_M = 0.010
	cutoff_value_N = 0.005
	cutoff_value_O = 0.000

def Set_CutOff_Values_Type_5():

	### GLOBAL ###
	global cutoff_value_1
	global cutoff_value_2
	global cutoff_value_3
	global cutoff_value_4
	global cutoff_value_5
	global cutoff_value_6

	global cutoff_value_7
	global cutoff_value_8
	global cutoff_value_9
	global cutoff_value_A
	global cutoff_value_B
	global cutoff_value_C

	global cutoff_value_D
	global cutoff_value_E
	global cutoff_value_F
	global cutoff_value_G
	global cutoff_value_H
	global cutoff_value_I

	global cutoff_value_J
	global cutoff_value_K
	global cutoff_value_L
	global cutoff_value_M
	global cutoff_value_N
	global cutoff_value_O

	### VALUES ###
	cutoff_value_1 = 0.200
	cutoff_value_2 = 0.190
	cutoff_value_3 = 0.180
	cutoff_value_4 = 0.170
	cutoff_value_5 = 0.160
	cutoff_value_6 = 0.150

	cutoff_value_7 = 0.140
	cutoff_value_8 = 0.130
	cutoff_value_9 = 0.120
	cutoff_value_A = 0.110
	cutoff_value_B = 0.100
	cutoff_value_C = 0.090

	cutoff_value_D = 0.080
	cutoff_value_E = 0.070
	cutoff_value_F = 0.060
	cutoff_value_G = 0.050
	cutoff_value_H = 0.040
	cutoff_value_I = 0.030

	cutoff_value_J = 0.025
	cutoff_value_K = 0.020
	cutoff_value_L = 0.015
	cutoff_value_M = 0.010
	cutoff_value_N = 0.005
	cutoff_value_O = 0.000

def Read_Data_File(in_name, out_name, max_cut, min_cut, dat_cut, frame_file_name):

	### CUTOFF VALUES ###
	global cutoff_value_1
	global cutoff_value_2
	global cutoff_value_3
	global cutoff_value_4
	global cutoff_value_5
	global cutoff_value_6

	global cutoff_value_7
	global cutoff_value_8
	global cutoff_value_9
	global cutoff_value_A
	global cutoff_value_B
	global cutoff_value_C

	global cutoff_value_D
	global cutoff_value_E
	global cutoff_value_F
	global cutoff_value_G
	global cutoff_value_H
	global cutoff_value_I

	global cutoff_value_J
	global cutoff_value_K
	global cutoff_value_L
	global cutoff_value_M
	global cutoff_value_N
	global cutoff_value_O
	#####################

	global max_diff_per_point

	global round_scale

	global print_matrix_file
	global print_pairwise_a
	global print_pairwise_b

	global item_list_sorting

	###################

	print "============================================="
	print "    RUN PARAMETERS:  "
	print " 1. INPUT  FILE:  " + in_name
	print " 2. OUTPUT FILE:  " + out_name
	print " 3. MAX VALUE:    " + str(max_cut)
	print " 4. MIN VALUE:    " + str(min_cut)
	print " 5. DATA CUTOFF:  " + str(dat_cut)
	print " 6. FRAME  LIST:  " + frame_file_name
	print " 7. LIST ORDER:   " + item_list_sorting

	in_file  = open(in_name,  "rb")
	if print_matrix_file == "TRUE":
		out_file1 = open(out_name + '.pairs_all', "wb")
	############# GOOD AND BAD SET ###############
	out_file4 = open(out_name + '.set_uniq', "wb")		; # NON-REDUNDANT ITEM LIST
	out_file4f = open(out_name + '.set_uniq.fixed', "wb")	; # FIXED VALUES
	out_file5 = open(out_name + '.set_dupl', "wb")		; # DUPLICATED IDs
	##############################################
	out_file6 = open(out_name + '.x_log_file', "wb")
	##############################################
	global out_file8
	out_file8 = open(out_name + '.x_tree_clust', "wb")
	##############################################
	###      TWO DIMENSIONAL MATRIX            ###
	##############################################
	out_file1_2d = open(out_name + '.pairs_all.2D_Matrix', "wb")
	out_file1_PWa = open(out_name + '.pairs_all.PW_Val.tab', "wb")
	out_file1_PWb = open(out_name + '.pairs_all.PW_0_1.tab', "wb")
	out_file1_2d_ascii = open(out_name + '.pairs_all.2D_ASCII', "wb")

	out_file6.write("=============================================" + '\n')
	out_file6.write("    RUN PARAMETERS: " + '\n')
	out_file6.write(" 1. INPUT  FILE:  " + in_name + '\n')
	out_file6.write(" 2. OUTPUT FILE:  " + out_name + '\n')
	out_file6.write(" 3. MAX CUTOFF:   " + str(max_cut) + '\n')
	out_file6.write(" 4. MIN CUTOFF:   " + str(min_cut) + '\n')
	out_file6.write(" 5. DATA CUTOFF:  " + str(dat_cut) + '\n')
	out_file6.write(" 6. FRAME  LIST:  " + frame_file_name + '\n')
	out_file6.write(" 7. LIST ORDER:   " + item_list_sorting + '\n')

	time.sleep(2)

	global id_list
	global id_array
	global pairs_array
	global frame_array
	global tree_clust_array
	global graph_depth
	global marker_depth
	global init_len

	id_list = []		; # LIST OF ALL NON-REDUNDANT IDs
	id_array = {}		; # TO CHECK AND CREATE NON-REDUNDANT LIST
	matrix_array = {}	; # PAIRWISE DATA FOR GOOD PAIRS ( <= 0.8 )

	########################### EXAMPLE DEFAULT VALUES 
	matrix_array_1 = {}	; # 0.800  01
	matrix_array_2 = {}	; # 0.600  02
	matrix_array_3 = {}	; # 0.500  03
	matrix_array_4 = {}	; # 0.400  04
	matrix_array_5 = {}	; # 0.300  05
	matrix_array_6 = {}	; # 0.250  06
	matrix_array_7 = {}	; # 0.200  07
	matrix_array_8 = {}	; # 0.180  08
	matrix_array_9 = {}	; # 0.160  09
	matrix_array_A = {}	; # 0.140  10  A
	matrix_array_B = {}	; # 0.120  11  B
	matrix_array_C = {}	; # 0.100  12  C
	matrix_array_D = {}	; # 0.090  13  D
	matrix_array_E = {}	; # 0.080  14  E
	matrix_array_F = {}	; # 0.070  15  F
	matrix_array_G = {}	; # 0.060  16  G
	matrix_array_H = {}	; # 0.050  17  H
	matrix_array_I = {}	; # 0.040  18  I
	matrix_array_J = {}	; # 0.030  19  J
	matrix_array_K = {}	; # 0.020  20  K
	matrix_array_L = {}	; # 0.015  21  L
	matrix_array_M = {}	; # 0.010  22  M
	matrix_array_N = {}	; # 0.005  23  N
	matrix_array_O = {}	; # 0.000  24  O
	##################################################

	pairs_array  = {}	; # ARRAY OF ALL POSITIVE PAIRS IN THE MATRIX
	pair_counter = 0	; # COUNTER OF POSITIVE PAIRS IN THE MATRIX
	sequence_array = {}	; # ARRAY OF ALL SEQUENCES
	sequence_array_bin = {}	; # ARRAY OF ALL SEQUENCES IN BINARY MODE: A or B and -
	frame_array = {}	; # ARRAY OF FRAME MARKERS
	frame_position = {}	; # ARRAY WITH FRAME MARKERS COORDINATES
	bit_sum_array = {}	; # SUMMARY OF BIT SCORES
	rec_sum_array = {}	; # SUMMARY OF REC SCORES
	tree_clust_array = {}	; # SUMMARY OF CLUSTERING ARRAY

	nr_scores_array = {}	; # ARRAY NON-REDUNDANT MARKER SCORES I
	nr_markers_array = {}	; # ARRAY NON-REDUNDANT MARKER SCORES II
	all_scores_array = {}	; # ARRAY WITH ALL MARKER SCORES

	pos_val_array = {}	; # ARRAY WITH POSITIVE VALUES
	neg_val_array = {}	; # ARRAY WITH NEGATIVE VALUES
	nul_val_array = {}	; # ARRAY WITH ZERO VALUES
	mis_val_array = {}	; # ARRAY WITH MISSING VALUES
	all_val_array = {}	; # ARRAY WITH ALL VALUES COUNT
	fix_val_array = {}	; # ARRAY WITH FIXED VALUES

	graph_depth  = {}	; # TYPE OF GRAPH: COMPLETE or CONNECTED
	marker_depth = {}	; # MARKER DEPTH: IS IT CONNECTED DIRECTLY TO ALL OTHER MARKERS IN A GROUP?

	global print_frame
	print_frame = "FALSE"

	### TRY TO READ FRAME MARKERS LIST ###
	try:
		frame_file = open(frame_file_name, "rb")
		print "USING FRAME MARKERS LIST"
		while 1:
			u = frame_file.readline()
			if u == '':
				break
			if '\n' in u:
				u = u[:-1]
			if '\r' in u:
				u = u[:-1]
			u = u.split('\t')
			fm = u[1]
			fl = u[0]
			try:
				fp = u[2]
			except:
				fp = "-"
			frame_array[fm] = fl
			frame_position[fm] = fp
			print_frame = "TRUE"
	except:
		print "DID NOT FIND FRAME MARKERS FILE:  " + frame_file_name
		print "WORKING WITHOUT FRAME MARKERS LIST"
		# continue

	### READ DATA FILE ###
	print "============================================="
	time.sleep(2)
	out_file6.write("=============================================" + '\n')
	n = 0
	d = 0
	l = 1
	init_len = 100000
	duplic_id = []
	while 1:
		t = in_file.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		tl = t.split('\t')
		curr_len = len(tl)
		if tl[0][0] == ";":
			print "============================================="
			print tl
			print "INFO LINE"
			print "============================================="
			out_file4.write(t + '\n')
			out_file4f.write(t + '\n')
			time.sleep(2)
		if l == 1 and curr_len >= 12 and tl[0][0] != ";":
			init_len = curr_len
		if curr_len == init_len and tl[0][0] != ";":
			####
			id = tl[0]
			# if id not in id_list:
			try:
				id_test = id_array[id]
				duplic_id.append(id)
				out_file5.write(t + '\n')
				print '\n'
				print id + "    IS DUPLICATED -=- CHECK DATA FOR DUPLICATION"
				out_file6.write( id + "    IS DUPLICATED -=- CHECK DATA FOR DUPLICATION" + '\n')
				d = d + 1
				time.sleep(2)
			except:
				id_array[id] = 1
				id_list.append(id)

				pos_val_array[id] = 0
				neg_val_array[id] = 0
				nul_val_array[id] = 0
				mis_val_array[id] = 0
				all_val_array[id] = 0
				fix_val_array[id] = 0
				
				out_file4.write(t + '\n')
				out_file4f.write(id + '\t')
				q = 1

				## COLLECTING DATAPOINTS
				while q < init_len:
					data_point = tl[q]
					try:
						point_value = float(data_point)
						test_float = "TRUE"
					except:
						test_float = "FALSE"

					if test_float == "TRUE":
						### DIGITAL DATA ###
						### sequence_array[id,q] = point_value
						if point_value > max_cut:
							point_value = max_cut
							fix_val_array[id] = fix_val_array[id] + 1
						if point_value < min_cut:
							point_value = min_cut
							fix_val_array[id] = fix_val_array[id] + 1
						if point_value == 0:
							nul_val_array[id] = nul_val_array[id] + 1
						if point_value  > 0:
							pos_val_array[id] = pos_val_array[id] + 1
						if point_value  < 0:
							neg_val_array[id] = neg_val_array[id] + 1
						all_val_array[id] = all_val_array[id] + 1
						### ARRAY WITH VALUES PER DATAPOINT
						sequence_array[id,q] = point_value
						point_val_fixed = str(round(point_value,round_scale))

					if test_float == "FALSE":
						### MISSING DATA ###
						sequence_array[id,q] = "-"	; # DASH IF NO DATA
						mis_val_array[id] = mis_val_array[id] + 1
						all_val_array[id] = all_val_array[id] + 1
						point_val_fixed = "-"

					if q < (init_len - 1):
						out_file4f.write(point_val_fixed + '\t')
					if q == (init_len - 1):
						out_file4f.write(point_val_fixed + '\n')

					q = q + 1
				sys.stdout.write(".")
				n = n + 1
			l = l + 1
		if curr_len != init_len and l > 1:
			print '\n'
			print "WRONG NUMBER OF DATA POINTS"
			print "CHECK LINE:  " + `l` + '\n'
			print t
			print "============================================="
			print "TERMINATED.........."
			print "============================================="
			out_file6.write("TERMINATED")
			out_file6.write(".........." + '\n')
			sys.exit()

	print '\n'
	print "============================================="
	print `n` + " UNIQ IDs IN THE SET FOUND"
	print `d` + " IDs ARE DUPLICATED"
	out_file6.write("=============================================" + '\n')
	out_file6.write(`n` + " UNIQ IDs IN THE SET FOUND" + '\n')
	out_file6.write(`d` + " IDs ARE DUPLICATED" + '\n')

	duplic_id.sort()
	print duplic_id
	out_file6.write("=============================================" + '\n')

	print "CONTINUE ANALYSIS WITH " + `n` + " SEQUENCES OUT OF " + `n+d`
	out_file6.write("CONTINUE ANALYSIS WITH " + `n` + " SEQUENCES OUT OF " + `n+d` + '\n')
	print "============================================="
	print "SUCCESS!!!"
	print "============================================="
	out_file6.write("=============================================" + '\n')

	time.sleep(2)

	######### PAIRWISE COMPARISON #########

	## SORTED ID LIST ##
	# id_list.sort()
	## IS IT BETTER NOT TO SORT ? ##
	print ""
	print ""
	print "ORIGINAL LIST: "
	print id_list
	print ""
	time.sleep(2)

	### WHAT TO DO WITH ITEM LIST - SORT - SHUFFLE - AS_IS ###
	if item_list_sorting == "ABC_ORDER":
		id_list.sort()
	if item_list_sorting == "DONT_TOUCH":
		mooba = "booba"
	if item_list_sorting == "SHUFFLE_ME":
		random.seed()
		random.shuffle(id_list)

	print item_list_sorting
	print ""

	print "WORKING LIST: "
	print id_list
	print ""
	time.sleep(2)

	################################

	### SET ZERO FOR MAX DIFF VALUE ###
	max_diff_per_point = 0

	for item in id_list:
		tree_clust_array[item] = []

	### START ###

	dummy_cat = 1
	dummy_len = len(id_list)
	dummy_value = dummy_len*dummy_len
	print `dummy_value` + "  PAIRWISE COMPARISONS HAVE TO BE DONE"

	count_item_a = 0
	for item_a in id_list:
		count_item_a = count_item_a + 1
		print item_a + "   -=-   " + `count_item_a` + "  MAX DIFF IS: " + str(round(max_diff_per_point,6))
		for item_b in id_list:
			p = 1
			diff_score = 0
			data_loss = 0
			data_points = 0

			line_tick_update = math.fmod(dummy_cat, 1000)
			if line_tick_update == 0:
				print `dummy_cat` + "  COMPARISONS DONE OUT OF  " + `dummy_value`

			while p < init_len:

				score_a = sequence_array[item_a,p]
				score_b = sequence_array[item_b,p]

				try:
					test_float_a = float(score_a)
				except:
					test_float_a = "-"
				try:
					test_float_b = float(score_b)
				except:
					test_float_b = "-"

				if test_float_a != "-" and test_float_b != "-":
					### PRINT DATA ###
					diff_score_ab = score_a - score_b
					diff_score = diff_score + abs(diff_score_ab)
					### print ""
					data_points = data_points + 1
				if test_float_a == "-" or test_float_b == "-":
					### print ".....MISSING DATA....."
					### print ""
					data_loss = data_loss + 1
				#####################################

				p = p + 1

			dummy_cat = dummy_cat + 1
			### DATA POINT TEST ###
			if data_points + data_loss != init_len - 1:
				print "...SOMETHING IS WRONG..."
				sys.exit()

			### ANALYZE ONLY VALUABLE CRAP ###
			data_fr = data_points*1.00/(data_points + data_loss)
			# data_fr = round(data_fr,round_scale)
			data_fr = round(data_fr,4)
			data_fr_str = str(data_fr)
			if data_points >= dat_cut:

				### MATRIX DATA  ###

				diff_per_point = diff_score/data_points
				diff_per_point = round(diff_per_point,round_scale)
				diff_per_point_str = str(diff_per_point)

				if diff_per_point > max_diff_per_point:
					max_diff_per_point = diff_per_point
					print "FOUND NEW MAX DIFF VALUE:  " + str(round(max_diff_per_point,6)) + "  FOR PAIR: " + item_a + " -=- " + item_b

				diff_score = round(diff_score,round_scale)
				diff_score_str = str(diff_score)

				pair_counter = pair_counter + 1
				current_pair = [item_a, item_b]
				pairs_array[pair_counter] = current_pair
				matrix_values = [diff_score, data_points, diff_per_point]
				matrix_array[item_a,item_b] = matrix_values
				### PAIRWISE MATRIX ###
				if print_matrix_file == "TRUE":
					out_file1.write(item_a + '\t' + item_b + '\t' + diff_score_str + '\t' + \
						"*[DFR>*" + '\t' + data_fr_str + '\t' + \
						"*[VPP>*" + '\t' + diff_per_point_str + '\t' + \
						"*[PCC>*" + '\t' + `data_points` + '\t' + `data_loss` + '\t' + \
						`data_points + data_loss` + '\n')
					### PW Matrix File ###
					if print_pairwise_a == "TRUE":
						out_file1_PWa.write(item_a + '\t' + item_b + '\t' + diff_per_point_str + '\n')
				########################################

			if data_points < dat_cut:
				print "DATA POINTS ARE LESS THAN CUTOFF FOR PAIR: " + item_a + " " + item_b + " " + \
										`data_points` + "   " + data_fr_str
			########################################

		#####################
		sys.stdout.write(".")
		#####################
	print ""
	print "============================================="
	print "GLOBAL MATRIX CREATED"
	print "============================================="

	print ""
	print "============================================="
	print "MAX DIFF VALUE:  " + str(round(max_diff_per_point,6))
	print "============================================="

	out_file6.write("MAX DIFF VALUE:  " + str(round(max_diff_per_point,6)) + '\n')

	### MATRIX PER ITERATION ###

	for item_a1 in id_list:
		for item_b2 in id_list:
			try:
				cur_diff_a1b2 = float(matrix_array[item_a1,item_b2][2])
				cur_diff_fr_a1b2 = cur_diff_a1b2/max_diff_per_point
				matrix_values_a1b2 = matrix_array[item_a1,item_b2]
				mooba = "OK"
			except:
				mooba = "booba"
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_1:
				matrix_array_1[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_2:
				matrix_array_2[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_3:
				matrix_array_3[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_4:
				matrix_array_4[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_5:
				matrix_array_5[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_6:
				matrix_array_6[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_7:
				matrix_array_7[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_8:
				matrix_array_8[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_9:
				matrix_array_9[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_A:
				matrix_array_A[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_B:
				matrix_array_B[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_C:
				matrix_array_C[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_D:
				matrix_array_D[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_E:
				matrix_array_E[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_F:
				matrix_array_F[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_G:
				matrix_array_G[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_H:
				matrix_array_H[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_I:
				matrix_array_I[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_J:
				matrix_array_J[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_K:
				matrix_array_K[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_L:
				matrix_array_L[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_M:
				matrix_array_M[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_N:
				matrix_array_N[item_a1,item_b2] = matrix_values_a1b2 
			if mooba == "OK" and cur_diff_fr_a1b2 <= cutoff_value_O:
				matrix_array_O[item_a1,item_b2] = matrix_values_a1b2 


	#### 24 ROUNDS OF CLUSTERING ####

	print ""
	print "STARTING GROUP SEARCH - 24 CLUSTERING ROUNDS"
	print "============================================="

	global out_file14

	out_file14 = open(out_name + '.group_info' + "_Summary", "wb")

	time.sleep(2)

	Seqs_Clustering(cutoff_value_1, "01", matrix_array_1)
	print "ROUND 1 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 01)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 01)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_2, "02", matrix_array_2)
	print "ROUND 2 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 02)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 02)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_3, "03", matrix_array_3)
	print "ROUND 3 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 03)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 03)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_4, "04", matrix_array_4)
	print "ROUND 4 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 04)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 04)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_5, "05", matrix_array_5)
	print "ROUND 5 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 05)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 05)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_6, "06", matrix_array_6)
	print "ROUND 6 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 06)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 06)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_7, "07", matrix_array_7)
	print "ROUND 7 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 07)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 07)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_8, "08", matrix_array_8)
	print "ROUND 8 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 08)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 08)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_9, "09", matrix_array_9)
	print "ROUND 9 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 09)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 09)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_A, "10", matrix_array_A)
	print "ROUND 10 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 10)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 10)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_B, "11", matrix_array_B)
	print "ROUND 11 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 11)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 11)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_C, "12", matrix_array_C)
	print "ROUND 12 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 12)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 12)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_D, "13", matrix_array_D)
	print "ROUND 13 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 13)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 13)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_E, "14", matrix_array_E)
	print "ROUND 14 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 14)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 14)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_F, "15", matrix_array_F)
	print "ROUND 15 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 15)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 15)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_G, "16", matrix_array_G)
	print "ROUND 16 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 16)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 16)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_H, "17", matrix_array_H)
	print "ROUND 17 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 17)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 17)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_I, "18", matrix_array_I)
	print "ROUND 18 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 18)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 18)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_J, "19", matrix_array_J)
	print "ROUND 19 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 19)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 19)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_K, "20", matrix_array_K)
	print "ROUND 20 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 20)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 20)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_L, "21", matrix_array_L)
	print "ROUND 21 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 21)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 21)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_M, "22", matrix_array_M)
	print "ROUND 22 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 22)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 22)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_N, "23", matrix_array_N)
	print "ROUND 23 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 23)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 23)" + '\n')

	time.sleep(2)

	Seqs_Clustering(cutoff_value_O, "24", matrix_array_O)
	print "ROUND 24 DONE"
	print `group_count` + "  GROUPS WERE FOUND (STEP 24)"
	out_file6.write(`group_count` + "  GROUPS WERE FOUND (STEP 24)" + '\n')

	##################

	time.sleep(2)

	out_file6.write("=============================================" + '\n')

	### GROUP LIST ###
	f = 0
	for item in id_list:
		### TRY FRAME MARKER ID ###
		if print_frame == "TRUE":
			try:
				frame_group = frame_array[item]
			except:
				frame_group = "-"
			try:
				frame_coords = frame_position[item]
			except:
				frame_coords = "-"
		if print_frame == "FALSE":
			frame_group  = "-"
			frame_coords = "-"
		tree_clust_array[item].append("***")
		tree_clust_array[item].append(frame_group)
		tree_clust_array[item].append(frame_coords)
		tree_clust_array[item].append("***")
		tree_clust_array[item].append(str(f))
		tree_clust_array[item].append("***")
		tree_clust_array[item].append("LG")
		tree_clust_array[item].append(item)
		tree_clust_array[item].append("*")
		# tree_clust_array[item].append("=1=")
		# tree_clust_array[item].append("=2=")
		# tree_clust_array[item].append("=3=")
		# tree_clust_array[item].append("=4=")
		# tree_clust_array[item].append("=5=")
		tree_clust_array[item].append(str(neg_val_array[item]))
		tree_clust_array[item].append(str(nul_val_array[item]))
		tree_clust_array[item].append(str(pos_val_array[item]))
		tree_clust_array[item].append("*")
		tree_clust_array[item].append(str(all_val_array[item]))
		tree_clust_array[item].append(str(mis_val_array[item]))
		tree_clust_array[item].append(str(fix_val_array[item]))
		tree_clust_array[item].append("*")
		f = f + 1
		print tree_clust_array[item]

	print ""
	print "============================================="	
	print "            FINAL STEP                       "
	print "       PROCESSING MEGALIST                   "
	print "         DENDRO SORTING                      "
	print "============================================="
	print ""

	time.sleep(2)

	### MEGA LIST ###
	mega_list = []
	for item in id_list:
		mega_list.append(tree_clust_array[item])
	# print mega_list
	mega_list.sort()
	f = 0
	for item in mega_list:
		current_list = []
		item.append(str(f))
		for subitem in item:
			subitem = str(subitem)
			current_list.append(subitem)
		current_list = "\t".join(current_list)
		print current_list
		out_file8.write(current_list + '\n')
		f = f + 1

	### TWO DIMENSIONAL MATRIX ###

	out_file1_2d.write(";" + '\t')

	mega_length = len(mega_list)

	### FIRST ROW ###
	counter = 0
	for item in mega_list:
		current_id = item[33]
		out_file1_2d.write(current_id)
		counter = counter + 1
		if counter < mega_length:
			out_file1_2d.write('\t')
		if counter == mega_length:
			out_file1_2d.write('\n')

	### PAIRWISE DATA ROWS ###
	for item in mega_list:
		current_id = item[33]
		out_file1_2d.write(current_id + '\t')
		out_file1_2d_ascii.write(current_id + '\t')
		counter = 0
		for other_item in mega_list:
			counter = counter + 1
			other_id = other_item[33]
			try:
				cur_diff = float(matrix_array[current_id,other_id][2])
				cur_diff_string = str(round(cur_diff,round_scale))
				cur_diff_fr = cur_diff/max_diff_per_point
				cur_value = str(round(cur_diff_fr,round_scale))
				cur_value_ascii = (cur_diff/max_diff_per_point)*10
				cur_value_ascii = round(cur_value_ascii,0)
				cur_value_ascii = int(cur_value_ascii)
				if cur_value_ascii >= 10:
					cur_value_ascii = 9
				cur_value_ascii = str(cur_value_ascii)
			except:
				cur_value = "-"
				cur_value_ascii = "-"
			out_file1_2d.write(cur_value)
			out_file1_2d_ascii.write(cur_value_ascii)
			if cur_value != "-" and print_pairwise_b == "TRUE":
				out_file1_PWb.write(current_id + '\t' + other_id + '\t' + cur_value + '\t' + "***" + '\t' + cur_diff_string + '\n')
			if counter < mega_length:
				out_file1_2d.write('\t')
			if counter == mega_length:
				out_file1_2d.write('\n')
				out_file1_2d_ascii.write('\n')


	##########################################

	print ""
	print "============================================="	
	print "            ANALYSIS DONE                    "
	print "       HAPPY ENDING IS OPTIONAL              "
	print "============================================="
	print ""

	in_file.close()
	if print_matrix_file == "TRUE":
		out_file1.close()
	out_file4.close()
	out_file4f.close
	out_file5.close()
	out_file6.close()
	out_file8.close()
	out_file14.close()
	out_file1_2d.close()
	out_file1_PWa.close()
	out_file1_PWb.close()
	out_file1_2d_ascii.close()

#### POSITIVE CLUSTERING ####

def Seqs_Clustering(cutoff_value, x_counter, matrix_array_current):

	global max_diff_per_point

	global id_list
	global id_array
	global adj_array
	global already_done
	global group_count
	global dfs_counter
	global working_group
	global group_depth
	global node_count
	global node_array
	global sequence_array
	global sequence_array_bin
	global tree_clust_array

	global pairs_array
	global frame_array
	global print_frame
	global round_scale

	global graph_depth
	global marker_depth

	global out_file14
	global out_file8

	global init_len

	global print_temp_clust_files

	node_array = {}
	adj_array    = {}
	already_done = []
	group_count = 0
	pair_matrix_count = 0
	pair_matrix_array = {}

	print ""
	print "ALREADY DONE"
	print already_done
	print ""
	print ""
	print "MATRIX ARRAY"
	print len(matrix_array_current)
	print ""
	time.sleep(2)

	out_file10 = open(out_name + '.matrix' + "_" + x_counter, "wb")
	if print_temp_clust_files == "TRUE":
		out_file12 = open(out_name + '.adj_list' + "_" + x_counter, "wb")
		out_file13 = open(out_name + '.group_info' + "_" + x_counter, "wb")

	out_file14.write("### CLUSTERING RESULTS ITERATION " + x_counter + " CUTOFF: " + str(cutoff_value) + " ###" + '\n')

	########################################
	###          CLUSTERING              ###
	########################################

	# print id_list
	print "------------------------------"
	print "FOUND " + `len(id_list)` + " UNIQ IDs"
	print "------------------------------"

	r = 0

	### CREATE NON-REDUNDANT MATRIX ###

	for key in pairs_array:
		id_a = pairs_array[key][0]
		id_b = pairs_array[key][1]
		try:
			cur_data1 = int(matrix_array_current[id_a,id_b][1])
			cur_diff1 = float(matrix_array_current[id_a,id_b][2])
			cur_diff1_fr = cur_diff1/max_diff_per_point
			query1 = 1
			print key
		except:
			# print "ALREADY PROCESSED"
			sys.stdout.write(".")
			query1 = 0

		try:
			cur_data2 = int(matrix_array_current[id_b,id_a][1])
			cur_diff2 = float(matrix_array_current[id_b,id_a][2])
			cur_diff2_fr = cur_diff2/max_diff_per_point
			r = r + 1
			# print `r` + " REVERSE PAIR FOUND"
			query2 = 1
		except:
			# print "NO REVERSE PAIR FOUND"
			query2 = 0

		if query1 == 1 and query2 == 0 and cur_diff1_fr <= cutoff_value and id_a != id_b:
			### STRING CONVERSION ###
			cur_diff1_str = str(round(cur_diff1,round_scale))
			cur_diff1_fr_str = str(round(cur_diff1_fr,round_scale))
			#########################
			out_file10.write(id_a + '\t' + id_b + '\t' + cur_diff1_str + '\t' + \
				`cur_data1` + '\t' + cur_diff1_fr_str + '\n')
			pair_matrix_count = pair_matrix_count + 1
			current_matrix_pair = [id_a, id_b]
			pair_matrix_array[id_a,id_b] = current_matrix_pair
			### NEED TO UNSET ARRAY
			try:
				del matrix_array_current[id_a,id_b]
			except:
				# print "ALREADY REMOVED"
				sys.stdout.write(".")

		if query1 == 1 and query2 == 1 and id_a != id_b:
			if cur_diff1_fr <= cur_diff2_fr:
				# print "CASE 1"
				if cur_diff1_fr <= cutoff_value:
					### STRING CONVERSION ###
					cur_diff1_str = str(round(cur_diff1,round_scale))
					cur_diff1_fr_str = str(round(cur_diff1_fr,round_scale))
					#########################
					out_file10.write(id_a + '\t' + id_b + '\t' + cur_diff1_str + '\t' + \
						`cur_data1` + '\t' + cur_diff1_fr_str + '\n')
					pair_matrix_count = pair_matrix_count + 1
					current_matrix_pair = [id_a, id_b]
					pair_matrix_array[id_a,id_b] = current_matrix_pair
					### NEED TO UNSET ARRAY
					try:
						del matrix_array_current[id_a,id_b]
						del matrix_array_current[id_b,id_a]
					except:
						# print "ALREADY REMOVED"
						sys.stdout.write(".")
			if cur_diff1_fr > cur_diff2_fr:
				# print "CASE 2"
				if cur_diff2_fr <= cutoff_value:
					### STRING CONVERSION ###
					cur_diff2_str = str(round(cur_diff2,round_scale))
					cur_diff2_fr_str = str(round(cur_diff2_fr,round_scale))
					#########################
					out_file10.write(id_a + '\t' + id_b + '\t' + cur_diff2_str + '\t' + \
						`cur_data2` + '\t' + cur_diff2_fr_str + '\n')
					pair_matrix_count = pair_matrix_count + 1
					current_matrix_pair = [id_b, id_a]
					pair_matrix_array[id_b,id_a] = current_matrix_pair
					### NEED TO UNSET ARRAY
					try:
						del matrix_array_current[id_a,id_b]
						del matrix_array_current[id_b,id_a]
					except:
						# print "ALREADY REMOVED"
						sys.stdout.write(".")

	print "-------------------------------------"
	print `pair_matrix_count` + " PAIRS IN REDUNDANT MATRIX"
	print "-------------------------------------"
	print "BEGIN CLUSTERING"

	time.sleep(2)

	item_count = 0
	id_list.sort()
	### CREATE ADJACENCY LIST ###
	for item in id_list:
		item_count = item_count + 1
		print `item_count` + '\t' + item

		item_list = [item]
		for key in pair_matrix_array:
			id_a = pair_matrix_array[key][0]
			id_b = pair_matrix_array[key][1]
			if id_a == item:
				item_list.append(id_b)
			if id_b == item:
				item_list.append(id_a)
		adj_array[item] = item_list
		item_string = " ".join(item_list)
		if print_temp_clust_files == "TRUE":
			out_file12.write(item_string + '\n')

	print "GROUP ANALYSIS"
	time.sleep(2)

	### GROUP ANALYSIS ###
	node_count = 0
	for item in id_list:
		# if item not in already_done:
		try:
			node_test = node_array[item]
		except:
			node_array[item] = 1
			group_count = group_count + 1
			already_done.append(item)
			node_count = node_count + 1
			working_group = [item]
			current_adj_list = adj_array[item]
			current_adj_len = len(current_adj_list)

			q = 0
			while q <= (current_adj_len - 1):
				current_adj_item = current_adj_list[q]
				# print `q` + '\t' + current_adj_item
				if current_adj_item in already_done:
					go_to_dfs = 0
				# if current_adj_item not in already_done:
				try:
					node_test = node_array[current_adj_item]
				except:
					node_array[current_adj_item] = 1
					already_done.append(current_adj_item)
					node_count = node_count + 1
					if current_adj_item not in working_group:
						working_group.append(current_adj_item)
					go_to_dfs = 1
					dfs_counter = 0
					# print 'Processing Group:  ' + `group_count`
					### GRAPH SEARCH ###
					DFS_procedure(current_adj_item)
				q = q + 1
			# if item not in already_done:
			#	already_done.append(item)
			working_group.sort()
			# print working_group
			print 'Processing Group:  ' + `group_count`
			print 'Number of processed nodes:  ' + `node_count`
			i = 0
			group_suffix = str(group_count)
			suffix_len = len(group_suffix)
			if suffix_len < 5:
				if suffix_len == 1:
					group_suffix = "0000" + group_suffix
				if suffix_len == 2:
					group_suffix = "000"  + group_suffix
				if suffix_len == 3:
					group_suffix = "00"   + group_suffix
				if suffix_len == 4:
					group_suffix = "0"    + group_suffix
				if suffix_len == 5:
					group_suffix = group_suffix
			## GRAPH TYPE: CONNECTED OR COMPLETE
			## NODE  TYPE: SATURATED OR DILUTED
			graph_type = "COMPLETE_GRAPH_" + group_suffix
			graph_depth[group_count] = graph_type
			for egg in working_group:
				current_adj_list7 = adj_array[egg]
				current_adj_len7 = len(current_adj_list7)
				current_group_len7 = len(working_group)
				# marker_depth[egg] = "SATURATED_NODE"
				marker_depth[egg] = "UNDEFINED_NODE"
				if current_adj_len7 == current_group_len7:
					marker_depth[egg] = "SATURATED_NODE"
				if current_adj_len7 != current_group_len7:
					graph_type = "LINKED___GROUP_" + group_suffix
					graph_depth[group_count] = graph_type
					marker_depth[egg] = "DILUTED___NODE"
				if current_group_len7 == 1:
					graph_type = "SINGLE____NODE_" + group_suffix
					graph_depth[group_count] = graph_type
					# marker_depth[egg] = "SINGLE____NODE"
			for element in working_group:
				current_adj_list1 = adj_array[element]
				current_adj_len1 = len(current_adj_list1)
				current_group_len1 = len(working_group)
				### TRY FRAME MARKER ID ###
				if print_frame == "TRUE":
					try:
						frame_marker = frame_array[element]
						frame_marker = ' ' + frame_marker + '_LinkageGroup '
					except:
						frame_marker = " __LinkageGroup "
				if print_frame == "FALSE":
					frame_marker = "_NONE_"
				###########################
				# tree_clust_array[element].append(str(group_count))
				tree_clust_array[element].append(group_count)
				if x_counter == "24":
					tree_clust_array[element].append(graph_depth[group_count])
					tree_clust_array[element].append(marker_depth[element])
				print tree_clust_array[element]
				###########################
				### WRITE DATA TO OUTPUT GROUP-INFO FILES ###
				if i == 0:
					if print_temp_clust_files == "TRUE":
						out_file13.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
							+ `current_group_len1` + '\t' + `group_count` + '\t' + '*****' \
							+ '\t' + frame_marker + '\t' + graph_depth[group_count] + '\t' \
							+ marker_depth[element] + '\n')
					out_file14.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '*****' \
						+ '\t' + frame_marker + '\t' + "*" + x_counter + "*" + '\t' \
						+ graph_depth[group_count] + '\t' + marker_depth[element] + '\n')
				if i != 0:
					if print_temp_clust_files == "TRUE":
						out_file13.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
							+ `current_group_len1` + '\t' + `group_count` + '\t' + '-----' \
							+ '\t' + frame_marker + '\t' + graph_depth[group_count] + '\t' \
							+ marker_depth[element] + '\n')
					out_file14.write(element + '\t' + `(current_adj_len1 - 1)` + '\t' \
						+ `current_group_len1` + '\t' + `group_count` + '\t' + '-----' \
						+ '\t' + frame_marker + '\t' + "*" + x_counter + "*" + '\t' \
						+ graph_depth[group_count] + '\t' + marker_depth[element] + '\n')
				i = i + 1
				###########################
	print '======================'
	# print already_done
	print len(already_done)
	print `group_count` + '\t' + 'GROUPS FOUND'
	print '======================'

	#######################################

	out_file10.close()
	if print_temp_clust_files == "TRUE":
		out_file12.close()
		out_file13.close()

#############################
### GRAPH SEARCH ###
def DFS_procedure(current_adj_item):

	global adj_array
	global already_done
	global group_count
	global dfs_counter
	global working_group
	global group_depth
	global node_count
	global node_array
	global print_matrix_file

	print "DFS" + '\t' + `dfs_counter`

	for key in adj_array:
		current_adj_list = adj_array[key]
		current_adj_len = len(current_adj_list)
		if current_adj_item in current_adj_list:
			q = 0
			while q <= (current_adj_len - 1):
				current_adj_item1 = current_adj_list[q]
				# print `q` + '\t' + current_adj_item1
				if current_adj_item1 in already_done:
					go_to_dfs = 0
				# if current_adj_item1 not in already_done:
				try:
					node_test = node_array[current_adj_item1]
				except:
					node_array[current_adj_item1] = 1
					already_done.append(current_adj_item1)
					node_count = node_count + 1
					if current_adj_item1 not in working_group:
						working_group.append(current_adj_item1)
					go_to_dfs = 1
					DFS_procedure(current_adj_item1)
					### HOW MANY TIMES LOOP INSIDE ITSELF ###
					dfs_counter = dfs_counter + 1
					print "DFS" + '\t' + `dfs_counter`
					group_depth = dfs_counter
				q = q + 1

	######################################################################

import math
import re
import sys
import string
import time
import random

global round_scale
global print_matrix_file
global print_pairwise_a
global print_pairwise_b
global print_temp_clust_files

# round_scale = 2
# round_scale = 4
round_scale = 6
# round_scale = 8

### WHAT TO DO WITH INITIAL ITEM LIST ###
###  1. RANDOM SHUFFLE - DEFAULT      ### 
###  2. NOTHING - ORDERED AS IS       ###
###  3. SORT BY ITEM ID IN ABC ORDER  ###
global item_list_sorting
# item_list_sorting = "SHUFFLE_ME"
# item_list_sorting = "DONT_TOUCH"
item_list_sorting = "ABC_ORDER"

if __name__ == "__main__":
	if len(sys.argv) <= 6 or len(sys.argv) > 8:
		print "                                                                     "
		print "   PROGRAM  USAGE:                                                   "
		print "   MAD MAPPER TAKES  6  ARGUMENTS/OPTIONS IN THE FOLLOWING ORDER:    "
		print "  (1)input_file[LOC_DATA/MARKER SCORES]        (2)output_file[NAME]  "
		print "  (3)max_value   (4)min_value   (5)paired_datapoints   (6)opt_frame  "
		print "                                                                     "
		input = raw_input("   TYPE \"HELP\" FOR HELP [ \"EXIT\" TO EXIT ] : ")
		if input == "HELP" or input == "help":
			print "                                                                             "
			print "              MAD MAPPER ARGUMENTS/OPTIONS - BRIEF EXPLANATION:              "
			print "                                                                             "
			print "       INPUT/OUTPUT FILES:                                                   "
			print "   [1] - Input File Name  (file with numerical values per item)              "
			print "   [2] - Output File Name (master name [prefix] for 80 or so output files)   "
			print "                                                                             "
			print "   [3] - MAX allowed value                                                   "
			print "   [4] - MIN allowed value                                                   "
			print "   [5] - min number of paired datapoints                                     "
			print "                                                                             "
			print "   [6] - optional list of FRAME items                                        "
			print "                                                                             "
			print "                                                                             "
			sys.exit()
		if input != "HELP" and input != "help":
			sys.exit()
		sys.exit()
	if len(sys.argv) == 7:
		in_name  = sys.argv[1]
		out_name = sys.argv[2]
		max_cut  = float(sys.argv[3])
		min_cut  = float(sys.argv[4])
		dat_cut  = int(sys.argv[5])
		frame_file_name = sys.argv[6]

		sys.setrecursionlimit(100000)

		# Set_CutOff_Values_Type_1()
		# Set_CutOff_Values_Type_2()
		# Set_CutOff_Values_Type_3()
		Set_CutOff_Values_Type_4()
		# Set_CutOff_Values_Type_5()

		# print_matrix_file = "TRUE"
		print_matrix_file = "FALSE"

		# print_temp_clust_files = "TRUE"
		print_temp_clust_files = "FALSE"

		# print_pairwise_a = "TRUE"
		print_pairwise_a = "FALSE"

		print_pairwise_b = "TRUE"
		# print_pairwise_b = "FALSE"

		Read_Data_File(in_name, out_name, max_cut, min_cut, dat_cut, frame_file_name)

#### THE END ####
