# SCRIPT: MC simulation
# MASTER STUDENT: Alberto Fontebasso
# PROJECT: Master Thesis
# SUPERVISOR: PD Dr. Meyer-Massetti
# LAB: IFD

import os

folder_count = 100
sim_count = 1000

nx = 272
ny = 256
l_1 = 775.2
l_2 = 729.6
b_c = 'h'
h_1 = 64.3663
h_2 = 60.4080
tol = 1e-12 

for k in range(folder_count):
	for i in range(sim_count):
		if (i < 10):
			commandline = f"printf '{nx} {ny} \n{l_1} {l_2} \n{b_c} \n{h_1} \n{h_2} \nConductivity_fields_2/Realization_{k}/Realization_{i} \n1 \nResults_2/Realization_{k}/Result_00{i}.txt \n{tol} \n' | ./darcy_amg"
			os.system(commandline)
		elif (i >= 10 and i < 100):
			commandline = f"printf '{nx} {ny} \n{l_1} {l_2} \n{b_c} \n{h_1} \n{h_2} \nConductivity_fields_2/Realization_{k}/Realization_{i} \n1 \nResults_2/Realization_{k}/Result_0{i}.txt \n{tol} \n' | ./darcy_amg"
			os.system(commandline)
		elif (i >= 100):
			commandline = f"printf '{nx} {ny} \n{l_1} {l_2} \n{b_c} \n{h_1} \n{h_2} \nConductivity_fields_2/Realization_{k}/Realization_{i} \n1 \nResults_2/Realization_{k}/Result_{i}.txt \n{tol} \n' | ./darcy_amg"
			os.system(commandline)
