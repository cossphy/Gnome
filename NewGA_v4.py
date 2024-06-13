import random
import os
import itertools
import re
import shutil
import subprocess
import time
import math
import graphviz
import concurrent.futures
from statistics import mean
import matplotlib.pyplot as plt
# from multiprocessing 
from tqdm import trange



# 1st step is creating DataBase:
#===================================#
def ReadUnfilteredDataBase(file_name):
    # Folder and file information
    folder_name = "Targets"
    # count = 1
    targets_file_path = os.path.join(folder_name, file_name)
    with open(targets_file_path, 'r', encoding="utf8") as file:
        lines = file.readlines()
        Data = []
        Row = []
        CountOfBase = 0
        RowsPerBase = 0
        # NumBases = 0
        for line in lines:
            line = line.strip()
            
            if line.isdigit():  # Check for a blank line
                RowsPerBase = int(line) + 2
                # print(helper)
            row = line.split()
            Row.append(row)
            CountOfBase += 1
            if CountOfBase == RowsPerBase:
                Data.append(Row)
                Row = []
                RowsPerBase = 0
                CountOfBase = 0
    return Data

def process_input_file(OldDataList, NewNameDataBase):
    with open(os.path.join('Targets', NewNameDataBase), 'w') as Output:
        for x in range(len(OldDataList)):
            for i in range(len(OldDataList[x])):
                if i == 0:
                    Output.write(f"{OldDataList[x][i][0]}\n")  
                else:  
                    # Output.write(f"{OldDataList[0][i]}\n")
                    if "Lattice=" in OldDataList[x][i][0]:
                        for j in range(len(OldDataList[x][i])):
                            if j == 0:
                                Output.write(f"{OldDataList[x][i][j][9:]} ")
                            elif j == 8:
                                Output.write(f"{OldDataList[x][i][j][:-1]} ")
                            elif j == 9:
                                Output.write(f"\n{OldDataList[x][i][j][12:]}\n")
                            else:
                                Output.write(f"{OldDataList[x][i][j]} ")
                    elif len(OldDataList[x][i]) == 1:
                        Output.write(f"\n{OldDataList[x][i][0]}\n")
                    else:
                        for k in range(len(OldDataList[x][i])):
                            if k == 6:
                                Output.write(f"{OldDataList[x][i][k]}\n")
                            else:
                                Output.write(f"{OldDataList[x][i][k]} ")
            Output.write("\n")        

NameofModifiedDataBase = "O_MODC_DASKALOS.xyz"
NameofUnfilteredDataBase = "database_Pd_O_20240318_adj.xyz"
unfilteredData = ReadUnfilteredDataBase(NameofUnfilteredDataBase)
process_input_file(unfilteredData, NameofModifiedDataBase)



# 2nd step is deleting all the past files used in previous runs
#===================================#
def delete_files_starting_with(prefixes):
    # Get the current working directory
    cwd = os.getcwd()
    # Get the list of files and folders in the current working directory
    files_and_folders = os.listdir(cwd)
    
    for item in files_and_folders:
        if item.startswith(prefixes):
            # Construct the full path to the file/folder
            full_path = os.path.join(cwd, item)
            shutil.rmtree(full_path)
            # print(f"Deleted folder and its contents: {full_path}")

with concurrent.futures.ThreadPoolExecutor() as executor0:
    executor0.map(delete_files_starting_with,['CalcFold_', 'Generation_'])


# 3rd step is reading the Input.txt and parameters Tersoff file.
#===================================#
def readInput():
    with open("input.txt") as Genetic:
        for line in Genetic:
            if re.search(r'BestFrac\s*=\s*([\d.]+)', line):
                BestFrac = float(re.search(r'BestFrac\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'FracGene\s*=\s*([\d.]+)', line):
                FracGene = float(re.search(r'FracGene\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'FracMut\s*=\s*([\d.]+)', line):
                FracMut = float(re.search(r'FracMut\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'FracRand\s*=\s*([\d.]+)', line):
                FracRand = float(re.search(r'FracRand\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'KeepBest\s*=\s*([\d.]+)', line):
                KeepBest = float(re.search(r'KeepBest\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'PopulationSize\s*=\s*([\d.]+)', line):
                PopulationSize = int(re.search(r'PopulationSize\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'MaxMutPerc\s*=\s*([\d.]+)', line):
                MaxMutPerc = int(re.search(r'MaxMutPerc\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'Generations\s*=\s*([\d.]+)', line):
                Generations = int(re.search(r'Generations\s*=\s*([\d.]+)', line).group(1))
            elif re.search(r'PotentialType\s*=\s*(\w+)', line):
                PotentialType = re.search(r'PotentialType\s*=\s*(\w+)', line).group(1)
            elif re.search(r'Constraints\s*=\s*(\w+)', line):
                Constraints = re.search(r'Constraints\s*=\s*(\w+)', line).group(1)
    return BestFrac, FracGene, FracMut, FracRand, KeepBest, PopulationSize, MaxMutPerc, Generations, PotentialType, Constraints

BestFrac, FracGene, FracMut, FracRand, KeepBest, PopulationSize, MaxMutPerc, Generations, PotentialType, Constraints = readInput()

# print(PotentialType)
# print(Constraints)

# import sys
# sys.exit()

total_fraction = FracGene + FracMut + FracRand + KeepBest
if total_fraction != 1:
    FracGene /= total_fraction
    FracMut /= total_fraction
    FracRand /= total_fraction
    KeepBest /= total_fraction

MUT_PERC = MaxMutPerc/100

def readParametersTersoff():
    with open('parameters-TERSOFF-MOD.txt', 'r') as f_param:
        lines = f_param.readlines()
        lower = []
        upper = [] 
        for line in lines[2:10]:
            line = line.split()
            lower.append(line)
        for line in lines[12:]:
            line = line.split()
            upper.append(line)

        combined_array = []

        for row1, row2 in zip(lower, upper):
            combined_row = [row1[:3]]  # First 3 elements from the first array
            combined_row.extend([[elem1, elem2] for elem1, elem2 in zip(row1[3:], row2[3:])])  # Combine elements from both arrays starting from index 3
            combined_array.append(combined_row)

    return combined_array
        
Matrix = readParametersTersoff()
# print(Matrix)
# import sys
# sys.exit()

def Folder_check(folder_path):
    # Create the folder if it doesn't exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    # Clear existing files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")

# 4th step is to create the first generation of potentials
#===================================#
if Constraints == "OFF":
    def FirstGenCreator(Matrix ,PopulationSize):
        Folder_check("Generation_1")
        for p in range(1, PopulationSize+1):
            file_name = os.path.join("Generation_1", f"001-{p:03d}.tersoff.modc")
            with open(file_name, 'a') as file:
                for row in range(len(Matrix)):
                    file.write(f"{Matrix[row][0][0]}\t{Matrix[row][0][1]}\t{Matrix[row][0][2]}\t")
                    for column in range(1, 19):
                        file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                    file.write("\n\n")
elif Constraints == "ON" and PotentialType == "TER":
    def FirstGenCreator(Matrix ,PopulationSize):
        Folder_check("Generation_1")
        for p in range(1, PopulationSize+1):
            file_name = os.path.join("Generation_1", f"001-{p:03d}.tersoff.modc")
            with open(file_name, 'a') as file:
                # 3 body parameters
                # gamma
                gamma_yellow = random.uniform(float(Matrix[0][2][0]), float(Matrix[0][2][1]))
                gamma_green = random.uniform(float(Matrix[2][2][0]), float(Matrix[2][2][1]))
                gamma_red = random.uniform(float(Matrix[6][2][0]), float(Matrix[6][2][1]))
                # lambda3
                lambda3_yellow = random.uniform(float(Matrix[0][3][0]), float(Matrix[0][3][1]))
                lambda3_green = random.uniform(float(Matrix[2][3][0]), float(Matrix[2][3][1]))
                lambda3_red = random.uniform(float(Matrix[6][3][0]), float(Matrix[6][3][1]))
                # c
                c_yellow = random.uniform(float(Matrix[0][4][0]), float(Matrix[0][4][1]))
                c_green = random.uniform(float(Matrix[2][4][0]), float(Matrix[2][4][1]))
                c_red = random.uniform(float(Matrix[6][4][0]), float(Matrix[6][4][1]))
                # d
                d_yellow = random.uniform(float(Matrix[0][5][0]), float(Matrix[0][5][1]))
                d_green = random.uniform(float(Matrix[2][5][0]), float(Matrix[2][5][1]))
                d_red = random.uniform(float(Matrix[6][5][0]), float(Matrix[6][5][1]))
                # costheta0
                costheta0_yellow = random.uniform(float(Matrix[0][6][0]), float(Matrix[0][6][1]))
                costheta0_green = random.uniform(float(Matrix[2][6][0]), float(Matrix[2][6][1]))
                costheta0_red = random.uniform(float(Matrix[6][6][0]), float(Matrix[6][6][1]))
                # 2 body parameters
                # n
                n = random.uniform(float(Matrix[3][7][0]), float(Matrix[3][7][1]))
                # lambda2
                lambda2 = random.uniform(float(Matrix[3][9][0]), float(Matrix[3][9][1]))
                # B
                B = random.uniform(float(Matrix[3][10][0]), float(Matrix[3][10][1]))
                # lambda1
                lambda1 = random.uniform(float(Matrix[3][13][0]), float(Matrix[3][13][1]))
                # A
                A = random.uniform(float(Matrix[3][14][0]), float(Matrix[3][14][1]))

                for row in range(len(Matrix)):
                    file.write(f"{Matrix[row][0][0]}\t{Matrix[row][0][1]}\t{Matrix[row][0][2]}\t")
                    if row <= 1:
                        file.write(f"{random.uniform(float(Matrix[row][1][0]), float(Matrix[row][1][1])):.9f}\t")
                        file.write(f"{(gamma_yellow):.9f}\t")
                        file.write(f"{(lambda3_yellow):.9f}\t")
                        file.write(f"{(c_yellow):.9f}\t")
                        file.write(f"{(d_yellow):.9f}\t")
                        file.write(f"{(costheta0_yellow):.9f}\t")
                        for column in range(7, 19):
                            file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                    elif row in [2,3,4,5]:
                        file.write(f"{random.uniform(float(Matrix[row][1][0]), float(Matrix[row][1][1])):.9f}\t")
                        file.write(f"{(gamma_green):.9f}\t")
                        file.write(f"{(lambda3_green):.9f}\t")
                        file.write(f"{(c_green):.9f}\t")
                        file.write(f"{(d_green):.9f}\t")
                        file.write(f"{(costheta0_green):.9f}\t")
                        for column in range(7, 19):
                            if (row == 3 or row == 4) and column == 7:
                                file.write(f"{n}\t")
                            elif (row == 3 or row == 4) and column == 9:
                                file.write(f"{lambda2}\t")
                            elif (row == 3 or row == 4) and column == 10:
                                file.write(f"{B}\t")
                            elif (row == 3 or row == 4) and column == 13:
                                file.write(f"{lambda1}\t")
                            elif (row == 3 or row == 4) and column == 14:
                                file.write(f"{A}\t")
                            else:
                                file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                    else:
                        file.write(f"{random.uniform(float(Matrix[row][1][0]), float(Matrix[row][1][1])):.9f}\t")
                        file.write(f"{(gamma_red):.9f}\t")
                        file.write(f"{(lambda3_red):.9f}\t")
                        file.write(f"{(c_red):.9f}\t")
                        file.write(f"{(d_red):.9f}\t")
                        file.write(f"{(costheta0_red):.9f}\t")
                        for column in range(7, 19):
                            file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                    file.write("\n\n")

elif Constraints == "ON" and PotentialType == "MOD":
    def FirstGenCreator(Matrix ,PopulationSize):
        Folder_check("Generation_1")
        for p in range(1, PopulationSize+1):
            file_name = os.path.join("Generation_1", f"001-{p:03d}.tersoff.modc")
            with open(file_name, 'a') as file:

                # 3 body parameters
                # alpha
                alpha_yellow = random.uniform(float(Matrix[0][2][0]), float(Matrix[0][2][1]))
                alpha_green = random.uniform(float(Matrix[2][2][0]), float(Matrix[2][2][1]))
                alpha_red = random.uniform(float(Matrix[6][2][0]), float(Matrix[6][2][1]))
                # h
                h_yellow = random.uniform(float(Matrix[0][3][0]), float(Matrix[0][3][1]))
                h_green = random.uniform(float(Matrix[2][3][0]), float(Matrix[2][3][1]))
                h_red = random.uniform(float(Matrix[6][3][0]), float(Matrix[6][3][1]))
                # c1
                c1_yellow = random.uniform(float(Matrix[0][13][0]), float(Matrix[0][13][1]))
                c1_green = random.uniform(float(Matrix[2][13][0]), float(Matrix[2][13][1]))
                c1_red = random.uniform(float(Matrix[6][13][0]), float(Matrix[6][13][1]))
                # c2
                c2_yellow = random.uniform(float(Matrix[0][14][0]), float(Matrix[0][14][1]))
                c2_green = random.uniform(float(Matrix[2][14][0]), float(Matrix[2][14][1]))
                c2_red = random.uniform(float(Matrix[6][14][0]), float(Matrix[6][14][1]))
                # c3
                c3_yellow = random.uniform(float(Matrix[0][15][0]), float(Matrix[0][15][1]))
                c3_green = random.uniform(float(Matrix[2][15][0]), float(Matrix[2][15][1]))
                c3_red = random.uniform(float(Matrix[6][15][0]), float(Matrix[6][15][1]))
                # c4
                c4_yellow = random.uniform(float(Matrix[0][16][0]), float(Matrix[0][16][1]))
                c4_green = random.uniform(float(Matrix[2][16][0]), float(Matrix[2][16][1]))
                c4_red = random.uniform(float(Matrix[6][16][0]), float(Matrix[6][16][1]))
                # c5
                c5_yellow = random.uniform(float(Matrix[0][17][0]), float(Matrix[0][17][1]))
                c5_green = random.uniform(float(Matrix[2][17][0]), float(Matrix[2][17][1]))
                c5_red = random.uniform(float(Matrix[6][17][0]), float(Matrix[6][17][1]))
                # c0
                c0_yellow = random.uniform(float(Matrix[0][18][0]), float(Matrix[0][18][1]))
                c0_green = random.uniform(float(Matrix[2][18][0]), float(Matrix[2][18][1]))
                c0_red = random.uniform(float(Matrix[6][18][0]), float(Matrix[6][18][1]))
                # 2 body parameters
                # eta
                eta = random.uniform(float(Matrix[3][4][0]), float(Matrix[3][4][1]))
                lambda2 = random.uniform(float(Matrix[3][6][0]), float(Matrix[3][6][1]))
                B = random.uniform(float(Matrix[3][7][0]), float(Matrix[3][7][1]))
                lambda1 = random.uniform(float(Matrix[3][10][0]), float(Matrix[3][10][1]))
                A = random.uniform(float(Matrix[3][11][0]), float(Matrix[3][11][1]))
                n = random.uniform(float(Matrix[3][12][0]), float(Matrix[3][12][1]))


                for row in range(len(Matrix)):
                    file.write(f"{Matrix[row][0][0]}\t{Matrix[row][0][1]}\t{Matrix[row][0][2]}\t")
                    if row <= 1:
                        file.write(f"{random.uniform(float(Matrix[row][1][0]), float(Matrix[row][1][1])):.9f}\t")
                        file.write(f"{(alpha_yellow):.9f}\t")
                        file.write(f"{(h_yellow):.9f}\t")
                        for column in range(4, 13):
                            file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                        file.write(f"{(c1_yellow):.9f}\t")
                        file.write(f"{(c2_yellow):.9f}\t")
                        file.write(f"{(c3_yellow):.9f}\t")
                        file.write(f"{(c4_yellow):.9f}\t")
                        file.write(f"{(c5_yellow):.9f}\t")
                        file.write(f"{(c0_yellow):.9f}\t")
                    elif row in [2,3,4,5]:
                        file.write(f"{random.uniform(float(Matrix[row][1][0]), float(Matrix[row][1][1])):.9f}\t")
                        file.write(f"{(alpha_green):.9f}\t")
                        file.write(f"{(h_green):.9f}\t")
                        for column in range(4, 13):
                            if (row == 3 or row == 4) and column == 4:
                                file.write(f"{eta}\t")
                            elif (row == 3 or row == 4) and column == 6:
                                file.write(f"{lambda2}\t")
                            elif (row == 3 or row == 4) and column == 7:
                                file.write(f"{B}\t")
                            elif (row == 3 or row == 4) and column == 10:
                                file.write(f"{lambda1}\t")
                            elif (row == 3 or row == 4) and column == 11:
                                file.write(f"{A}\t")
                            elif (row == 3 or row == 4) and column == 12:
                                file.write(f"{n}\t")
                            else:
                                file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                        file.write(f"{(c1_green):.9f}\t")
                        file.write(f"{(c2_green):.9f}\t")
                        file.write(f"{(c3_green):.9f}\t")
                        file.write(f"{(c4_green):.9f}\t")
                        file.write(f"{(c5_green):.9f}\t")
                        file.write(f"{(c0_green):.9f}\t")
                    else:
                        file.write(f"{random.uniform(float(Matrix[row][1][0]), float(Matrix[row][1][1])):.9f}\t")
                        file.write(f"{(alpha_red):.9f}\t")
                        file.write(f"{(h_red):.9f}\t")
                        for column in range(4, 13):
                            file.write(f"{random.uniform(float(Matrix[row][column][0]), float(Matrix[row][column][1])):.9f}\t")
                        file.write(f"{(c1_red):.9f}\t")
                        file.write(f"{(c2_red):.9f}\t")
                        file.write(f"{(c3_red):.9f}\t")
                        file.write(f"{(c4_red):.9f}\t")
                        file.write(f"{(c5_red):.9f}\t")
                        file.write(f"{(c0_red):.9f}\t")
                    file.write("\n\n")

    # print("First Generation is ready!")
                    
FirstGenCreator(Matrix, PopulationSize)            

# import sys
# sys.exit()
# 5th step is creating the atom_definition files from the data base
#===================================#
Folder_check("Definitions")

def read_database(file_name):
    folder_name = "Targets"
    blank_count = 0
    targets_file_path = os.path.join(folder_name, file_name)
    with open(targets_file_path, 'r') as file:
        lines = file.readlines()
        data = []
        current_list = []
        for line in lines:
            line = line.strip()
            if not line:  # Check for a blank line
                blank_count += 1
                data.append(current_list)
                current_list = []
                continue
            else:
                line = line.split()
                current_list.append(line)
                # data.append(current_list)
    return data, blank_count

data, count = read_database(NameofModifiedDataBase)

def Get_elements(matrix):
    elements = []
    for row in matrix:
        first_element = row[0][0]
        if first_element not in elements:
            elements.append(first_element)
    return elements

elements = Get_elements(Matrix)
# print(elements)

def write_definitions(data_list, number_of_bases, element_list):
    masses = []
    list_of_defs = []
    # Open and read the "Library.txt" file
    with open('Specific/Library.txt', 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split()
            for element in element_list:
                if line[0] == element:
                    masses.append(line)
    extension = '.txt'
    folder = 'Definitions'
    atoms = len(element_list)
    for i in range(number_of_bases):
        new_filename = f"atom_definition_{i+1}" + extension
        list_of_defs.append(new_filename)
        file_name = os.path.join(folder, new_filename)
        with open(file_name, 'w') as file:
            file.write(" Start File for LAMMPS\n")
            # file.write(f"#free_energy = {float(data_list[i][2][0])}\n")
            file.write(f" {int(data_list[i][0][0])} atoms\n\n")
            file.write(f" {atoms} atom types\n\n")
            file.write(f"   0   {float(data_list[i][1][0])} xlo xhi\n")
            file.write(f"   0   {float(data_list[i][1][4])} ylo yhi\n")
            file.write(f"   0   {float(data_list[i][1][8])} zlo zhi\n\n")
            file.write(" Masses\n\n")
            for j in range(atoms):
                file.write(f" {j+1}  {float(masses[j][1])} # {masses[j][0]}\n")
            file.write("\n\n")
            file.write("Atoms\n\n")
            num_atoms = int(data_list[i][0][0])
            for k in range(3, num_atoms+3):
                for l in range(atoms):
                    # if data_list[i][k][0] == element_list[l]:
                    if data_list[i][k][0] == masses[l][0]:
                        file.write(f"{k-2}\t{l+1}\t{float(data_list[i][k][1]):.9f}\t{float(data_list[i][k][2]):.9f}\t{float(data_list[i][k][3]):.9f}\n")
        # print(f"atom_definition_{i+1} is done!")
    # print("Atom Definition files are done!")
    return list_of_defs, masses

atom_defs, masses = write_definitions(data, count, elements)


# print(atom_defs)
# print(masses)
# 6th step is geting the expected forces and energies from database:
#===================================#
def getForcesAndEnergy_fromData(Data, counter):
    Energies = []
    FxFyFz = []
    for i in range(counter):
        FofData = []
        Energies.append(float(Data[i][2][0]))
        for j in range(int(Data[i][0][0])):
            ForceLine = []
            for k in [4, 5, 6]:
                ForceLine.append(float(Data[i][j+3][k]))
            FofData.append(ForceLine)
        FxFyFz.append(FofData)
    return Energies, FxFyFz

ExpectedEnergies_from_DataBase, ExpectedForces_from_DataBase = getForcesAndEnergy_fromData(data, count)

# print(ExpectedEnergies_from_DataBase)
# print(ExpectedForces_from_DataBase)

# import sys
# sys.exit()

def write_inputs(gen_number, List_of_definitions, masses, PotentialType):
    Potentials = os.listdir(f"Generation_{gen_number}")
    # Generate all possible combinations
    combinations = list(itertools.product(Potentials, List_of_definitions))
    for i in range(len(combinations)):
        new_filename = f"input_{i+1}.txt"
        file_name = os.path.join("Inputs", new_filename)
        with open(file_name, 'w') as file:
            file.write("# Initialization\n\n")
            file.write("units		  metal\ndimension 	      3\nboundary	  p p p\natom_style	 atomic\n\n\n")
            file.write("# Atom Definition\n\n")
            file.write(f"read_data {combinations[i][1]}\n\n\n")
            file.write("# Settings\n\n")
            if PotentialType == "MOD":
                file.write("pair_style tersoff/mod/c\n")
            elif PotentialType == "TER":
                file.write("pair_style tersoff\n")
            file.write(f"pair_coeff	* * {combinations[i][0]} ")
            for j in range(len(masses)):
                file.write(f"{masses[j][0]} ")
            file.write("\n\n")
            # file.write("compute pe all pe/atom pair\n\n")
            file.write("dump myDump all custom 10 positions_and_forces.txt id type x y z fx fy fz\n")
            file.write("dump_modify myDump sort id\n\n")
            file.write("thermo 10\n\n\n# Run a simulation\n\nrun              1\n\n")
            file.write('variable pe equal "pe"\n')
            file.write('print "Potential energy (eV) = ${pe};"\n')
            file.write('print "All done!"')

# 7th step is going to be the main loop:
#===================================#

MeanFitnessOfGen = [0 for i in range(Generations)]
# s_t_d = [0 for i in range(Generations)]
BestFitnessList = [0 for i in range(Generations)]

ALL_ener = []
ALL_forc = []

# Folder_check("FoldWithResults")
xCoord = []
WE = 1
WF = 1

# this function copies a file from folder 1 to folder 2
def copy_file(source_folder, source_file, destination_folder):
    shutil.copy(os.path.join(source_folder, source_file), os.path.join(destination_folder, source_file))
    

# Collect the energy from log.lammps file for a potential
def getEnergy_fromLog(folder):
    Energy = None
    with open(os.path.join(os.getcwd(), folder, "log.lammps")) as logg:
        lines = logg.readlines()
        for targetLine in lines:
            if re.search(r'Potential\s*energy\s*\(eV\)\s*=\s*([-+]?\d*\.\d+|\d+)', targetLine):
                Energy = float(re.search(r'Potential\s*energy\s*\(eV\)\s*=\s*([-+]?\d*\.\d+|\d+)', targetLine).group(1))
    return Energy

# save the forces of the positions_and_forces.txt in a list
def getForces_fromOutput(folder):
    with open(os.path.join(os.getcwd(), folder, "positions_and_forces.txt"), 'r') as file:
        lines = file.readlines()[9:]
        forcesList = []
        for line in lines:
            forcesFloat = []
            # Split the line into elements
            elements = line.split()
            # Extract the last 3 elements
            last_three_elements = elements[-3:]
            # Convert the elements to floats and append to the list
            forcesFloat = [float(element) for element in last_three_elements]
            forcesList.append(forcesFloat)
    return forcesList

# This function computes the fitness
def computeFitness(DataBaseEnergies, DataBaseForces, EnergyFromLog, ForciesFromOutput, D):
    d = int(D)
    DataBaseEnergies = DataBaseEnergies[d-1]
    DataBaseForces = DataBaseForces[d-1]
    energy_error = (EnergyFromLog - float(DataBaseEnergies))**2
    force_error = []
    for i in range(len(ForciesFromOutput)):
        force_error.append((float(DataBaseForces[i][0])-ForciesFromOutput[i][0])**2 
                                    + (float(DataBaseForces[i][1])-ForciesFromOutput[i][1])**2
                                    + (float(DataBaseForces[i][2])-ForciesFromOutput[i][2])**2)
    TotalforceError = sum(force_error)
    return energy_error, TotalforceError

# Create the results.txt by reading the Origins_Gen_#.txt from Origins Folder
def make_result(We, Wf, gen, FinalEnergies, FinalForces):
    name_list = []
    origins = []
    bestbest = float('inf')
    with open(os.path.join(os.getcwd(), "Origins", f"Origins_Gen_{gen}.txt"), 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            columns = line.split()
            name_list.append(columns[0])
            origins.append(columns[1])
    with open("results.txt", 'w') as ff:
        ff.write(" ID\tOrigin\tForce_fitness\tEnergy_fitness\n")
        for i in range(len(name_list)):
            if float(We * FinalEnergies[i] + Wf * FinalForces[i]) < bestbest:
                bestbest = float(We * FinalEnergies[i] + Wf * FinalForces[i])
            ff.write(f"{name_list[i]}\t{origins[i]}\t{FinalForces[i]:.9f}\t{FinalEnergies[i]:.9f}\n")
    return bestbest

# Computing the final fitness (forces and energy together)
def fit_list(wE, wF, final_PopulationSize):
    Fitness_Results = []
    with open("results.txt") as results_file:
        lines = results_file.readlines()
        for line in lines[1:final_PopulationSize+1]:
            data = re.findall(r'\S+', line)
            if len(data) >= 4:
                Fitness_Results.append((data[0] ,wF*float(data[2]) + wE*float(data[3])))
    sorted_list = sorted(Fitness_Results, key=lambda x: x[1])
    
    return sorted_list


def read_potential_values(filename, genA):
    # genA = get_current_generation()
    file_path = os.path.join(f"Generation_{genA}", filename)
    # print(filename)
    with open(file_path, 'r') as f:
        lines = f.read().split('\n\n')
        all_lines = []
        for line in lines:
            if not line:
                continue
            elements, values = line.strip().split('\t')[:3], line.strip().split('\t')[3:]
            # print(values)
            float_values = []
            for i in range(len(values)):
                if i == 0 or i == 4 or i == 7 or i == 8:
                    float_values.append(values[i])
                else:
                    float_values.append(round(float(values[i]), 9))
            line_list = elements + float_values
            all_lines.append(line_list)
            if all_lines[-1] == ['']:
                all_lines.pop()

    return all_lines

def save_potential(folder, filename, values):
    file_path = os.path.join(folder, filename)
    with open(file_path, 'w') as f:
        for i in range(len(values)):
            for j in range(len(values[0])):
                f.write(f"{values[i][j]}\t")
            f.write("\n\n")




def select_random_parents(sorted_list, BestFracCount, gen):
    remaining_potentials = sorted_list[:BestFracCount]

    canditates = random.sample(remaining_potentials, 8)
    canditates = sorted(canditates, key=lambda x: x[1])
    #canditates = canditates[:4]
    canditates = random.sample(canditates, 4)
    canditates = random.sample(canditates, 2)


    parent1, parent2 = canditates

    parent_A = parent1[0] + '.tersoff.modc'
    parent_B = parent2[0] + '.tersoff.modc'
    values_parentA = read_potential_values(parent_A, gen)
    values_parentB = read_potential_values(parent_B, gen)

    return values_parentA, values_parentB, parent1[0], parent2[0]





def crossover(parent1, parent2):
    child = []
    splitter = random.randint(5,18)
    # splitter = 9
    for i in range(len(parent1)):
        small_list = []
        for j in range(len(parent1[0])):
            if j == 0 or j == 1 or j == 2:
                small_list.append(parent1[i][j])
            elif j > 2 and j <= splitter:
                small_list.append(parent1[i][j])
            else:
                small_list.append(parent2[i][j])
        child.append(small_list)
    return child


if Constraints == "OFF":
    def mutate(values, MUt_PERC, Matrix):
        for i in range(len(values)):
            for j in range(len(values[0])):
                if j >= 3:
                    if Matrix[i][j-2][0] != Matrix[i][j-2][1]:
                        
                        values[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                        values[i][j] = round(values[i][j], 9)
                    else:
                        continue
        mutated_values = values
        return mutated_values
elif Constraints == "ON" and PotentialType == "TER":
    def mutate(values, MUt_PERC, Matrix):
        gamma_yellow = random.uniform(float(values[0][4]) * (1 - MUt_PERC), float(values[0][4]) * (1 + MUt_PERC))
        gamma_green = random.uniform(float(values[2][4]) * (1 - MUt_PERC), float(values[2][4]) * (1 + MUt_PERC))
        gamma_red = random.uniform(float(values[6][4]) * (1 - MUt_PERC), float(values[6][4]) * (1 + MUt_PERC))

        lambda3_yellow = random.uniform(float(values[0][5]) * (1 - MUt_PERC), float(values[0][5]) * (1 + MUt_PERC))
        lambda3_green = random.uniform(float(values[2][5]) * (1 - MUt_PERC), float(values[2][5]) * (1 + MUt_PERC))
        lambda3_red = random.uniform(float(values[6][5]) * (1 - MUt_PERC), float(values[6][5]) * (1 + MUt_PERC))

        c_yellow = random.uniform(float(values[0][6]) * (1 - MUt_PERC), float(values[0][6]) * (1 + MUt_PERC))
        c_green = random.uniform(float(values[2][6]) * (1 - MUt_PERC), float(values[2][6]) * (1 + MUt_PERC))
        c_red = random.uniform(float(values[6][6]) * (1 - MUt_PERC), float(values[6][6]) * (1 + MUt_PERC))

        d_yellow = random.uniform(float(values[0][7]) * (1 - MUt_PERC), float(values[0][7]) * (1 + MUt_PERC))
        d_green = random.uniform(float(values[2][7]) * (1 - MUt_PERC), float(values[2][7]) * (1 + MUt_PERC))
        d_red = random.uniform(float(values[6][7]) * (1 - MUt_PERC), float(values[6][7]) * (1 + MUt_PERC))

        costheta0_yellow = random.uniform(float(values[0][8]) * (1 - MUt_PERC), float(values[0][8]) * (1 + MUt_PERC))
        costheta0_green = random.uniform(float(values[2][8]) * (1 - MUt_PERC), float(values[2][8]) * (1 + MUt_PERC))
        costheta0_red = random.uniform(float(values[6][8]) * (1 - MUt_PERC), float(values[6][8]) * (1 + MUt_PERC))

        n = random.uniform(float(values[3][9]) * (1 - MUt_PERC), float(values[3][9]) * (1 + MUt_PERC))
        lambda2 = random.uniform(float(values[3][11]) * (1 - MUt_PERC), float(values[3][11]) * (1 + MUt_PERC))
        B = random.uniform(float(values[3][12]) * (1 - MUt_PERC), float(values[3][12]) * (1 + MUt_PERC))
        lambda1 = random.uniform(float(values[3][15]) * (1 - MUt_PERC), float(values[3][15]) * (1 + MUt_PERC))
        A = random.uniform(float(values[3][16]) * (1 - MUt_PERC), float(values[3][16]) * (1 + MUt_PERC))

        Rows = len(values)
        Columns = len(values[0])
        Helper_array = [["nothing"] * Columns for _ in range(Rows)]

        # Yellows
        Helper_array[0][4] = Helper_array[1][4] = round(gamma_yellow, 9)
        Helper_array[0][5] = Helper_array[1][5] = round(lambda3_yellow, 9)
        Helper_array[0][6] = Helper_array[1][6] = round(c_yellow, 9)
        Helper_array[0][7] = Helper_array[1][7] = round(d_yellow, 9)
        Helper_array[0][8] = Helper_array[1][8] = round(costheta0_yellow, 9)

        # Greens
        Helper_array[2][4] = Helper_array[3][4] = Helper_array[4][4] = Helper_array[5][4] = round(gamma_green, 9)
        Helper_array[2][5] = Helper_array[3][5] = Helper_array[4][5] = Helper_array[5][5] = round(lambda3_green, 9)
        Helper_array[2][6] = Helper_array[3][6] = Helper_array[4][6] = Helper_array[5][6] = round(c_green, 9)
        Helper_array[2][7] = Helper_array[3][7] = Helper_array[4][7] = Helper_array[5][7] = round(d_green, 9)
        Helper_array[2][8] = Helper_array[3][8] = Helper_array[4][8] = Helper_array[5][8] = round(costheta0_green, 9)
        
        # Reds
        Helper_array[6][4] = Helper_array[7][4] = round(gamma_red, 9)
        Helper_array[6][5] = Helper_array[7][5] = round(lambda3_red, 9)
        Helper_array[6][6] = Helper_array[7][6] = round(c_red, 9)
        Helper_array[6][7] = Helper_array[7][7] = round(d_red, 9)
        Helper_array[6][8] = Helper_array[7][8] = round(costheta0_red, 9)

        # Blue
        Helper_array[3][9] = Helper_array[4][9] = round(n, 9)
        Helper_array[3][11] = Helper_array[4][11] = round(lambda2, 9)
        Helper_array[3][12] = Helper_array[4][12] = round(B, 9)
        Helper_array[3][15] = Helper_array[4][15] = round(lambda1, 9)
        Helper_array[3][16] = Helper_array[4][16] = round(A, 9)

        for i in range(Rows):
            for j in range(Columns):
                if j >= 3:
                    if Matrix[i][j-2][0] != Matrix[i][j-2][1]:
                        if Helper_array[i][j] != "nothing":
                            values[i][j] = Helper_array[i][j]
                            # Helper_array[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                            # Helper_array[i][j] = round(Helper_array[i][j], 9)
                        else:
                            values[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                            values[i][j] = round(values[i][j], 9)
                        # if i in [0,1] and j == 4:
                        #     values[i][j] = round(gamma_yellow, 9)
                        # elif i in [0,1] and j == 5:
                        #     values[i][j] = round(lambda3_yellow, 9)
                        # elif i in [0,1] and j == 6:
                        #     values[i][j] = round(c_yellow, 9)
                        # elif i in [0,1] and j == 7:
                        #     values[i][j] = round(d_yellow, 9)
                        # elif i in [0,1] and j == 8:
                        #     values[i][j] = round(costheta0_yellow, 9)

                        # if i in [2,3,4,5] and j == 4:
                        #     values[i][j] = round(gamma_green, 9)
                        # elif i in [2,3,4,5] and j == 5:
                        #     values[i][j] = round(lambda3_green, 9)
                        # elif i in [2,3,4,5] and j == 6:
                        #     values[i][j] = round(c_green, 9)
                        # elif i in [2,3,4,5] and j == 7:
                        #     values[i][j] = round(d_green, 9)
                        # elif i in [2,3,4,5] and j == 8:
                        #     values[i][j] = round(costheta0_green, 9)

                        # if i in [6,7] and j == 4:
                        #     values[i][j] = round(gamma_red, 9)
                        # elif i in [6,7] and j == 5:
                        #     values[i][j] = round(lambda3_red, 9)
                        # elif i in [6,7] and j == 6:
                        #     values[i][j] = round(c_red, 9)
                        # elif i in [6,7] and j == 7:
                        #     values[i][j] = round(d_red, 9)
                        # elif i in [6,7] and j == 8:
                        #     values[i][j] = round(costheta0_red, 9)

                        # if i in [3,4] and j == 9:
                        #     values[i][j] = round(n, 9)
                        # elif i in [3,4] and j == 11:
                        #     values[i][j] = round(lambda2, 9)
                        # elif i in [3,4] and j == 12:
                        #     values[i][j] = round(B, 9)
                        # elif i in [3,4] and j == 15:
                        #     values[i][j] = round(lambda1, 9)
                        # elif i in [3,4] and j == 16:
                        #     values[i][j] = round(A, 9)
                        # else:
                        #     values[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                        #     values[i][j] = round(values[i][j], 9)
                    else:
                        continue

        
        mutated_values = values
        return mutated_values
elif Constraints == "ON" and PotentialType == "MOD":
    def mutate(values, MUt_PERC, Matrix):
        
        alpha_yellow = random.uniform(float(values[0][4]) * (1 - MUt_PERC), float(values[0][4]) * (1 + MUt_PERC))
        alpha_green = random.uniform(float(values[2][4]) * (1 - MUt_PERC), float(values[2][4]) * (1 + MUt_PERC))
        alpha_red = random.uniform(float(values[6][4]) * (1 - MUt_PERC), float(values[6][4]) * (1 + MUt_PERC))

        h_yellow = random.uniform(float(values[0][5]) * (1 - MUt_PERC), float(values[0][5]) * (1 + MUt_PERC))
        h_green = random.uniform(float(values[2][5]) * (1 - MUt_PERC), float(values[2][5]) * (1 + MUt_PERC))
        h_red = random.uniform(float(values[6][5]) * (1 - MUt_PERC), float(values[6][5]) * (1 + MUt_PERC))

        c1_yellow = random.uniform(float(values[0][15]) * (1 - MUt_PERC), float(values[0][15]) * (1 + MUt_PERC))
        c2_yellow = random.uniform(float(values[0][16]) * (1 - MUt_PERC), float(values[0][16]) * (1 + MUt_PERC))
        c3_yellow = random.uniform(float(values[0][17]) * (1 - MUt_PERC), float(values[0][17]) * (1 + MUt_PERC))
        c4_yellow = random.uniform(float(values[0][18]) * (1 - MUt_PERC), float(values[0][18]) * (1 + MUt_PERC))
        c5_yellow = random.uniform(float(values[0][19]) * (1 - MUt_PERC), float(values[0][19]) * (1 + MUt_PERC))
        c0_yellow = random.uniform(float(values[0][20]) * (1 - MUt_PERC), float(values[0][20]) * (1 + MUt_PERC))

        c1_green = random.uniform(float(values[2][15]) * (1 - MUt_PERC), float(values[2][15]) * (1 + MUt_PERC))
        c2_green = random.uniform(float(values[2][16]) * (1 - MUt_PERC), float(values[2][16]) * (1 + MUt_PERC))
        c3_green = random.uniform(float(values[2][17]) * (1 - MUt_PERC), float(values[2][17]) * (1 + MUt_PERC))
        c4_green = random.uniform(float(values[2][18]) * (1 - MUt_PERC), float(values[2][18]) * (1 + MUt_PERC))
        c5_green = random.uniform(float(values[2][19]) * (1 - MUt_PERC), float(values[2][19]) * (1 + MUt_PERC))
        c0_green = random.uniform(float(values[2][20]) * (1 - MUt_PERC), float(values[2][20]) * (1 + MUt_PERC))

        c1_red = random.uniform(float(values[6][15]) * (1 - MUt_PERC), float(values[6][15]) * (1 + MUt_PERC))
        c2_red = random.uniform(float(values[6][16]) * (1 - MUt_PERC), float(values[6][16]) * (1 + MUt_PERC))
        c3_red = random.uniform(float(values[6][17]) * (1 - MUt_PERC), float(values[6][17]) * (1 + MUt_PERC))
        c4_red = random.uniform(float(values[6][18]) * (1 - MUt_PERC), float(values[6][18]) * (1 + MUt_PERC))
        c5_red = random.uniform(float(values[6][19]) * (1 - MUt_PERC), float(values[6][19]) * (1 + MUt_PERC))
        c0_red = random.uniform(float(values[6][20]) * (1 - MUt_PERC), float(values[6][20]) * (1 + MUt_PERC))
        

        eta = random.uniform(float(values[3][6]) * (1 - MUt_PERC), float(values[3][6]) * (1 + MUt_PERC))
        lambda2 = random.uniform(float(values[3][8]) * (1 - MUt_PERC), float(values[3][8]) * (1 + MUt_PERC))
        B = random.uniform(float(values[3][9]) * (1 - MUt_PERC), float(values[3][9]) * (1 + MUt_PERC))
        lambda1 = random.uniform(float(values[3][12]) * (1 - MUt_PERC), float(values[3][12]) * (1 + MUt_PERC))
        A = random.uniform(float(values[3][13]) * (1 - MUt_PERC), float(values[3][13]) * (1 + MUt_PERC))
        n = random.uniform(float(values[3][14]) * (1 - MUt_PERC), float(values[3][14]) * (1 + MUt_PERC))


        Rows = len(values)
        Columns = len(values[0])
        Helper_array = [["nothing"] * Columns for _ in range(Rows)]

        # Yellows
        Helper_array[0][4] = Helper_array[1][4] = round(alpha_yellow, 9)
        Helper_array[0][5] = Helper_array[1][5] = round(h_yellow, 9)
        Helper_array[0][15] = Helper_array[1][15] = round(c1_yellow, 9)
        Helper_array[0][16] = Helper_array[1][16] = round(c2_yellow, 9)
        Helper_array[0][17] = Helper_array[1][17] = round(c3_yellow, 9)
        Helper_array[0][18] = Helper_array[1][18] = round(c4_yellow, 9)
        Helper_array[0][19] = Helper_array[1][19] = round(c5_yellow, 9)
        Helper_array[0][20] = Helper_array[1][20] = round(c0_yellow, 9)
        

        # Greens
        Helper_array[2][4] = Helper_array[3][4] = Helper_array[4][4] = Helper_array[5][4] = round(alpha_green, 9)
        Helper_array[2][5] = Helper_array[3][5] = Helper_array[4][5] = Helper_array[5][5] = round(h_green, 9)
        Helper_array[2][15] = Helper_array[3][15] = Helper_array[4][15] = Helper_array[5][15] = round(c1_green, 9)
        Helper_array[2][16] = Helper_array[3][16] = Helper_array[4][16] = Helper_array[5][16] = round(c2_green, 9)
        Helper_array[2][17] = Helper_array[3][17] = Helper_array[4][17] = Helper_array[5][17] = round(c3_green, 9)
        Helper_array[2][18] = Helper_array[3][18] = Helper_array[4][18] = Helper_array[5][18] = round(c4_green, 9)
        Helper_array[2][19] = Helper_array[3][19] = Helper_array[4][19] = Helper_array[5][19] = round(c5_green, 9)
        Helper_array[2][20] = Helper_array[3][20] = Helper_array[4][20] = Helper_array[5][20] = round(c0_green, 9)
        
        
        # Reds
        Helper_array[6][4] = Helper_array[7][4] = round(alpha_red, 9)
        Helper_array[6][5] = Helper_array[7][5] = round(h_red, 9)
        Helper_array[6][15] = Helper_array[7][15] = round(c1_red, 9)
        Helper_array[6][16] = Helper_array[7][16] = round(c2_red, 9)
        Helper_array[6][17] = Helper_array[7][17] = round(c3_red, 9)
        Helper_array[6][18] = Helper_array[7][18] = round(c4_red, 9)
        Helper_array[6][19] = Helper_array[7][19] = round(c5_red, 9)
        Helper_array[6][20] = Helper_array[7][20] = round(c0_red, 9)

        # Blue
        Helper_array[3][6] = Helper_array[4][6] = round(eta, 9)
        Helper_array[3][8] = Helper_array[4][8] = round(lambda2, 9)
        Helper_array[3][9] = Helper_array[4][9] = round(B, 9)
        Helper_array[3][12] = Helper_array[4][12] = round(lambda1, 9)
        Helper_array[3][13] = Helper_array[4][13] = round(A, 9)
        Helper_array[3][14] = Helper_array[4][14] = round(n, 9)



        for i in range(Rows):
            for j in range(Columns):
                if j >= 3:
                    if Matrix[i][j-2][0] != Matrix[i][j-2][1]:
                        if Helper_array[i][j] != "nothing":
                            values[i][j] = Helper_array[i][j]
                            # Helper_array[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                            # Helper_array[i][j] = round(Helper_array[i][j], 9)
                        else:
                            values[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                            values[i][j] = round(values[i][j], 9)
                    # if Matrix[i][j-2][0] != Matrix[i][j-2][1]:
                    #     if (i == 0 or i == 1) and j == 4:
                    #         values[i][j] = round(alpha_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 5:
                    #         values[i][j] = round(h_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 15:
                    #         values[i][j] = round(c1_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 16:
                    #         values[i][j] = round(c2_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 17:
                    #         values[i][j] = round(c3_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 18:
                    #         values[i][j] = round(c4_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 19:
                    #         values[i][j] = round(c5_yellow, 9)
                    #     elif (i == 0 or i == 1) and j == 20:
                    #         values[i][j] = round(c0_yellow, 9)

                    #     elif i in [2,3,4,5] and j == 4:
                    #         values[i][j] = round(alpha_green, 9)
                    #     elif i in [2,3,4,5] and j == 5:
                    #         values[i][j] = round(h_green, 9)
                    #     elif i in [2,3,4,5] and j == 15:
                    #         values[i][j] = round(c1_green, 9)
                    #     elif i in [2,3,4,5] and j == 16:
                    #         values[i][j] = round(c2_green, 9)
                    #     elif i in [2,3,4,5] and j == 17:
                    #         values[i][j] = round(c3_green, 9)
                    #     elif i in [2,3,4,5] and j == 18:
                    #         values[i][j] = round(c4_green, 9)
                    #     elif i in [2,3,4,5] and j == 19:
                    #         values[i][j] = round(c5_green, 9)
                    #     elif i in [2,3,4,5] and j == 20:
                    #         values[i][j] = round(c0_green, 9)

                    #     elif i in [3,4] and j == 6:
                    #         values[i][j] = round(eta, 9)
                    #     elif i in [3,4] and j == 8:
                    #         values[i][j] = round(lambda2, 9)
                    #     elif i in [3,4] and j == 9:
                    #         values[i][j] = round(B, 9)
                    #     elif i in [3,4] and j == 12:
                    #         values[i][j] = round(lambda1, 9)
                    #     elif i in [3,4] and j == 13:
                    #         values[i][j] = round(A, 9)
                    #     elif i in [3,4] and j == 14:
                    #         values[i][j] = round(n, 9)


                    #     elif (i == 6 or i == 7) and j == 4:
                    #         values[i][j] = round(alpha_red, 9)
                    #     elif (i == 6 or i == 7) and j == 5:
                    #         values[i][j] = round(h_red, 9)
                    #     elif (i == 6 or i == 7) and j == 15:
                    #         values[i][j] = round(c1_red, 9)
                    #     elif (i == 6 or i == 7) and j == 16:
                    #         values[i][j] = round(c2_red, 9)
                    #     elif (i == 6 or i == 7) and j == 17:
                    #         values[i][j] = round(c3_red, 9)
                    #     elif (i == 6 or i == 7) and j == 18:
                    #         values[i][j] = round(c4_red, 9)
                    #     elif (i == 6 or i == 7) and j == 19:
                    #         values[i][j] = round(c5_red, 9)
                    #     elif (i == 6 or i == 7) and j == 20:
                    #         values[i][j] = round(c0_red, 9)

                    #     else:
                    #         values[i][j] = random.uniform(float(values[i][j]) * (1 - MUt_PERC), float(values[i][j]) * (1 + MUt_PERC))
                    #         values[i][j] = round(values[i][j], 9)
                    else:
                        continue
        mutated_values = values
        return mutated_values


def create_next_generation(We, Wf, final_PopulationSize, BestFrac, KeepBest, FracGene, FracMut, MUt_PERC, current_generation):
    ancestors = []
    origin = []
    # current_generation = get_current_generation()
    next_gen_folder = f"Generation_{current_generation + 1}"
    os.makedirs(next_gen_folder, exist_ok=True)
    keep_best_count = round(KeepBest * final_PopulationSize)
    BestFracCount = round(BestFrac * final_PopulationSize)
    sorted_list = fit_list(We, Wf, final_PopulationSize)
    # print()
    # print("Sorted list:\n\n")
    # print(sorted_list)
    # print(sorted_list[0][1])
    # print(f"Elites: {sorted_list[:keep_best_count]}")

        # Initialize variables
    elites = []
    fitness_values = set()

    # Iterate through the sorted individuals and select elites
    for individual in sorted_list:
        if individual[1] not in fitness_values:
            elites.append(individual)
            fitness_values.add(individual[1])

        # Stop when you have enough elites
        if len(elites) == keep_best_count:  # Set desired_elites_count to your desired number of elites
            break
        
    # print(elites)

    # import sys
    # sys.exit()
    
    for i in range(len(elites)):
        ancestors.append([elites[i][0], '-'])
        origin.append('Elite')
        old_filename = os.path.join(f"Generation_{current_generation}", elites[i][0] + '.tersoff.modc')
        new_filename = f"{current_generation + 1:03d}-{i+1:03d}" + '.tersoff.modc'
        shutil.copy(old_filename, os.path.join(next_gen_folder, new_filename))

    crossover_count = round(final_PopulationSize * FracGene)
    for i in range(crossover_count):
        parent1, parent2, parentA, parentB = select_random_parents(sorted_list, BestFracCount, current_generation)
        ancestors.append([parentA, parentB])
        origin.append('Crossover')
        child_values = crossover(parent1, parent2)
        new_filename = f"{current_generation + 1:03d}-{keep_best_count+i+1:03d}" + '.tersoff.modc'
        save_potential(next_gen_folder, new_filename, child_values)

    mutation_count = round(final_PopulationSize*FracMut)
    for i in range(mutation_count):
        mutant_canditate, _ , parent_mut, _= select_random_parents(sorted_list, BestFracCount, current_generation)
        # print(mutant_canditate)
        # import sys
        # sys.exit()
        ancestors.append([parent_mut, '-'])
        origin.append('Mutant')
        mutated_values = mutate(mutant_canditate, MUt_PERC, Matrix)
        # print(mutated_values)
        
        new_filename = f"{current_generation + 1:03d}-{keep_best_count+crossover_count+i+1:03d}" + '.tersoff.modc'
        save_potential(next_gen_folder, new_filename, mutated_values)

    
    return ancestors, origin

    # create_next_generation()

def find_origin(ancestors, origin, final_PopulationSize, current_generation):
    # current_generation = get_current_generation()
    file_name = os.path.join('Origins', f"Origins_Gen_{current_generation+1}.txt")
    with open(file_name, 'w') as ff:
        ff.write("Name\t\tOrigin\t\t\tParent 1\tParent 2\n")
        for i in range(1, final_PopulationSize+1):
            ff.write(f"{current_generation+1:03d}-{i:03d}\t\t{origin[i-1]}\t\t{ancestors[i-1][0]}\t\t{ancestors[i-1][1]}\n")

# LAMMPS simulation through python
def RunLAMMPS(N):
    def check_phrase_in_file(file_path):
        with open(file_path, 'r') as file:
            contents = file.read()
            return "All done!" in contents

    lmpDir = f"CalcFold_{N}"
    OriginalDir = os.getcwd()
    lmp_path = os.path.join(OriginalDir, lmpDir)
    # Run LAMMPS using subprocess
    lammps_command = f"lmp -in input_{N}.txt"

    
    subprocess.run(lammps_command, shell=True, check=True, cwd=lmp_path,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    
    Path = os.path.join(OriginalDir, lmpDir, "log.lammps")

    # Loop until the file is found and the phrase is present
    while True:
        if os.path.exists(Path):
            if check_phrase_in_file(Path):
                # print("Phrase 'All done!' found in the file.")
                break
        else:
            continue
            # print("File not found. Waiting...")
        # Wait for 1 second before checking again
        time.sleep(0.1)
    # print(f"Process completed for input_{N}.\n")

def FinalResults(Pop, ListOfMeans, BestFitnessList, ALL_ener, ALL_forc, we, wf):
        OriginData = []
        BestPotentialID = 'None'
        SmallestTotalFitness = float('inf')
        # Get a list of file names and sort them numerically
        folder_path = os.path.join(os.getcwd(), "Origins")
        file_names = os.listdir(folder_path)
        file_names.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
        for fileName in file_names:
            Path = os.path.join(os.getcwd(), "Origins", fileName)
            # print(Path)
            with open(Path, 'r') as file:
                lines = file.readlines()
                for line in lines[1:]:
                    columns = line.strip().split()[:]
                    OriginData.append(columns)

        # print(OriginData)
        with open("Final_Results.txt", 'w') as file:
            gen = 0
            file.write("  ID\t\tOrigin\t\t\tEnergy Fitness\t\t\tForce Fitness\t\t\tTotal Fitness\t\t\tParents\n")
            file.write("  --\t\t-------\t\t\t--------------\t\t\t--------------\t\t\t--------------\t\t\t--------------\n")
            for i in range(len(OriginData)):
                if float(we*ALL_ener[i]+wf*ALL_forc[i]) < SmallestTotalFitness:
                    SmallestTotalFitness = float(we*ALL_ener[i]+wf*ALL_forc[i])
                    BestPotentialID = OriginData[i][0]
                file.write(f"[{OriginData[i][0]}]\t[{OriginData[i][1]}]\t\t{ALL_ener[i]}\t\t{ALL_forc[i]}\t\t{we*ALL_ener[i]+wf*ALL_forc[i]}\t\t{[OriginData[i][2]],[OriginData[i][3]]}\n")
                if i > 0 and (i+1) % Pop == 0:
                    file.write("------------------------------\n")
                    file.write(f"Mean Fitness of Generation: {ListOfMeans[gen]}\n")
                    # file.write(f"Std of Fitness of the Generation: {ListOfDevs[gen]}\n")
                    file.write(f"Best Fitness of Generation: {BestFitnessList[gen]}\n")
                    file.write("------------------------------\n")

                    gen += 1
            file.write(f"Best Potential was: [{BestPotentialID}] with Total fitness: [{SmallestTotalFitness}] \n")
            file.write("Optimization done!")

def read_origins(folder_path = 'Origins'):
    data_list = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r') as file:
                lines = file.readlines()
                lines = lines[1:]
                for line in lines:
                    values = line.strip().split('\t\t')
                    data_list.append(values)

    return data_list


def create_genetic_graph(data):
    genetic_graph = graphviz.Digraph('john', format='png',graph_attr={'size': '100%,100%'}, node_attr={'rankdir': 'TB'})
    genetic_graph.body.append('rankdir=TB')
    for entry in data:
        name, origin, parent1, parent2 = entry
        if origin == "Elite":
            genetic_graph.node(name, label=f'Name: {name}\nOrigin: {origin}', style='filled', fillcolor='royalblue', shape='box')
        elif origin == "Random":
            genetic_graph.node(name, label=f'Name: {name}\nOrigin: {origin}', style='filled', fillcolor='green', shape='box')
        elif origin == "Mutant":
            genetic_graph.node(name, label=f'Name: {name}\nOrigin: {origin}', style='filled', fillcolor='red', shape='box')
        elif origin == "Crossover":
            genetic_graph.node(name, label=f'Name: {name}\nOrigin: {origin}', style='filled', fillcolor='lightgrey', shape='box')
        elif origin == "Seed":
            genetic_graph.node(name, label=f'Name: {name}\nOrigin: {origin}', style='filled', fillcolor='yellow', shape='box')
        else:
            genetic_graph.node(name, label=f'Name: {name}\nOrigin: {origin}', style='filled', fillcolor='white', shape='box')
        if parent1 != '-':
            if parent1 == 'Matrix':
                genetic_graph.edge(parent1, name, style='dashed')
            else:
                genetic_graph.edge(parent1, name)
        if parent2 != '-':
            if parent2 == 'Matrix':
                continue
            else:
                genetic_graph.edge(parent2, name)
    return genetic_graph



Folder_check("Origins")
# Origin for the first generation
file_name = os.path.join('Origins', 'Origins_Gen_1.txt')
with open(file_name, 'w') as origin_1:
    origin_1.write("Name\t\tOrigin\t\t\tParent 1\tParent 2\n")
    for i in range(1, PopulationSize+1):
        
        origin_1.write(f"001-{i:03d}\t\tProtoplast\t\tMatrix\t\tMatrix\n")



# triplets = list(itertools.product(elements, repeat=3))

for g in trange(1, Generations + 1, desc = "Loading", unit="Generation"):
    
    # start_gen = time.time()
    # First we need to create the inputs for the generation
    Folder_check("Inputs")
    write_inputs(g, atom_defs, masses, PotentialType)
    # Create the CalcFolds for the LAMMPS to run
    
    
    combs = list(itertools.product(sorted(os.listdir(f"Generation_{g}"), key=lambda x: int(x.split('-')[-1].split('.')[0])), sorted(os.listdir("Definitions"), key=lambda x: int(x.split('_')[-1].split('.')[0]))))
    Input_files = sorted(os.listdir("Inputs"), key=lambda x: int(x.split('_')[1].split('.')[0]))
    CombinationsOfFiles = [(a, b, c) for (a, b), c in list(zip(combs, Input_files))]
    ParNum = len(Input_files)

    Definitions_files = os.listdir("Definitions")
    how_many_defs = len(Definitions_files)
    EnergyFitness_of_Pots = [[0 for _ in range(how_many_defs)] for _ in range(PopulationSize)]
    ForceFitness_of_Pots = [[0 for _ in range(how_many_defs)] for _ in range(PopulationSize)]
    FinalEnergies = [0 for _ in range(PopulationSize)]   
    FinalForces = [0 for _ in range(PopulationSize)]

    def process_file(i):
        folder_name = f"CalcFold_{i}"
        os.chdir(os.getcwd())
        Folder_check(folder_name)
        copy_file(f"Generation_{g}", f"{CombinationsOfFiles[i-1][0]}", folder_name)
        copy_file("Definitions", f"{CombinationsOfFiles[i-1][1]}", folder_name)
        copy_file("Inputs", f"{CombinationsOfFiles[i-1][2]}", folder_name)
        RunLAMMPS(i)
        CurrentAtomDef = combs[i-1][1]
        CurrentPot = combs[i-1][0]
        # Extract the number using regular expression
        match1 = re.search(r'\d+', CurrentAtomDef)
        if match1:
            BaseNumber = match1.group()
            BaseNumber = int(BaseNumber)
        match2 = re.search(r'-(0*(\d+))\.', CurrentPot)  
        if match2:
            PotNumber = match2.group(2)
            PotNumber = int(PotNumber)

        energyError_i, forceError_i = computeFitness(ExpectedEnergies_from_DataBase, ExpectedForces_from_DataBase,
                                                      getEnergy_fromLog(folder_name), getForces_fromOutput(folder_name)
                                                      , BaseNumber)
        EnergyFitness_of_Pots[PotNumber-1][BaseNumber-1] = energyError_i
        ForceFitness_of_Pots[PotNumber-1][BaseNumber-1] = forceError_i

    with concurrent.futures.ThreadPoolExecutor() as executor1:
        executor1.map(process_file, list(range(1, ParNum+1)))

    
    for i in range(PopulationSize):
        FinalEnergies[i] = math.sqrt(sum(EnergyFitness_of_Pots[i][:]))
        FinalForces[i] = math.sqrt(sum(ForceFitness_of_Pots[i][:]))
        ALL_ener.append(FinalEnergies[i])
        ALL_forc.append(FinalForces[i])

    MeanOfGen = mean([WE*e+WF*f for e, f in zip(FinalEnergies, FinalForces)])
    MeanFitnessOfGen[g-1] = MeanOfGen
    # s_t_d[g-1] = stdev([WE*e+WF*f for e, f in zip(FinalEnergies, FinalForces)])
    # write the result.txt
    best = make_result(WE, WF, g, FinalEnergies, FinalForces)
    # shutil.copy("results.txt", os.path.join("FoldWithResults", f"results_{g}.txt"))
    BestFitnessList[g-1] = best
    FinalResults(PopulationSize, MeanFitnessOfGen[:g], BestFitnessList[:g], ALL_ener, ALL_forc, WE, WF)
    
    if g != Generations:
        ancestors, origin = create_next_generation(WE, WF, PopulationSize, BestFrac, KeepBest, FracGene, FracMut, MUT_PERC, g)
        # print(f"Generation_{g+1} has been created!\n")
        find_origin(ancestors, origin, PopulationSize, g)


FinalResults(PopulationSize, MeanFitnessOfGen, BestFitnessList, ALL_ener, ALL_forc, WE, WF)


# Generate x values (assuming you want a simple index-based x-axis)
x = range(1, len(MeanFitnessOfGen) + 1)
plt.scatter(x, MeanFitnessOfGen, color='red')
plt.xlabel('Number of Generation')
plt.ylabel('Mean Fitness')
plt.title('Mean Fitness per Generation')
# Save the plot to a file (e.g., a PNG file)
plt.savefig('FinalMeans.png')
plt.close()


plt.scatter(x, BestFitnessList)
plt.xlabel('Number of Generation')
plt.ylabel('Best Fitness')
plt.title('Best Fitness value per Generation')
# Save the plot to a file (e.g., a PNG file)
plt.savefig('FinalBests.png')
plt.close()

data = read_origins()

# print(data)
genetic_graph = create_genetic_graph(data)
# genetic_graph.view()

genetic_graph.render(os.path.join(os.getcwd(), 'Genetic_Graph'), format='png', cleanup=True, engine='dot')

print("ALL DONE!")   



    
    










