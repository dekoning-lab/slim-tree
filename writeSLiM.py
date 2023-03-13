#Program to write SLiM Scripts from set functions
#Required packages:
#random
#csv
#numpy
#os
#pandas
#sys

import random, csv, os, pandas, sys
import numpy as np

#Base class from which other classes inherit
class writeSLiM:

    #Initialize required parameters
    def __init__(self, start_para_dict):

        #Set up variables that remain constant for every part of the simulation
        self.general_output_filename = start_para_dict["output_file"]
        self.genome_length = start_para_dict["genome_length"]
        self.fasta_filename = start_para_dict["fasta_filename"]
        self.user_provided_sequence = start_para_dict["user_provided_sequence"]
        self.fitness_profile_calc = start_para_dict["fitness_profile_calc"]
        if(self.user_provided_sequence or (not self.fitness_profile_calc)):
            self.ancestral_sequence = start_para_dict["ancestral_sequence"]

        #Set up the fitness profile and starting distribution of amino acids
        if(self.fitness_profile_calc):
            self.fitness_profile_nums = start_para_dict["fitness_profile_nums"]
            self.fitness_profiles = start_para_dict["fitness_profiles"]
            self.starting_allele_dist = start_para_dict["stationary_distributions"]
            self.amino_acids = start_para_dict["amino_acids"]
            self.min_fitness = str(start_para_dict["min_fitness"])
        else:
            self.dist_pdb_count = start_para_dict["dist_pdb_count"]
            self.main_pdb = os.getcwd() + "/cmaps/main_contact_mat.csv"
            self.distribution_pdbs = os.getcwd() + "/cmaps/distribution_contacts.csv"
            self.max_contacts = start_para_dict["max_contacts"]
            self.max_contact_string = start_para_dict["max_contact_string"]

        #Set up type of model
        self.model_type = start_para_dict["wf_model"]
        if(self.model_type == False and self.fitness_profile_calc):
            self.scaling_factor = start_para_dict["scaling_value"]

        self.coding_regions = start_para_dict["coding_seqs"]

        self.haploidy = start_para_dict["haploidy"]

        #Set up conversion from amino acid to codon
        if(self.fitness_profile_calc):
            with open(os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/slim_codon_nums.csv', newline='') as slim_codon_nums:
                reader = csv.reader(slim_codon_nums)
                slim_codons = list(reader)[1:]

            slim_codon_nums.close()

            self.slim_codon_dict = {}

            for codons in slim_codons :
                amino_acid = codons[2]
                slim_codon_number = int(codons[0])

                if(amino_acid in self.slim_codon_dict.keys()):
                    self.slim_codon_dict[amino_acid].append(slim_codon_number)
                else:
                    self.slim_codon_dict[amino_acid] = [slim_codon_number]



        #Write the initialize function for the SLiM script
    def write_initialize(self, population_parameters):

        initialize_string = ("initialize() {")

        if (self.model_type == False):
            initialize_string += "\n\tinitializeSLiMModelType(\"nonWF\");"

        initialize_string += ("\n\tsetSeed(" + str(random.randint(0,1000000000)) + ");" + "\n\tinitializeSLiMOptions(nucleotideBased=T);")

        #Starting population does not inherit parent sequence, other populations do
        if(population_parameters["parent_pop_name"] == None):

            #Initialize with codons if random, otherwise initialize with sequence given by the user
            if (self.user_provided_sequence or (not self.fitness_profile_calc)):
                initialize_string += "\n\tinitializeAncestralNucleotides(\"" + self.ancestral_sequence + "\");"
            else:
                aa_codon_sequence = str(self.create_codon_seq())
                aa_codon_sequence_str = "c(" + aa_codon_sequence[1: len(aa_codon_sequence) -1] + ")" #Remove brackets and add parentheses
                initialize_string += ("\n\tdefineConstant(\"codons\", " + aa_codon_sequence_str + ");" +
                                "\n\tinitializeAncestralNucleotides(codonsToNucleotides(codons, format=\"char\"));")

        else:
            initialize_string += ("\n\tinitializeAncestralNucleotides(\"" +
                                  population_parameters["parent_pop_name"] + ".fasta\");")

        #If Jukes-Cantor model set mutation rate to Jukes-Cantor model, otherwise set to the mutational matrix
        initialize_string += "\n\tmm = "

        if(population_parameters["jukes_cantor"]):
            initialize_string += "mmJukesCantor(" + str(population_parameters ["mutation_rate"]/3) + ");"
        else:
            initialize_string += population_parameters["mutation_matrix"] +";"


        initialize_string += ("\n\tinitializeMutationTypeNuc(\"m1\", 0.5, \"f\", 0.0);" +
                        "\n\tm1.convertToSubstitution = F;" +
                        "\n\tinitializeGenomicElementType(\"g1\", m1, 1.0, mm);" +
                        "\n\tinitializeGenomicElementType(\"g2\", m1, 1.0, mm);")
        if (self.haploidy == True):
            initialize_string += "\n\tinitializeRecombinationRate(0);"
        else:
            initialize_string += "\n\tinitializeRecombinationRate("+ str(population_parameters ["recombination_rate"])+");"

        #Initialize Genomic Elements according to number of genes for easy visualization in SLiMgui. g1 = coding region, g2 = non-coding region
        if(len(self.coding_regions) == 1):
            initialize_string += "\n\tinitializeGenomicElement(g1, " + str(self.coding_regions[0, 0] * 3) + ", " + str(self.coding_regions[0, 1] * 3 - 1) + ");"
        else:
            for region_num in range(len(self.coding_regions)):

                initialize_string += "\n\tinitializeGenomicElement(g1, " + str(self.coding_regions[region_num, 0] * 3) + ", " + str(self.coding_regions[region_num, 1] * 3 -1) + ");"

                if(region_num == len(self.coding_regions)-1):
                    initialize_string += "\n\tinitializeGenomicElement(g2, " + str((self.coding_regions[region_num, 1] * 3) ) + ", " + str(int(self.genome_length*3)-1) + ");"
                else:
                    initialize_string += "\n\tinitializeGenomicElement(g2, " + str((self.coding_regions[region_num, 1] * 3) ) + ", " + str((self.coding_regions[region_num+1, 0] * 3) -1) + ");"


        initialize_string +="\n}\n\n\n"

        self.output_file.write(initialize_string)



    #Create an initial codon sequence to put into SLiM based on the fitness profile
    def create_codon_seq(self):
        start_codon_nums = self.coding_regions[:,0]
        stop_codon_nums = self.coding_regions[:,1]
        stop_codon_nums = stop_codon_nums - 1

        #Methionine - start codon
        start_codon = 14

        #Stop codons
        stop_codons = [48, 50, 56]

        #Middle codons - chosen according to distribution of alleles
        amino_acids = []

        for dist_num in self.fitness_profile_nums:
                weights = list(self.starting_allele_dist.iloc[:,dist_num])
                amino_acids += random.choices(self.amino_acids, weights = weights, k = 1)


        codons = list(map(self.convert_amino_acid, amino_acids))

        #Replace start and stop codons with start and stop codons
        for codon_num in start_codon_nums: codons[codon_num] = start_codon
        for codon_num in stop_codon_nums: codons[codon_num] = random.choice(stop_codons)

        return (codons)




    #Convert an amino acid to a codon in SLiM by choosing a random available codon
    def convert_amino_acid(self, amino_acid):
        codon_list = self.slim_codon_dict[amino_acid]
        selected_codon = codon_list[random.randint(0,len(codon_list)-1)]

        return selected_codon





    #Write the fitness callback for the SLiM script according to the distribution of fitness effects
    def write_fitness(self):
        #Set up a dictionary in SLiM which takes in amino acids as keys and returns vector of fitnesses
        set_up_fitness = "function (void) setup_fitness(void){\n"

        for key_value in self.fitness_profiles:
            aa_fitnesses = str(self.fitness_profiles[key_value])
            aa_fitnesses = "c(" + aa_fitnesses[1:len(aa_fitnesses)-1] + ")" #Remove first and last parts of strings (the brackets) and add parentheses
            set_up_fitness += "\tsim.setValue(\"" + key_value + "\", " + aa_fitnesses + ");\n"



        #Define required constants
        count = 0;
        profile_num = 0;

        textfile = open("a_file.txt", "w")
        for element in self.fitness_profile_nums:
            textfile.write(str(element) + "\n")
        textfile.close()

        for coding_seq in self.coding_regions:
            fitness_vector = str(self.fitness_profile_nums[coding_seq[0]:(coding_seq[1]+1)])

            fitness_vector = "c(" + fitness_vector[1:len(fitness_vector)-1] + ")"
            set_up_fitness += "\n\tsim.setValue(\"fitness_profiles" + str(count) +"\"," + fitness_vector + ");"
            count += 1


        #List of start and stop codons, remove first and last parts of strings (the brackets) and add parentheses
        start_stop_codons = list(self.coding_regions.flatten()*3)
        #start_stop_codons[-1] = start_stop_codons[-1] -1
        start_stop_codons = str(start_stop_codons)
        print(start_stop_codons)
        if (len(start_stop_codons.split(",")) == 2) :
            value = start_stop_codons.split(",")[1]
            value = value[1:len(value)-1]
            set_up_fitness += ("\n\tdefineConstant(\"start_stop_codon_positions\",matrix(c(0, "+ str(int(value)-3) +"), ncol = 2, byrow = T));\n")
        else:
            set_up_fitness += ("\n\tdefineConstant(\"start_stop_codon_positions\",matrix(c(" +
                        start_stop_codons[1: len(start_stop_codons) -1] + "), ncol = 2, byrow = T));\n")

        #Set up the starting fitnesses
        set_up_fitness += ("\n\tdefineConstant(\"seq_length\", " + str(len(self.fitness_profile_nums)) + ");" +
                        "\n\tget_fitness();\n")



        #Write code to start a fixed state from the starting nucleotide sequence
        set_up_fitness += "\n\tsim.setValue(\"fixations_p1\", sim.chromosome.ancestralNucleotides(format = \"integer\"));"


        #At the start of the sim there are no fixations counted and no non-synonymous or synonymous mutations
        set_up_fitness += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
        set_up_fitness += "\n\tsim.setValue(\"dN_p1\", 0);"
        set_up_fitness += "\n\tsim.setValue(\"dS_p1\", 0);"
        set_up_fitness += "\n}\n\n\n"


        self.output_file.write(set_up_fitness)


        #Defining a function in SLiM which returns the fitness of the ancestral amino acid sequence
        fitness_function_string = ("function (void) get_fitness (void){" +
                                "\n\tposes = start_stop_codon_positions;" +
                                "\n\n\tfor (row_num in (0:(nrow(start_stop_codon_positions)-1))){" +
                                "\n\t\tfitnesses = c();" +
                                "\n\t\taas = codonsToAminoAcids(sim.chromosome.ancestralNucleotides(" +
                                "drop(poses[row_num,0]), drop(poses[row_num, 1])+2, \"codon\"), " +
                                "paste = F);" +
                                "\n\n\t\tsim.setValue(\"ancestral_aa_seq\" + asString(row_num),aas); " +
                                "\n\n\t\tcount = 0;" +
                                "\n\t\tfor (aa in aas){" +
                                "\n\t\t\tfitnesses = c(fitnesses, sim.getValue(aa)[sim.getValue(\"fitness_profiles\" + row_num)[count]]);" +
                                "\n\t\t\tcount = count + 1; \n\t\t}\n" +
                                "\n\t\tsim.setValue(\"ancestral_fitnesses\" + asString(row_num), fitnesses);" +
                                "\n\t\tsim.setValue(\"ancestral_fitness_value\" + asString(row_num), product(fitnesses));"+
                                "\n\t\tsim.setValue(\"ancestral_aas\" + asString(row_num), aas);\n\n\t}\n}\n\n\n")

        self.output_file.write(fitness_function_string)



        #Defining a function in SLiM which returns the fitness of an individual genome

        genome_fitness_function_string = ("function (float) get_genome_fitness (object nucs){" +
                                    "\n\tfitness_value = 1.0;" +
                                    "\n\tfor (row_num in (0:(nrow(start_stop_codon_positions) -1))){" +
                                    "\n\t\tstarting_pos = drop(start_stop_codon_positions[row_num,0]);" +
                                    "\n\t\tending_pos = drop(start_stop_codon_positions[row_num,1])+2;" +
                                    "\n\t\taa_stop_pos = (ending_pos - starting_pos)/3;" +
                                    "\n\t\taa_seq = codonsToAminoAcids(nucs.nucleotides(start = starting_pos, " +
                                    "end = ending_pos, format = \"codon\"), paste = F);" +
                                    "\n\t\tposes = (aa_seq != sim.getValue(\"ancestral_aas\" + row_num));" +
                                    "\n\n\t\tif(sum(poses) == 0){" +
                                    "\n\t\t\tfitness_value = fitness_value * sim.getValue(\"ancestral_fitness_value\" + asString(row_num));" +
                                    "\n\t\t\t next;\n\t\t}" +
                                    "\n\n\t\tfitnesses = sim.getValue(\"ancestral_fitnesses\"+row_num);" +
                                    "\n\t\tfitnesses[poses] = sapply(which(poses), \"sim.getValue(aa_seq[applyValue]);\")" +
                                    "[sim.getValue(\"fitness_profiles\" + row_num)[poses]];")
        if (not self.user_provided_sequence):
            genome_fitness_function_string += ("\n\n\t\tif(any(poses[0] | poses[aa_stop_pos])){" +
                                    "\n\t\t\tfitness_value = fitness_value * product(" + self.min_fitness + "/fitnesses);" +
                                    "\n\t\t\tnext;\n\t\t}")
        genome_fitness_function_string += ("\n\n\t\tif(any(aa_seq[poses] == \"X\")){" +
                                    "\n\t\t\tpos_stop = match(\"X\", aa_seq[0:(length(aa_seq)-1)]);"+
                                    "\n\t\t\tif(pos_stop == 0){fitnesses = " + self.min_fitness + "/fitnesses;}" +
                                    "\n\t\t\telse if (pos_stop + 1 < aa_stop_pos) " +
                                    "{fitnesses[(pos_stop+1):aa_stop_pos] = " + self.min_fitness + "/fitnesses[(pos_stop+1):aa_stop_pos];}\n\t\t}" +
                                    "\n\n\t\tfitness_value = fitness_value * product(fitnesses);\n\t}"+
                                    "\n\n\treturn fitness_value;\n}\n\n\n")

        self.output_file.write(genome_fitness_function_string)




       #Now write out the fitness callback based on the fitness distribution
        if (self.haploidy == False):
            fitness_callback_string = ("fitness(NULL) {return(get_genome_fitness(genome1)" +
                            "*get_genome_fitness(genome2));//If error says total fitness < 0.0, mutation rate is lethal\n}\n\n\n")
        else:
            fitness_callback_string = ("fitness(NULL) {return(get_genome_fitness(genome1));//If error says total fitness < 0.0, mutation rate is lethal\n}\n\n\n")

        self.output_file.write(fitness_callback_string)


    #Writes the fitness functions for protein based contact maps
    def write_fitness_protein_contact(self):
        #Write a void set up fitness function to satisfy existing pipelines for site-heterogeneous
        setup_fitness = "function (void) setup_fitness(void){}\n\n\n"
        self.output_file.write(setup_fitness)

        #Write function to get the fitness of all individuals
        fitness_calc_function = ("function (void) get_fitnesses (No sub_pop, string sub_pop_name) {" +
            "\n\tto_write = codonsToAminoAcids(nucleotidesToCodons(sim.getValue(\"fixations_\" + sub_pop_name)))" +
            "+ codonsToAminoAcids(sub_pop.genomes.nucleotides(format = \"codon\"));" +
            "\n\tfilename = writeTempFile(sub_pop_name, \".txt\", to_write);" +
            "\n\tto_call = \"" + sys.path[0]+"/GetEnergy " +
            str(self.genome_length) + " " + str(self.dist_pdb_count) + " " +  self.distribution_pdbs + " " +
            self.main_pdb + " \" + filename + \" " + str(self.max_contacts) + " " + str(self.max_contact_string) +
            "\";"
            "\n\tfitnesses = asFloat(strsplit(system(to_call), sep = \",\"));" +
            "\n\tsim.setValue(sub_pop_name + \"_fitnesses\", fitnesses);\n}\n\n\n")

        self.output_file.write(fitness_calc_function)

        #Write fitness function to be run for each individual
        fitness_function = ("fitness(NULL){\n\tind = genome2.individual;" +
            "\n\tindex = ind.index * 2;" +
            "\n\tsubpop_id = \"p\" + ind.subpopulation.id;" +
            "\n\treturn product(sim.getValue(subpop_id + \"_fitnesses\")[index :(index + 1)]);\n}\n\n\n")
        self.output_file.write(fitness_function)


    #Write the reproduction callback for non-Wright-Fisher models
    def write_reproduction(self):
        #Basic reproduction callback for now; more functionality could be added later if necessary.
        if (self.haploidy == False):
            reproduction_string = ("reproduction() { " +
                            "\n\tsubpop.addCrossed(individual, subpop.sampleIndividuals(1));" +
                            "\n }\n\n\n")
        else:
            reproduction_string = ("reproduction(){\n\tsubpop.addRecombinant(genome1, NULL, NULL, NULL, NULL, NULL);\n}\n\n\n")

        self.output_file.write(reproduction_string)


    #Write code to count substitutions, make a backup and count generations
    def write_repeated_commands(self, population_parameters):
        #Set up variables for repeated commands
        start_dist = int(population_parameters["dist_from_start"])+1
        end_dist = int(population_parameters["end_dist"])
        pop_name =  population_parameters["pop_name"]

        repeated_commands_string = str(start_dist) +":" + str(end_dist) + "late () {"

        #Write a command to count the substitutions (identity by state) and calculate selection
        if (population_parameters["count_subs"] | population_parameters["calculate_selection"]):
            repeated_commands_string += ("\n\tif(length(sim.mutations)!= 0){"
                        "\n\t\tancestral_genome = sim.getValue(\"fixations_" + pop_name + "\");" +
                        "\n\t\trow_num = " + pop_name + ".individualCount")
            if (self.haploidy):
                repeated_commands_string += ";\n\t\tmuts_mat = integer(row_num*1500);\n\t\tmuts_mat = " + pop_name + ".individuals.genome1.nucleotides(NULL, NULL, \"integer\");"
            else:
                repeated_commands_string += "* 2;\n\t\tmuts_mat = integer(row_num*1500);\n\t\tmuts_mat = " + pop_name + ".genomes.nucleotides(NULL, NULL, \"integer\");"

            repeated_commands_string += ("\n\t\tmuts_mat = matrix(muts_mat, nrow = row_num, byrow = T);" +
                            "\n\t\tcompare_seq = c(muts_mat[0,]);"+
                            "\n\n\t\tfixed_nucs = c(matrixMult(matrix(rep(1, row_num), ncol = " +
                            "row_num), muts_mat)% row_num == 0);" +
                            "\n\n\t\tdifferent_muts = (ancestral_genome != compare_seq);" +
                            "\n\t\tnew_fixations = different_muts & fixed_nucs;" +
                            "\n\t\tsim.setValue(\"fixations_counted_" + pop_name +
                            "\", sim.getValue(\"fixations_counted_" + pop_name+ "\") + sum(new_fixations));" +
                            "\n\n\t\tancestral_genome[new_fixations] = compare_seq[new_fixations];" +
                            "\n\t\tsim.setValue(\"fixations_" + pop_name + "\", ancestral_genome);")
            
            if(population_parameters["calculate_selection"]):
                repeated_commands_string += ("\n\t\tnew_fixations_space = which(new_fixations);" +
                                "\n\t\tfor(fix in new_fixations_space){" +
                                "\n\t\t\tfix_pos = fix % 3;" +
                                "\n\t\t\tif(fix_pos == 0 | fix_pos == 1){" +
                                "\n\t\t\t\told_codon = nucleotidesToCodons(ancestral_genome[(fix-2+2*fix_pos):(fix+2*fix_pos)]);" +
                                "\n\t\t\t\tnew_codon = nucleotidesToCodons(compare_seq[(fix-2):fix]);" +
                                "\n\t\t\t\tif (old_codon == new_codon){" +
                                "\n\t\t\t\t\tsim.setValue(\"dN_" + pop_name +"\", sim.getValue(\"dN" + population_parameters["pop_name"] + "\") + 1);" +
                                "\n\t\t\t\t} else {" +
                                "\n\t\t\t\t\tsim.setValue(\"dS_" + pop_name +"\", sim.getValue(\"dS" + population_parameters["pop_name"] + "\") + 1);" +
                                "\n\t\t\t\t};\n\t\t\t} else {" +
                                "\n\t\t\t\tsim.setValue(\"dN_" + pop_name + "\", sim.getValue(\"dN" + population_parameters["pop_name"] + "\") + 1);" +
                                "\n\t\t\t};\n\t\t};")
            
            
            repeated_commands_string += "\n\t};"

        #Write a command to output when every 100th generation has passed
        if(population_parameters["output_gens"]):
            repeated_commands_string += "\n\n\tif (sim.generation%100 == 0) {\n\t\tcatn(sim.generation);\n\t};"


        #Write a command to write a backup of all individuals after every 100 generations
        if (population_parameters["backup"]):
             repeated_commands_string += ("\n\n\tif (sim.generation%100 == 0) {" +
                        "\n\t\twriteFile(\"" + os.getcwd()+ "/backupFiles/" + pop_name + ".fasta\"," +
                        "(\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));" +
                        "\n\t\tsim.outputFull(\"" + os.getcwd()+ "/backupFiles/" + pop_name + ".txt\");\n\t};")
                        

        repeated_commands_string += "\n}\n\n\n"

        self.output_file.write(repeated_commands_string)



    #Write code to add first population, subpopulation or completely remove population and replace with another
    def write_subpop(self, population_parameters):


        if(population_parameters["parent_pop_name"] == None):
                self.set_up_sim(population_parameters)
        else:
            #Not the starting population, break off from existing population
            define_population_string = (str(int(population_parameters["dist_from_start"])+1) + " { \n" +
                    "\tsim.addSubpopSplit(\""+ population_parameters["pop_name"] + "\"," +
                    str(population_parameters["population_size"]) + ", " + population_parameters["parent_pop_name"]+ ");"+
                    "\n\n\tsim.setValue(\"fixations_" + population_parameters["pop_name"] + "\", sim.getValue(\"fixations_"+
                    population_parameters["parent_pop_name"] +"\"));" +
                    "\n\tsim.setValue(\"fixations_counted_"+ population_parameters["pop_name"]+"\", 0);" +
                    "\n\tsim.setValue(\"dN_"+ population_parameters["parent_pop_name"]+ "\", 0);")

            if(population_parameters["last_child_clade"] == True):
                define_population_string += "\n\t" + population_parameters["parent_pop_name"]+".setSubpopulationSize(0);"

            if (self.haploidy):
                define_population_string += "\n\t" + population_parameters["pop_name"] + ".setCloningRate(1.0);"


            define_population_string += "\n}\n\n\n"

            self.output_file.write(define_population_string)


        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(population_parameters)

        #Write the end of each population
        self.write_end_pop(population_parameters)




    #Write code to add first population, subpopulation or completely remove population and replace with another with non-Wright-Fisher models
    def write_subpop_nonwf(self, population_parameters):
        pop_name = population_parameters["pop_name"]
        if(population_parameters["parent_pop_name"] == None):
                self.set_up_sim(population_parameters)
        else:
            #Not the starting population, break off from existing population
            define_population_string = (str(int(population_parameters["dist_from_start"])) + " late() { \n" +
                                    "\tsim.addSubpop(\"" + pop_name + "\", 0);")

            #If this is the last population broken off, take the remainder of the parent population
            if (population_parameters["last_child_clade"] == True):
                define_population_string += str("\n\tcatn(" + population_parameters["parent_pop_name"] + ".individualCount);"+
                "\n\t" + population_parameters["pop_name"] + ".takeMigrants(" + population_parameters["parent_pop_name"] + ".individuals);" +
                "\n\t" + population_parameters["parent_pop_name"] + ".removeSubpopulation();")
            else:
                #Take proportion of the parent population
                define_population_string += str("\n\tmigrants = sample(" + population_parameters["parent_pop_name"] + ".individuals, asInteger("
                                    + population_parameters["parent_pop_name"] + ".individualCount * " + str(population_parameters["split_ratio"]) + "));\n\t"
                                    + pop_name + ".takeMigrants(migrants);\n\tcatn(" + pop_name + ".individualCount);")

            define_population_string += str("\n\n\tsim.setValue(\"fixations_" + pop_name + "\", sim.getValue(\"fixations_"+
                                    population_parameters["parent_pop_name"] +"\"));" +
                                    "\n\tsim.setValue(\"fixations_counted_"+ pop_name+"\", 0);")
                                    
            define_population_string += "\n\tsim.setValue(\"dN_"+ population_parameters["parent_pop_name"] +"\", 0);"


            define_population_string += "\n}\n\n\n"

            self.output_file.write(define_population_string)

        #Write the early commands - this may need tweaking w/ the fitness algorithm
        early_event = (str(int(population_parameters["dist_from_start"]) + 1) + ":" + str(int(population_parameters["end_dist"])) +
                        " early(){\n\t" + pop_name + ".fitnessScaling = " +
                        str(int(population_parameters["population_size"])) + "/ (" + pop_name +
                        ".individualCount")
        if(self.fitness_profile_calc):
            early_event += ( " * " + str(self.scaling_factor) + ");" )
        else:
            early_event += (");" + "\n\tget_fitnesses(" + pop_name + ", \"" + pop_name + "\");")

        early_event+= "\n}\n\n\n"

        self.output_file.write(early_event)


        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(population_parameters)

        #Write the end of each population
        self.write_end_pop(population_parameters)


    #Set up the simulation by initializing everything
    def set_up_sim(self, population_parameters):
        self.output_file = open(self.general_output_filename + "_" + population_parameters["pop_name"] + ".slim" , "w")

        #Set up the initialize and fitness functions for the new script
        self.write_initialize(population_parameters)
        if(self.fitness_profile_calc):
            self.write_fitness()
        else:
            self.write_fitness_protein_contact()

        #Write reproduction callback if this is a non-WF model
        if (self.model_type == False):
            self.write_reproduction()

        #Make the population and set up fitness effects
        pop_string = ("1 early() {" +
                    "\n\tsetup_fitness();" +
                    "\n\twriteFile(\"" + self.fasta_filename + "_aa.fasta\", \"\", append = F);" +
                    "\n\twriteFile(\"" + self.fasta_filename + "_nuc.fasta\", \"\", append = F);" +
                    "\n\tsim.addSubpop(\"p1\", " + str(population_parameters["population_size"]) + ");")

        if (self.haploidy and self.model_type == True):
            pop_string += "\n\tp1.setCloningRate(1.0);"

        #Write code to start a fixed state from the starting nucleotide sequence
        pop_string += "\n\tsim.setValue(\"fixations_p1\", sim.chromosome.ancestralNucleotides(format = \"integer\"));"

        #At the start of the sim there are no fixations counted
        pop_string += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
        pop_string += "\n}\n\n\n"

        if (self.haploidy and self.model_type == True):
            pop_string += "late(){\n\tsim.subpopulations.individuals.genome2.removeMutations();\n}\n\n\n\n\n"

        self.output_file.write(pop_string)


    #Write the end of a population to save the number of substitutions and output sequence data
    def write_end_pop (self, population_parameters):
        end_population_string = str(int(population_parameters["end_dist"])) + " late() {"

        #If terminal clade output data
        if(population_parameters["terminal_clade"]):
            end_population_string += self.write_terminal_output(population_parameters, pop = population_parameters["pop_name"])

        #Write file with the substitution counts
        if(population_parameters["count_subs"]):
            end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_fixed_mutation_counts.txt\"," +
                "asString(sim.getValue(\"fixations_counted_" + population_parameters["pop_name"] + "\")));" +
                "\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_fixed_mutations.txt\"," +
                " paste(codonsToNucleotides(nucleotidesToCodons(sim.getValue(\"fixations_" + population_parameters["pop_name"] + "\"))), sep = \"\"));")

        #Calculate dN/dS for the population and write into parameters file
        if (self.fitness_profile_calc):
            #end_population_string+= ("\n\tsystem(paste(\"Rscript " + sys.path[0] +"/dNdSCalculations.R\","+ str(population_parameters["population_size"]) +", "+
             #       str(population_parameters["mutation_rate"]) +", \""+ population_parameters["pop_name"] + "\", \""+ sys.path[0]+ "\", \""+ os.getcwd()+ "\", sep = \" \"));" +
              #      "\n\tdNdSFile = readFile(\"" + os.getcwd() + "/"+population_parameters["pop_name"]+"_dNdSDistributions.csv\");\n\tdNdSValues = c();" +
               #     "for (i in 1:(length(sim.getValue(\"X\"))-1)){\n\t\tdNdSValues = c(dNdSValues, asFloat(strsplit(dNdSFile[i], \",\")[1]));}\n\tvalues = c(")

            for i in range(len(self.coding_regions)):
                end_population_string += "sim.getValue(\"fitness_profiles"+ str(i) +"\")[sim.getValue(\"fitness_profiles"+ str(i) +"\") < max(sim.getValue(\"fitness_profiles"+str(i)+"\"))],"
            end_population_string = end_population_string[0:len(end_population_string)-1] #Remove comma at end
            end_population_string += ");\n\twriteFile(\""+ self.fasta_filename +"_parameters.txt\", paste(\"\\n"+ population_parameters["pop_name"] +" estimated dNdS: \", sum(dNdSValues[values])/length(values), sep = \"\"), append = T);"

        #Write out counted dN and dS if required
        if(population_parameters["calculate_selection"]):
            end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_dNdSCounted\"," +
                "paste(\"dN: \", sim.getValue(\"dN_" + population_parameters["pop_name"] + "\")," +
                "\" dS: \", sim.getValue(\"dS_ " + population_parameters["pop_name"] + "\"),  sep = \"\"));" )


        #Write files containing polymorphisms in each population and relative proportions
        if(population_parameters["polymorphisms"]):
            if (self.haploidy):
                end_population_string += ("\n\tpop_seq = sample("+ population_parameters["pop_name"] +".individuals.genome1, 1).nucleotides();\n\tpop_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(pop_seq)), sep = \"\");" +
                            "\n\tpolymorph_str = c();\n\tfixed_str=c();\n\tfor (a in 0:(length(pop_seq)-1)) {\n\t\tdiffs = c();\n\t\tfor (g in " + population_parameters["pop_name"] + ".individuals.genome1.nucleotides()){")

            else:
                end_population_string += ("\n\tpop_seq = sample("+ population_parameters["pop_name"] +".individuals.genomes, 1).nucleotides();\n\tpop_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(pop_seq)), sep = \"\");" +
                            "\n\tpolymorph_str = c();\n\tfixed_str=c();\n\tfor (a in 0:(length(pop_seq)-1)) {\n\t\tdiffs = c();\n\t\tfor (g in " + population_parameters["pop_name"] + ".individuals.genomes.nucleotides()){")

            end_population_string += ("\n\t\t\taa_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(g)), sep = \"\");\n\t\t\tdiffs = c(diffs, aa_seq[a]);\n\t\t}" +
            "\n\t\tunique_diffs = unique(diffs);\n\t\tif (length(unique_diffs) > 1) {\n\t\t\tpolymorph_str = c(polymorph_str, a, \": \");\n\t\t\tfor (p in unique_diffs) {" +
            "\n\t\t\t\tpolymorph_str = c(polymorph_str, p, \": \", length(which(diffs == p)) / length(diffs), \" \");\n\t\t\t}\n\t\tpolymorph_str = c(polymorph_str, \"\\n\");\n\t\t}" +
            " else if (length(unique_diffs) == 1) {\n\t\t\tfixed_str = c(fixed_str, a, \": \", unique_diffs, \"\\n\");\n\t\t}" +
            "\n\t}\n\twriteFile(\"" + os.getcwd() + "/" + population_parameters["pop_name"] + "_polymorphisms.txt\", paste(polymorph_str, sep = \"\"));" +
            "\n\twriteFile(\"" + os.getcwd() + "/" + population_parameters["pop_name"] + "_fixed_sites.txt\", paste(fixed_str, sep = \"\"));")


        if(population_parameters["terminal_clade"] and not self.model_type):
            end_population_string += "\n\t" + population_parameters["pop_name"] + ".removeSubpopulation();"

        end_population_string += "\n}\n\n\n"

        self.output_file.write(end_population_string)






    #Write code to write the output for terminal populations after they have reached their population
    def write_terminal_output(self, population_parameters, pop = "p1"):

        #Set up the names of the 3 fasta files to be output to
        nuc_filename = self.fasta_filename + "_nuc.fasta"
        aa_filename =  self.fasta_filename + "_aa.fasta"
        ancestral_filename = self.fasta_filename + "_fixed.fasta"


        #Set up sampling of the population
        pop_name = population_parameters["pop_name"]
        pop_size = population_parameters["population_size"]
        samp_size = population_parameters["sample_size"]


        #Sample according to number given by user
        terminal_output_string = ""

        #Formulate and output consensus sequence
        if (samp_size == "consensus"):
            terminal_output_string += ("\n\n\tconsensus = \"\";" +
                                    "\n\tfor (i in 0:" + str(self.genome_length * 3 - 1) + "){" +
                                    "\n\t\tconsensus = consensus+ c(\"A\", \"C\", \"G\", \"T\")[whichMax(nucleotideCounts(paste0(matrix(sapply(" + pop_name +
                                    ".genomes.nucleotides(), \"strsplit(applyValue, sep = '');\"), ncol = " + str(self.genome_length * 3) + ", byrow = T)[,i])))];\n\t}" +
                                    "\n\n\tfasta_string_nuc = paste0(\">" + pop_name + ": \\n\", consensus);" +
                                    "\n\twriteFile(\"" + nuc_filename + "\", fasta_string_nuc,append = T);" )
            if (not self.user_provided_sequence):
                terminal_output_string +=("\n\n\tfasta_string_prot = paste0(\">" + pop_name + ": \\n\", codonsToAminoAcids(nucleotidesToCodons(consensus)));" +
                                        "\n\twriteFile(\"" + aa_filename + "\", fasta_string_prot,append = T);")

        else:
            if(samp_size == "all"):
                terminal_output_string += "\n\tgenomes = " + pop + ".genomes;"
            else:
                terminal_output_string += ("\n\tgenomes = sample(" + pop + ".genomes, min(" + str(int(samp_size)) +
                                    ", 2*" + pop + ".individualCount), replace=F);")



            #Iterate through each random sample to write script to output samples of amino acids and nucleotides to fasta files
            terminal_output_string += ("\n\n\tfor (g in genomes){" +
                                        "\n\t\tfasta_string_nuc = paste0(\">\", g.individual, \", " + pop_name + ": \\n\", g.nucleotides());" +
                                        "\n\t\twriteFile(\"" + nuc_filename + "\", fasta_string_nuc,append = T);" )

            if(self.user_provided_sequence):
                terminal_output_string += "}"
            else:
                terminal_output_string += ("\n\t\tfasta_string_prot = paste0(\">\", g.individual, \", " + pop_name + ": \\n\", codonsToAminoAcids(nucleotidesToCodons(g.nucleotides())));" +
                                        "\n\t\twriteFile(\"" + aa_filename + "\", fasta_string_prot,append = T);}" )




        return terminal_output_string





    #Closes the file that has been appended to
    def close_file(self):
        self.output_file.close()
