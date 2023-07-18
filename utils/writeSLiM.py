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
    def __init__(self, start_params):
        self.start_params = start_params


    #Set up the simulation by initializing everything
    def set_up_sim(self, population_parameters):
        self.output_file = open(self.start_params["filenames"][0] + "_" + population_parameters["pop_name"] + ".slim" , "w")
   
        #Set up the initialize and fitness functions for the new script
        self.write_initialize(population_parameters)
        self.write_fitness()
        
        
        #Write reproduction callback if this is a non-WF model
        if (self.start_params["nonWF"]):
            self.write_reproduction()

        #Make the population and set up fitness effects
        pop_string = ("1 early() {" +
                    "\n\tsetup_fitness();" +
                    "\n\twriteFile(\"" + self.start_params["filenames"][1] + "_aa.fasta\", \"\", append = F);" +
                    "\n\twriteFile(\"" + self.start_params["filenames"][1] + "_nuc.fasta\", \"\", append = F);" +
                    "\n\tsim.addSubpop(\"p1\", " + str(population_parameters["population_size"]) + ");" + 
                    "\n\tsim.setValue(\"fixations_p1\", sim.chromosome.ancestralNucleotides(format = \"integer\"));" +    #Write code to start a fixed state from the starting nucleotide sequence
                    "\n\tsim.setValue(\"fixations_counted_p1\", 0);" + #At the start of the sim there are no fixations counted
                    "\n}\n\n\n") #final commands

        self.output_file.write(pop_string)
        


    #Write the initialize function for the SLiM script
    def write_initialize(self, population_parameters):

        initialize_string = ("initialize() {")

        if (self.start_params["nonWF"]):
            initialize_string += "\n\tinitializeSLiMModelType(\"nonWF\");"

        initialize_string += ("\n\tsetSeed(" + str(random.randint(0,1000000000)) + ");" + "\n\tinitializeSLiMOptions(nucleotideBased=T);")

        #Starting population does not inherit parent sequence, other populations do
        if(population_parameters["parent_pop_name"] == None):
        
            #Initialize with codons with ancestral sequence
            initialize_string += ("\n\tinitializeAncestralNucleotides (codonsToNucleotides(c(" + 
                ",".join(self.start_params["ancestral_sequence"]) + "),format=\"char\"));")

        else:
            initialize_string += ("\n\tinitializeAncestralNucleotides(\"" +
                                  population_parameters["parent_pop_name"] + ".fasta\");")

        #If Jukes-Cantor model set mutation rate to Jukes-Cantor model, otherwise set to the mutational matrix
        initialize_string += "\n\tmm = "

        if(self.start_params["jukes_cantor"]):
            initialize_string += "mmJukesCantor(" + str(population_parameters ["mutation_rate"]/3) + ");"
        else:
            initialize_string += population_parameters["mutation_matrix"] +";"


        initialize_string += ("\n\tinitializeMutationTypeNuc(\"m1\", 0.5, \"f\", 0.0);" +
                        "\n\tm1.convertToSubstitution = F;" +
                        "\n\tinitializeGenomicElementType(\"g1\", m1, 1.0, mm);" +
                        "\n\tinitializeGenomicElementType(\"g2\", m1, 1.0, mm);" +
                        "\n\tinitializeRecombinationRate("+ str(population_parameters ["recombination_rate"])+");")

        #Initialize Genomic Elements according to number of genes for easy visualization in SLiMgui. g1 = coding region, g2 = non-coding region
        for region_num in range(len(self.start_params["coding_seqs"])):

            initialize_string += ("\n\tinitializeGenomicElement(g1, " + str(self.start_params["coding_seqs"][region_num][0] * 3) + 
                        ", " + str((self.start_params["coding_seqs"][region_num][1] + 1) * 3 - 1) + ");")

            if(region_num == len(self.start_params["coding_seqs"])-1 and 
                            self.start_params["coding_seqs"][region_num][1] + 1 != self.start_params["genome_length"]):
                initialize_string += ("\n\tinitializeGenomicElement(g2, " + str((self.start_params["coding_seqs"][region_num] [1] + 1) * 3)  + 
                            ", " + str(self.start_params["genome_length"]*3 - 1) + ");")
            elif(self.start_params["coding_seqs"][region_num][1] + 1 != self.start_params["genome_length"]):
                initialize_string += ("\n\tinitializeGenomicElement(g2, " + str((self.start_params["coding_seqs"][region_num][1] + 1) * 3 )  + 
                            ", " + str(self.start_params["coding_seqs"][region_num + 1][0] * 3 - 1) + ");")


        initialize_string +="\n}\n\n\n"

        self.output_file.write(initialize_string)




    #Write the fitness callback for the SLiM script according to the distribution of fitness effects
    def write_fitness(self):
        #Set up a dictionary in SLiM which takes in amino acids as keys and returns vector of fitnesses
        set_up_fitness = "function (void) setup_fitness(void){\n"
        fitness_profiles = self.start_params["fitness_profiles"]

        for key_value in fitness_profiles:
            aa_fitnesses = str(fitness_profiles[key_value])
            aa_fitnesses = "c(" + aa_fitnesses[1:len(aa_fitnesses)-1] + ")" #Remove first and last parts of strings (the brackets) and add parentheses
            set_up_fitness += "\tsim.setValue(\"" + key_value + "\", " + aa_fitnesses + ");\n"



        #Define required constants
        count = 0;
        profile_num = 0;

        textfile = open("a_file.txt", "w")
        for element in self.start_params["fitness_profile_nums"]:
            textfile.write(str(element) + "\n")
        textfile.close()

        for coding_seq in self.start_params["coding_seqs"]:
            fitness_vector = str(self.start_params["fitness_profile_nums"][coding_seq[0]:(coding_seq[1]+1)])

            fitness_vector = "c(" + fitness_vector[1:len(fitness_vector)-1] + ")"
            set_up_fitness += "\n\tsim.setValue(\"fitness_profiles" + str(count) +"\"," + fitness_vector + ");"
            count += 1


        #List of start and stop codons, remove first and last parts of strings (the brackets) and add parentheses
        start_stop_codons = list(self.start_params["coding_seqs"].flatten()*3)
        start_stop_codons = str(start_stop_codons)
        
        if (len(start_stop_codons.split(",")) == 2) :
            value = start_stop_codons.split(",")[1]
            value = value[1:len(value)-1]
            set_up_fitness += ("\n\tdefineConstant(\"start_stop_codon_positions\",matrix(c(0, "+ str(int(value)) +"), ncol = 2, byrow = T));\n")
        else:
            set_up_fitness += ("\n\tdefineConstant(\"start_stop_codon_positions\",matrix(c(" +
                        start_stop_codons[1: len(start_stop_codons)-1] + "), ncol = 2, byrow = T));\n")

        #Set up the starting fitnesses
        set_up_fitness += ("\n\tdefineConstant(\"seq_length\", " + str(self.start_params["genome_length"]*3) + ");" +
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
                                    "\n\t\taa_stop_pos = (ending_pos - starting_pos + 1)/3 - 1;" +
                                    "\n\t\taa_seq = codonsToAminoAcids(nucs.nucleotides(start = starting_pos, " +
                                    "end = ending_pos, format = \"codon\"), paste = F);" +
                                    "\n\t\tposes = (aa_seq != sim.getValue(\"ancestral_aas\" + row_num));" +
                                    "\n\n\t\tif(sum(poses) == 0){" +
                                    "\n\t\t\tfitness_value = fitness_value * sim.getValue(\"ancestral_fitness_value\" + asString(row_num));" +
                                    "\n\t\t\t next;\n\t\t}" +
                                    "\n\n\t\tfitnesses = sim.getValue(\"ancestral_fitnesses\"+row_num);" +
                                    "\n\t\tfitnesses[poses] = sapply(which(poses), \"sim.getValue(aa_seq[applyValue]);\")" +
                                    "[sim.getValue(\"fitness_profiles\" + row_num)[poses]];")
        #Calculate the impact of the stop codon
        if (self.start_params["fasta_file"] ==  None):
            genome_fitness_function_string += ("\n\n\t\tif(any(poses[0] | poses[aa_stop_pos])){" +
                                    "\n\t\t\tfitness_value = fitness_value * product(" + str(self.start_params["min_fitness"]) + "/fitnesses);" +
                                    "\n\t\t\tnext;\n\t\t}")
        genome_fitness_function_string += ("\n\n\t\tif(any(aa_seq[poses] == \"X\")){" +
                                    "\n\t\t\tpos_stop = match(\"X\", aa_seq[0:(length(aa_seq)-1)]);"+
                                    "\n\t\t\tif(pos_stop == 0){fitnesses = " + str(self.start_params["min_fitness"]) + "/fitnesses;}" +
                                    "\n\t\t\telse if (pos_stop + 1 < aa_stop_pos) " +
                                    "{fitnesses[(pos_stop+1):aa_stop_pos] = " + str(self.start_params["min_fitness"]) + "/fitnesses[(pos_stop+1):aa_stop_pos];}\n\t\t}" +
                                    "\n\n\t\tfitness_value = fitness_value * product(fitnesses);\n\t}"+
                                    "\n\n\treturn fitness_value;\n}\n\n\n")

        self.output_file.write(genome_fitness_function_string)




        #Now write out the fitness callback based on the fitness distribution

        fitness_callback_string = ("fitnessEffect() {return(get_genome_fitness(individual.genome1));" + 
                    "//If error says total fitness < 0.0, mutation rate is lethal\n}\n\n\n")

        self.output_file.write(fitness_callback_string)


    #Write the reproduction callback for non-Wright-Fisher models
    def write_reproduction(self):
        #Basic reproduction callback for now; more functionality could be added later if necessary.
        reproduction_string = ("reproduction(){\n\tsubpop.addRecombinant(individual.genome1, NULL, NULL, NULL, NULL, NULL);\n}\n\n\n")

        self.output_file.write(reproduction_string)


    #Write code to count substitutions, make a backup and count generations
    def write_repeated_commands(self, population_parameters, pop_name = None, out = None):
        #Set up variables for repeated commands
        start_dist = int(population_parameters["dist_from_start"])+1
        end_dist = int(population_parameters["end_dist"])
        
        if(pop_name == None):
            pop_name = population_parameters["pop_name"]
        backup_name =  population_parameters["pop_name"]
        
        if(out == None):
            out = self.output_file

        repeated_commands_string = str(start_dist) +":" + str(end_dist) + "late () {"

        #Write a command to count the substitutions (identity by state) and calculate selection
        if (population_parameters["count_subs"] | population_parameters["calculate_selection"]):
            repeated_commands_string += ("\n\tif(length(sim.mutations)!= 0){"
                        "\n\t\tancestral_genome = sim.getValue(\"fixations_" + pop_name + "\");" +
                        "\n\t\trow_num = " + pop_name + ".individualCount* 2;" +
                        "\n\t\tmuts_mat = integer(row_num*1500);\n\t\tmuts_mat = " + pop_name + 
                        ".genomes.nucleotides(NULL, NULL, \"integer\");"+
                        "\n\t\tmuts_mat = matrix(muts_mat, nrow = row_num, byrow = T);" +
                        "\n\t\tcompare_seq = c(muts_mat[0,]);"+
                        "\n\n\t\tfixed_nucs = c(matrixMult(matrix(rep(1, row_num), ncol = " +
                        "row_num), muts_mat)% row_num == 0);" +
                        "\n\n\t\tdifferent_muts = (ancestral_genome != compare_seq);" +
                        "\n\t\tnew_fixations = different_muts & fixed_nucs;" +
                        "\n\n\t\tif(any(new_fixations)){" +
                        "\n\t\t\tnew_fixed = ancestral_genome;"+
                        "\n\t\t\tnew_fixed[new_fixations] = compare_seq[new_fixations];")
                            
                            
            #If there is a flag to also calculate selection - ie. count fixed dN/dS, figure out if dN or dS
            if(population_parameters["calculate_selection"]):
                repeated_commands_string += ("\n\t\t\tnew_fixations_space = which(new_fixations);" +
                                "\n\n\t\t\tdN_name = \"dN_" + pop_name + "\";" +
                                "\n\t\t\tdS_name = \"dS_" + pop_name + "\";"
                                "\n\t\t\tfor(fix in new_fixations_space){" +
                                "\n\t\t\t\tfix_pos = (fix + 1) % 3;" +
                                "\n\t\t\t\tif (fix_pos == 0) {" +
                                "\n\t\t\t\t\told_codon = codonsToAminoAcids(nucleotidesToCodons(ancestral_genome[(fix-2):fix]));" +
                                "\n\t\t\t\t\tnew_codon = codonsToAminoAcids(nucleotidesToCodons(new_fixed[(fix-2):fix]));" +
                                "\n\t\t\t\t\tif (old_codon == new_codon){" +
                                "\n\t\t\t\t\t\tsim.setValue(dS_name, sim.getValue(dS_name) + 1);" +
                                "\n\t\t\t\t\t} else {" +
                                "\n\t\t\t\t\t\tsim.setValue(dN_name, sim.getValue(dN_name) + 1);" +
                                "\n\t\t\t\t\t};\n\t\t\t} else if (fix_pos == 1) {" +
                                "\n\t\t\t\t\told_codon = codonsToAminoAcids(nucleotidesToCodons(ancestral_genome[fix:(fix+2)]));" +
                                "\n\t\t\t\t\tnew_codon = codonsToAminoAcids(nucleotidesToCodons(new_fixed[fix:(fix+2)]));" +
                                "\n\t\t\t\t\tif (old_codon == new_codon){" +
                                "\n\t\t\t\t\t\tsim.setValue(dS_name, sim.getValue(dS_name) + 1);" +
                                "\n\t\t\t\t\t} else {" +
                                "\n\t\t\t\t\t\tsim.setValue(dN_name, sim.getValue(dN_name) + 1);" +
                                "\n\t\t\t\t\t};\n\t\t\t} else {" +
                                "\n\t\t\t\t\tsim.setValue(dN_name, sim.getValue(dN_name) + 1);" +
                                "\n\t\t\t\t};\n\t\t\t};")
            
            #If there is a flag to count substitutions, save fixed substitutions to file
            if(population_parameters["count_subs"]):
                repeated_commands_string += ("\n\n\t\t\tsim.setValue(\"fixations_counted_" + pop_name +
                                "\", sim.getValue(\"fixations_counted_" + pop_name+ "\") + sum(new_fixations));")
                                
            repeated_commands_string += ("\n\t\t\tancestral_genome = new_fixed;" +
                                        "\n\t\t\tsim.setValue(\"fixations_" + pop_name + "\", ancestral_genome);" +
                                        "\n\t\t};\n\t};")

        #Write a command to output when every 100th generation has passed
        if(population_parameters["output_gens"]):
            repeated_commands_string += "\n\n\tif (sim.cycle%100 == 0) {\n\t\tcatn(sim.cycle);\n\t};"


        #Write a command to write a backup of all individuals after every 100 generations
        if (population_parameters["backup"]):
             repeated_commands_string += ("\n\n\tif (sim.cycle%100 == 0) {" +
                        "\n\t\twriteFile(\"" +self.start_params["filenames"][2] +"/" + backup_name + ".fasta\"," +
                        "(\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));" +
                        "\n\t\tsim.outputFull(\"" + self.start_params["filenames"][2]  + "/" + backup_name + ".txt\");\n\t};")
                        

        repeated_commands_string += "\n}\n\n\n"

        out.write(repeated_commands_string)



    #Write code to add first population, subpopulation or completely remove population and replace with another
    def write_subpop(self, population_parameters):

        if(population_parameters["parent_pop_name"] == None):
                self.set_up_sim(population_parameters)
        else:
            # #Not the starting population, break off from existing population
            define_population_string = (str(int(population_parameters["dist_from_start"])+1) + " early() { \n" +
                    "\tsim.addSubpopSplit(\""+ population_parameters["pop_name"] + "\"," +
                    str(population_parameters["population_size"]) + ", " + population_parameters["parent_pop_name"]+ ");"+
                    "\n\n\tsim.setValue(\"fixations_" + population_parameters["pop_name"] + "\", sim.getValue(\"fixations_"+
                    population_parameters["parent_pop_name"] +"\"));" +
                    "\n\tsim.setValue(\"fixations_counted_"+ population_parameters["pop_name"]+"\", 0);" +
                    "\n\tsim.setValue(\"dN_"+ population_parameters["pop_name"]+"\", 0);" +
                    "\n\tsim.setValue(\"dS_"+ population_parameters["pop_name"]+"\", 0);")

            if(population_parameters["last_child_clade"] == True):
                define_population_string += "\n\t" + population_parameters["parent_pop_name"]+".setSubpopulationSize(0);"

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
                define_population_string += str("\n\t" + population_parameters["pop_name"] + ".takeMigrants(" + 
                population_parameters["parent_pop_name"] + ".individuals);" +
                "\n\t" + population_parameters["parent_pop_name"] + ".removeSubpopulation();")
            else:
                #Take proportion of the parent population
                define_population_string += str("\n\tmigrants = sample(" + population_parameters["parent_pop_name"] + ".individuals, asInteger("
                                    + population_parameters["parent_pop_name"] + ".individualCount * " + str(population_parameters["split_ratio"]) + "));\n\t"
                                    + pop_name + ".takeMigrants(migrants);")

            define_population_string += str("\n\n\tsim.setValue(\"fixations_" + pop_name + "\", sim.getValue(\"fixations_"+
                                    population_parameters["parent_pop_name"] +"\"));" +
                                    "\n\tsim.setValue(\"fixations_counted_"+ pop_name+"\", 0);")


            define_population_string += "\n}\n\n\n"

            self.output_file.write(define_population_string)

        #Write the early commands - this may need tweaking w/ the fitness algorithm
        early_event = (str(int(population_parameters["dist_from_start"]) + 1) + ":" + str(int(population_parameters["end_dist"])) +
                        " early(){\n\t" + pop_name + ".fitnessScaling = " +
                        str(int(population_parameters["population_size"])) + "/ (" + pop_name +
                        ".individualCount * " + str(self.start_params["scaling_value"]) + ");" )

        early_event+= "\n}\n\n\n"

        self.output_file.write(early_event)


        #Write the commands that are run for every simulation and the starting population
        self.write_repeated_commands(population_parameters)

        #Write the end of each population
        self.write_end_pop(population_parameters)




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

        #Write file with the number of synonymous and synonymous mutations
        if(population_parameters["calculate_selection"]):
            end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_dNdS_mutations.txt\"," +
                "paste(\"dN: \", sim.getValue(\"dN_" + population_parameters["pop_name"] + "\")/" + str(self.start_params["dn_denom"]) + ", " +
                "\"\\ndS: \", sim.getValue(\"dS_" + population_parameters["pop_name"] + "\")/" + str(self.start_params["ds_denom"]) + 
                ", sep = \"\"));" )


        #Write files containing polymorphisms in each population and relative proportions
        if(population_parameters["polymorphisms"]):
            end_population_string += ("\n\tpop_seq = sample("+ population_parameters["pop_name"] +".individuals.genomes, 1).nucleotides();" +
                                        "\n\tpop_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(pop_seq)), sep = \"\");" +
                                        "\n\tpolymorph_str = c();\n\tfixed_str=c();" +
                                        "\n\tfor (a in 0:(length(pop_seq)-1)) " +
                                        "{\n\t\tdiffs = c();" + 
                                        "\n\t\tfor (g in " + population_parameters["pop_name"] + ".individuals.genomes.nucleotides()){" +
                                        "\n\t\t\taa_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(g)), sep = \"\");" +
                                        "\n\t\t\tdiffs = c(diffs, aa_seq[a]);\n\t\t}" +
                                        "\n\t\tunique_diffs = unique(diffs);" +
                                        "\n\t\tif (length(unique_diffs) > 1) {" +
                                        "\n\t\t\tpolymorph_str = c(polymorph_str, a, \": \");" +
                                        "\n\t\t\tfor (p in unique_diffs) {" +
                                        "\n\t\t\t\tpolymorph_str = c(polymorph_str, p, \": \", length(which(diffs == p)) / length(diffs), \" \");" +
                                        "\n\t\t\t}\n\t\tpolymorph_str = c(polymorph_str, \"\\n\");\n\t\t}" +
                                        " else if (length(unique_diffs) == 1) {" +
                                        "\n\t\t\tfixed_str = c(fixed_str, a, \": \", unique_diffs, \"\\n\");\n\t\t}" +
                                        "\n\t}\n\twriteFile(\"" + os.getcwd() + "/" + population_parameters["pop_name"] + "_polymorphisms.txt\", " +
                                        "paste(polymorph_str, sep = \"\"));" +
                                        "\n\twriteFile(\"" + os.getcwd() + "/" + population_parameters["pop_name"] + 
                                        "_fixed_sites.txt\", paste(fixed_str, sep = \"\"));")


        if(population_parameters["terminal_clade"] and self.start_params["nonWF"]):
            end_population_string += "\n\t" + population_parameters["pop_name"] + ".removeSubpopulation();"

        end_population_string += "\n}\n\n\n"

        self.output_file.write(end_population_string)






    #Write code to write the output for terminal populations after they have reached their population
    def write_terminal_output(self, population_parameters, pop = "p1"):

        #Set up the names of the 3 fasta files to be output to
        nuc_filename = self.start_params["filenames"][1] + "_nuc.fasta"
        aa_filename =  self.start_params["filenames"][1] + "_aa.fasta"
        ancestral_filename = self.start_params["filenames"][1] + "_fixed.fasta"


        #Set up sampling of the population
        pop_name = population_parameters["pop_name"]
        pop_size = population_parameters["population_size"]
        samp_size = population_parameters["sample_size"]


        #Sample according to number given by user
        terminal_output_string = ""

        #Formulate and output consensus sequence
        if (samp_size == "consensus"):
            terminal_output_string += ("\n\n\tconsensus = \"\";" +
                                    "\n\tfor (i in 0:" + str(self.start_params["genome_length"] * 3 - 1) + "){" +
                                    "\n\t\tconsensus = consensus+ c(\"A\", \"C\", \"G\", \"T\")[whichMax(nucleotideCounts(paste0(matrix(sapply(" + pop_name +
                                    ".genomes.nucleotides(), \"strsplit(applyValue, sep = '');\"), ncol = " + str(self.start_params["genome_length"] * 3) + 
                                    ", byrow = T)[,i])))];\n\t}" +
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
                                        "\n\t\twriteFile(\"" + nuc_filename + "\", fasta_string_nuc,append = T);" + 
                                        "\n\t\tfasta_string_prot = paste0(\">\", g.individual, \", " + pop_name + 
                                        ": \\n\", codonsToAminoAcids(nucleotidesToCodons(g.nucleotides())));" +
                                        "\n\t\twriteFile(\"" + aa_filename + "\", fasta_string_prot,append = T);}" )


        return terminal_output_string





    #Closes the file that has been appended to
    def close_file(self):
        self.output_file.close()
