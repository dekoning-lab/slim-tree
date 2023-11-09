#Program to write SLiMTree code for a high performance computing cluster
#Required packages:
#random
#csv
#os


import random, csv, os, sys
from utils.writeSLiM import writeSLiM



class writeSLiMHPC(writeSLiM):

    #Write code to add a new population to the simulation by writing a new script for that population
    def write_subpop(self, population_parameters):
        self.create_scripts(population_parameters)       

        #Set up the initialize and fitness functions for the new script
        super().write_initialize(population_parameters)
        super().write_fitness()

        #Write the commands that are run for every simulation and the starting population
        self.write_start_pop(population_parameters)
        super().write_repeated_commands(population_parameters, pop_name = "p1", out = self.output_file)


        #Finish writing the script
        self.write_end_sim(population_parameters)
        self.output_file.close()




    def write_subpop_nonwf(self, population_parameters):
        self.create_scripts(population_parameters)   
        
        #Set up the initialize and fitness functions for the new script
        super().write_initialize(population_parameters)
        super().write_fitness()

        super().write_reproduction()

        #Write the commands that are run for every simulation and the starting population
        self.write_start_pop(population_parameters)
        super().write_repeated_commands(population_parameters, pop_name = "p1", out = self.output_file)

        self.write_early_function(int(population_parameters["dist_from_start"]) +1, int(population_parameters["end_dist"]), population_parameters)


        #Finish writing the script
        self.write_end_sim(population_parameters)
        self.output_file.close()
        
    
    ##Create a new slim script and batch file for the population        
    def create_scripts(self, population_parameters):
        
        pop_name = population_parameters["pop_name"]
    
        #Write out batch file to run slim script
        batch_file = open(self.start_params["filenames"][0] + "_" + pop_name + ".sh", "w")
        batch_file.write("#!/bin/sh\n\n#SBATCH -J SLiM_Simulation_" + pop_name + "\n#SBATCH -t " + population_parameters["time"] +
                "\n#SBATCH -p "  + population_parameters["partition"] + "\n#SBATCH -o " + pop_name + ".out\n#SBATCH -e " + pop_name +
                ".err\n#SBATCH -n 1" + "\n\nslim \'" + self.start_params["filenames"][0] + "_" + pop_name +".slim\'")
        batch_file.close()
        
        #Open slim script
        self.output_file = open(self.start_params["filenames"][0]+ "_" + pop_name + ".slim" , "w")


    #Write code to set up the starting population for each simulation. If first population, population established, otherwise starting population is loaded
    def write_start_pop(self, population_parameters):

        pop_string = (str(int(population_parameters["dist_from_start"] + 1)) + " late() {" +
                    "\n\tsetup_fitness();")

        #If first population make the population, otherwise load from the parent
        if(population_parameters["parent_pop_name"] == None):
            pop_string += ("\n\twriteFile(\"" + self.start_params["filenames"][1] + "_aa.fasta\", \"\", append = F);" +
                                   "\n\twriteFile(\"" + self.start_params["filenames"][1] + "_nuc.fasta\", \"\", append = F);" +
                                   "\n\tsim.addSubpop(\"p1\", " + str(population_parameters["population_size"]) + ");" )

            #Write code to start a fixed state from the starting nucleotide sequence
            pop_string += "\n\tsim.setValue(\"fixations\", strsplit(sim.chromosome.ancestralNucleotides(),sep = \"\"));"
            
            #Create a second population to dump stuff into at the end of the population simulation for splitting
            if (self.start_params["nonWF"]):
                pop_string += "\n\tsim.addSubpop(\"p2\", 0);"
        else:

            #Set appropriate starting population size
            if (self.start_params["nonWF"]):
                #If a non-WF model, take half of the individuals from the parent population to represent the population split according to tags assigned in previous generation.
                if (population_parameters["last_child_clade"]):
                    #Have population tag 1 have fitness 0.0 so they won't influence next generation
                    pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + "_2.txt\");")
                else:
                    #Have population tag 2 have fitness 0.0 so they won't influence next generation.
                    pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + "_1.txt\");")
            else:
                pop_string += ("\n\tsim.readFromPopulationFile(\"" + population_parameters["parent_pop_name"]  + ".txt\");")
                pop_string += ("\n\tp1.setSubpopulationSize(" + str(population_parameters["population_size"]) + ");")

            #Write code to import in the prevouisly fixed state
            pop_string += ("\n\n\tsim.setValue(\"fixations_p1\", codonsToNucleotides(nucleotidesToCodons(readFile(\""+ 
                            population_parameters["parent_pop_name"] + "_fixed_mutations.txt\")), format = \"integer\"));")


            #Load population into the end of the parent population's script to start this script when parent's finishes
            parent_output_file = open(self.start_params["filenames"][0] + "_" + population_parameters["parent_pop_name"] + ".slim" , "a")
            parent_output_file.write("\n\tsystem(\"sbatch \\\"" + self.start_params["filenames"][0] + "_" + population_parameters["pop_name"] + ".sh\\\"\");")

            if(population_parameters["last_child_clade"]):
                parent_output_file.write("\n}")

            parent_output_file.close()

        #At the start of the sim there are no fixations (synonymous or non-synonymous)
        pop_string += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
        pop_string += "\n\tsim.setValue(\"dN_p1\", 0);"
        pop_string += "\n\tsim.setValue(\"dS_p1\", 0);"
        pop_string += "\n}\n\n\n"

        self.output_file.write(pop_string)


    #Write code for early functions in nonWF models.
    def write_early_function(self, start_dist, end_dist, population_parameters):
        #Write the early commands - this may need tweaking w/ the fitness algorithm
        #Write the early commands - this may need tweaking w/ the fitness algorithm
        pop_name = population_parameters["pop_name"]
        early_event = (str(int(population_parameters["dist_from_start"]) + 2) + ":" + str(int(population_parameters["end_dist"]) + 1) +
                        " early(){\n\t" + "p1.fitnessScaling = " +
                        str(int(population_parameters["population_size"])) + "/" +
                        "p1.individualCount;" +
                        "\n\t p2.fitnessScaling = 0;")

        early_event+= "\n}\n\n\n"

        self.output_file.write(early_event)


    #Write the closing statements to end the simulation and either allow for starting of subsequent simulations or output data (depending on terminal status)
    def write_end_sim(self, population_parameters):

        end_population_string = str(int(population_parameters["end_dist"]) + 1) + " late() {"

        #If terminal clade output data otherwise create data to be loaded into scripts of the clades children
        if(population_parameters["terminal_clade"]):
            end_population_string += super().write_terminal_output(population_parameters, "p1")
        else:
            if (self.start_params["nonWF"]):
                #Tag each individual with either 1 or 2 to go into different subpopulations. Should be split according to proportions.
                end_population_string += ("\n\tp1.individuals.tag = 0;\n\tsample(p1.individuals, asInteger(p1.individualCount* "+
                        str(population_parameters["split_ratio"]) +")).tag = 1;\n\tp1.individuals[p1.individuals.tag == 0].tag = 2;"
                        "\n\tp2.takeMigrants(p1.individuals[p1.individuals.tag == 2]);\n\tsim.outputFull(\""+ 
                        population_parameters["pop_name"] +"_1.txt\");\n\tp1.takeMigrants(p2.individuals);\n\tp2.takeMigrants" +
                        "(p1.individuals[p1.individuals.tag == 1]);"
                        "\n\tsim.outputFull(\""+ population_parameters["pop_name"] +
                        "_2.txt\");\n\tp1.takeMigrants(p2.individuals);\n\tp2.removeSubpopulation();")
            else:
                end_population_string += "\n\tsim.outputFull(\"" + population_parameters["pop_name"] + ".txt\");"

            end_population_string += ("\n\twriteFile(\"" + population_parameters["pop_name"] +
                                      ".fasta\", (\">parent_ancestral_to_load\\n\" + sim.chromosome.ancestralNucleotides()));")
        
        #Write out the fixed mutations - this is different than in single computer because we need the old fixations to run the next pop
        end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_fixed_mutations.txt\"," +
                " paste(codonsToNucleotides(nucleotidesToCodons(sim.getValue(\"fixations_p1\"))), sep = \"\"));")
        
        #Write file with the substitution counts
        if(population_parameters["count_subs"]):
            end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_fixed_mutation_counts.txt\"," +
                "asString(sim.getValue(\"fixations_counted_p1\")));" )

        #Write file with the number of synonymous and synonymous mutations
        if(population_parameters["calculate_selection"]):
            end_population_string += ("\n\twriteFile(\"" + os.getcwd()+ "/" + population_parameters["pop_name"] + "_dNdS.txt\"," +
                "paste(\"dN: \", sim.getValue(\"dN_p1\")/" + str(self.start_params["dn_denom"]) + ", " +
                "\"\\ndS: \", sim.getValue(\"dS_p1\") /" + str(self.start_params["ds_denom"]) +
                ", sep = \"\"));" )



        #Write files containing polymorphisms in each population and relative proportions
        if(population_parameters["polymorphisms"]):
            end_population_string += ("\n\tpop_seq = sample(p1.individuals.genomes, 1).nucleotides();"+
                                    "\n\tpop_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(pop_seq)), sep = \"\");" +
                                    "\n\tpolymorph_str = c();\n\tfixed_str=c();"+
                                    "\n\tfor (a in 0:(length(pop_seq)-1)) {" +
                                    "\n\t\tdiffs = c();\n\t\tfor (g in p1.individuals.genomes.nucleotides()){"+
                                    "\n\t\t\taa_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(g)), sep = \"\");"+
                                    "\n\t\t\tdiffs = c(diffs, aa_seq[a]);\n\t\t}" +
                                    "\n\t\tunique_diffs = unique(diffs);\n\t\tif (length(unique_diffs) > 1) {"+
                                    "\n\t\t\tpolymorph_str = c(polymorph_str, a, \": \");"+
                                    "\n\t\t\tfor (p in unique_diffs) {" +
                                    "\n\t\t\t\tpolymorph_str = c(polymorph_str, p, \": \", length(which(diffs == p)) / length(diffs), \" \");" +
                                    "\n\t\t\t}\n\t\tpolymorph_str = c(polymorph_str, \"\\n\");\n\t\t}" +
                                    " else if (length(unique_diffs) == 1) {" +
                                    "\n\t\t\tfixed_str = c(fixed_str, a, \": \", unique_diffs, \"\\n\");\n\t\t}" +
                                    "\n\t}\n\twriteFile(\"" + os.getcwd() + "/" + population_parameters["pop_name"] + "_polymorphisms.txt\", " +
                                    "paste(polymorph_str, sep = \"\"));" +
                                    "\n\twriteFile(\"" + os.getcwd() + "/" + population_parameters["pop_name"] + 
                                    "_fixed_sites.txt\", paste(fixed_str, sep = \"\"));")


        #If this is the last clade from a certain parent, write script to destroy that parent's temporary files
        if(population_parameters["last_child_clade"]):
            end_population_string += ("\n\n\tsystem(\"rm " + population_parameters["parent_pop_name"] + ".txt\");" +
                                      "\n\tsystem(\"rm " + population_parameters["parent_pop_name"] + ".fasta\");")


        
        

        end_population_string += "\n\tsim.outputFixedMutations();"
        if(population_parameters["terminal_clade"]):
            end_population_string += "\n}"

        self.output_file.write(end_population_string)
