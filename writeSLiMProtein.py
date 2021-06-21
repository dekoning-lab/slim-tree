#For proteins

from writeSLiM import writeSLiM
import random, csv, os, pandas
import numpy as np

class writeSLiMProtein(writeSLiM):


    #Create a codon sequence to be read into SLiM based on the randomly generated contact map.
    def create_codon_seq(self):


        start_codon_nums = self.coding_regions[:,0]
        stop_codon_nums = self.coding_regions[:,1]

        #Methionine - start codon
        start_codon = 14

        #Stop codons
        stop_codons = [48, 50, 56]

        #Middle codons - chosen according to distribution of alleles
        amino_acids = []

        #Open the contact map, and make amino acids according to size.
        contactmap = self.contact_map

        for dist_num in self.fitness_profile_nums:
                weights = list(self.starting_allele_dist.iloc[:,dist_num])
                amino_acids += random.choices(self.amino_acids, weights = weights, k = 1)


        codons = list(map(self.convert_amino_acid, amino_acids))

        #Replace start and stop codons with start and stop codons - only one
        codons[0] = start_codon
        codons[len(codons) - 1] = random.choice(stop_codons)

        # for codon_num in start_codon_nums: codons[codon_num] = start_codon
        # for codon_num in stop_codon_nums: codons[codon_num] = random.choice(stop_codons)


        return (codons)


    #Write fitness based on protein structure.
    def write_fitness(self):
        #also stuff here

        #Set up a dictionary in SLiM which takes in amino acids as keys and returns vector of fitnesses
        set_up_fitness = "function (void) setup_fitness(void){\n"

        #Map amino acids to indexes in the Miyazawa-Jernigan matrix
        set_up_fitness += "\n\tsim.setValue(\"C\", 0);\n\tsim.setValue(\"M\", 1);\n\tsim.setValue(\"F\", 2);\n\tsim.setValue(\"I\", 3);\n\tsim.setValue(\"L\", 4);\n\tsim.setValue(\"V\", 5);"
        set_up_fitness += "\n\tsim.setValue(\"W\", 6);\n\tsim.setValue(\"Y\", 7);\n\tsim.setValue(\"A\", 8);\n\tsim.setValue(\"G\", 9);\n\tsim.setValue(\"T\", 10);\n\tsim.setValue(\"S\", 11);\n\tsim.setValue(\"Q\", 12);"
        set_up_fitness += "\n\tsim.setValue(\"N\", 13);\n\tsim.setValue(\"E\", 14);\n\tsim.setValue(\"D\", 15);\n\tsim.setValue(\"H\", 16);\n\tsim.setValue(\"R\", 17);\n\tsim.setValue(\"K\", 18);\n\tsim.setValue(\"P\", 19);"


        #Read in the Miyazawa-Jerigan matrix
        set_up_fitness += "\n\n\tmj_matrix_file = readFile(\"fitnessDataFiles/mj_test_matrix.csv\");"

        #Convert into a matrix in SLiM
        set_up_fitness += ("\n\tmj_matrix_temp = matrix(rnorm(400, 0, 1), nrow = 20, ncol = 20);" +
                           "\n\tfor (i in 0:19) {\n\t\tnums = strsplit(mj_matrix_file[i], \",\");" +
                           "\n\t\tfor (j in 0:19) {\n\t\t\t	mj_matrix_temp[i,j] = asFloat(nums[j]);\n\t\t}\n\t}"+
                           "\n\tdefineConstant(\"mj_matrix\", mj_matrix_temp);")


        #Define the contacts. Read in the contact map and then write it into the file.
        map = self.contact_map

        contacts = [] #List of contacts
        i = 0
        j = 1
        while (i < len(map)):
            while (j < len(map)):
                if (map[i][j]):
                    contacts.append([i,j])
                j+=1
            i+=1
            j = i+1

        set_up_fitness += "\n\tdefineConstant(\"contacts\", c("

        for c in contacts:
            set_up_fitness += (str(c[0]) + ", "+ str(c[1]) + ", ")

        set_up_fitness = set_up_fitness[0:len(set_up_fitness)-2]

        set_up_fitness += "));"

        set_up_fitness += ("\n\tdefineConstant(\"contact_matrix\", matrix(contacts, nrow = 2));");


        #Define required constants
        set_up_fitness += ("\n\tdefineConstant(\"seq_length\", " + str(len(map) +1) + ");" +
                           "\n\tdefineConstant(\"ancestral_aa_seq\", strsplit(codonsToAminoAcids(codons), sep = \"\"));" +
                           "\n\tdefineConstant(\"ancestral_fitnesses\", get_fitness(codonsToAminoAcids(codons)));" +
                           "\n\tdefineConstant(\"ancestral_fitness_value\", product(ancestral_fitnesses));" )

        start_stop_codons = str(list(self.coding_regions.flatten()))
        set_up_fitness += "\n\tdefineConstant(\"start_stop_codon_positions\",c(" + start_stop_codons[1: len(start_stop_codons) -1] + "));"

        #Write code to start a fixed state from the starting nucleotide sequence
        set_up_fitness += "\n\tsim.setValue(\"fixations_p1\", strsplit(sim.chromosome.ancestralNucleotides(),sep = \"\"));"


        #At the start of the sim there are no fixations counted
        set_up_fitness += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
        set_up_fitness += "\n}\n\n\n"


        self.output_file.write(set_up_fitness)


        #Defining a function in SLiM which returns the fitness of the ancestral amino acid sequence

        fitness_function_string = ("function (float) get_fitness (string aa_seq_string) {" +
                                "\n\taa_seq = strsplit(aa_seq_string, sep = \"\");" +
                                "\n\tfitnesses = c();" +
                                "\n\tc=0;" +
                                "\n\twhile (c < length(contact_matrix[0,])) {" +
                                "\n\t\tif (aa_seq[drop(contact_matrix[0, c])] == \"X\") {\n\t\t\t return fitnesses;\n\t\t}" +
                                "\n\t\tif (aa_seq[drop(contact_matrix[1, c])] != \"X\") {" +
                                "\n\t\t\tnew_fitness = (mj_matrix[sim.getValue(aa_seq[drop(contact_matrix[0, c])]), sim.getValue(aa_seq[drop(contact_matrix[1, c])])]*-1) + 1;" +
                                "\n\t\t\tfitnesses = c(fitnesses, new_fitness);\n\t\t}" +
                                "\n\t\telse {\n\t\t\t fitnesses = c(fitnesses, 1);\n\t\t}\n\t\tc = c+1;\n\t}" +
                                "\n\treturn fitnesses;\n}\n\n\n\n\n")


        self.output_file.write(fitness_function_string)


        genome_fitness_function_string = ("function (float) get_genome_fitness (string aa_seq_string, string population_aa_seq, float population_fitnesses) {" +
                                        "\n\taa_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(aa_seq_string)), sep = \"\");" +
                                        "\n\tpop_aa_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(population_aa_seq)), sep = \"\");" +
                                        "\n\tposes = which(aa_seq != pop_aa_seq);" +
                                        "\n\tif (length(poses) == 0) {\n\t\treturn product(population_fitnesses);\n\t}" +
                                        "\n\tseqs = aa_seq[poses];\n\tfitnesses = population_fitnesses;" +
                                        "\n\tif (aa_seq[length(aa_seq)-1] != \"X\"){\n\t\treturn 0.0;\n\t}" +
                                        "\n\tif (aa_seq[0] != \"M\") {\n\t\treturn 0.0;\n\t}" +
                                        "\n\tfor (p in poses) {" +
                                        "\n\t\tif (aa_seq[p] == \"X\") {" +
                                        "\n\t\t\tfitnesses[which(contact_matrix[0, ] >= p)] = 0.1/fitnesses[which(contact_matrix[0, ] >= p)];" +
                                        "\n\t\t\treturn(product(fitnesses));\n\t\t}" +
                                        "\n\t\tinit_contacts = contact_matrix[1, which(contact_matrix[0, ] == p)];" +
                                        "\n\t\tif (length(init_contacts) != 0) {" +
                                        "\n\t\t\tfor (i in init_contacts) {\n\t\t\t\tif (aa_seq[i] != \"X\") {" +
                                        "\n\t\t\t\t\tnew_fitness = (mj_matrix[sim.getValue(aa_seq[p]), sim.getValue(aa_seq[i])] * -1) + 1;" +
                                        "\n\t\t\t\t\tfitnesses[setIntersection(which(contact_matrix[0,] == p), which(contact_matrix[1,] == i))] = new_fitness;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t}" +
                                        "\n\t\tsecond_contacts = contact_matrix[0, which(contact_matrix[1, ] == p)];" +
                                        "\n\t\tif (length(second_contacts) != 0) {\n\t\t\tfor (s in second_contacts) {" +
                                        "\n\t\t\t\tif (!any(poses[] == s)) {\n\t\t\t\tif (aa_seq[s] != \"X\") {" +
                                        "\n\t\t\t\t\t\tnew_fitness = (mj_matrix[sim.getValue(aa_seq[p]), sim.getValue(aa_seq[s])] *-1) + 1;" +
                                        "\n\t\t\t\t\t\tfitnesses[setIntersection(which(contact_matrix[1,] == p), which(contact_matrix[0,] == s))] = new_fitness;" +
                                        "\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t}" +
                                        "\n\t\tstop_pos = which(aa_seq[] == \"X\")[0];"+
                                        "\n\t\tif (stop_pos != length(aa_seq) - 1) {\n\t\t\tfitnesses[which(contact_matrix[0, ] >= stop_pos)] = 0.1/fitnesses[which(contact_matrix[0, ] >= stop_pos)];\n\t\t}" +
                                        "\n\treturn product(fitnesses);\n}\n\n\n\n")

        self.output_file.write(genome_fitness_function_string)


       #Now write out the fitness callback based on the fitness distribution
        fitness_callback_string = ("fitness(NULL) {" +
                                   "\n\tind_fitness = get_genome_fitness(genome1.nucleotides(), subpop.getValue(\"population_aa_seq\"), subpop.getValue(\"population_fitness\")) * get_genome_fitness(genome2.nucleotides(), subpop.getValue(\"population_aa_seq\"), subpop.getValue(\"population_fitness\"));" +
                                   "\n\tvalue = ind_fitness / subpop.getValue(\"mean\");\n\treturn value;\n}\n\n\n\n")

        self.output_file.write(fitness_callback_string)

    #A function to write an early function needed for protein scaling, then closing the file.
    def close_file(self):
        scaling_str = ("early() {\n\tfor (s in sim.subpopulations) {" +
		              "\n\t\tpop_seq = sample(s.individuals.genomes, 1).nucleotides();"+
                      "\n\t\tif (all(pop_seq == s.individuals.genomes.nucleotides())){" +
                      "\n\t\t\ts.setValue(\"population_aa_seq\", pop_seq);" +
                      "\n\t\t\ts.setValue(\"population_fitness\", get_fitness(codonsToAminoAcids(nucleotidesToCodons(pop_seq))));\n\t\t}" +
                      "\n\t\tif (isNULL(s.getValue(\"population_aa_seq\"))) {" +
                      "\n\t\t\ts.setValue(\"population_aa_seq\", pop_seq);" +
                      "\n\t\t\ts.setValue(\"population_fitness\", get_fitness(codonsToAminoAcids(nucleotidesToCodons(pop_seq))));\n\t\t}" +
                      "\n\t\tfitnesses = c();\n\t\tfor (g in s.individuals.genomes) {" +
                      "\n\t\t\tfitnesses = c(fitnesses, get_genome_fitness(g.nucleotides(), s.getValue(\"population_aa_seq\"), s.getValue(\"population_fitness\")));\n\t\t}" +
                      "\n\t\ts.setValue(\"mean\", mean(fitnesses)^2);\n\t\ts.setValue(\"sd\", sd(fitnesses));\n\t}\n}\n\n\n")

        self.output_file.write(scaling_str)
        self.output_file.close()
