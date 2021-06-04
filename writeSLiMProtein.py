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

        set_up_fitness = set_up_fitness[0:len(set_up_fitness)-3]

        set_up_fitness += "));"


        #Define required constants
        set_up_fitness += ("\n\tdefineConstant(\"seq_length\", " + str(len(map)) + ");" +
                           "\n\tdefineConstant(\"ancestral_aa_seq\", strsplit(codonsToAminoAcids(codons), sep = \"\"));" +
                           "\n\tdefineConstant(\"ancestral_fitnesses\", get_fitness(codonsToAminoAcids(codons)));" +
                           "\n\tdefineConstant(\"ancestral_fitness_value\", product(ancestral_fitnesses));" )

        start_stop_codons = str(list(self.coding_regions.flatten()))
        set_up_fitness += "\n\tdefineConstant(\"start_stop_codon_positions\",c(" + start_stop_codons[1: len(start_stop_codons) -1] + "));" #Remove first and last parts of strings (the brackets) and add parentheses

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
                                "\n\twhile (c < length(contacts)) {" +
                                "\n\t\tif (aa_seq[contacts[c]] == \"X\") {\n\t\t\t break;\n\t\t}" +
                                "\n\t\tif (aa_seq[contacts[c+1]] != \"X\") {" +
                                "\n\t\t\tnew_fitness = mj_matrix[sim.getValue(aa_seq[contacts[c]]), sim.getValue(aa_seq[contacts[c+1]])] + 1;" +
                                "\n\t\t\tfitnesses = c(fitnesses, new_fitness);\n\t\t}" +
                                "\n\t\tc = c+2;\n\t}" +
                                "\n\treturn product(fitnesses);\n}\n\n\n\n\n")


        self.output_file.write(fitness_function_string)





       #Now write out the fitness callback based on the fitness distribution
        fitness_callback_string = ("fitness(NULL) {" +
                                   "\n\tfor (g in individual.genomes){" +
                                   "\n\t\taa_seq = codonsToAminoAcids(nucleotidesToCodons(g.nucleotides()));"
                                   "\n\t\tind_fitness = get_fitness(aa_seq);" +
                                   "\n\t\tvalue = ind_fitness / subpop.getValue(\"mean\");" +
                                   "\n\t\treturn value;\n\t}\n}\n\n\n\n")

        self.output_file.write(fitness_callback_string)

    #A function to write an early function needed for protein scaling, then closing the file.
    def close_file(self):
        scaling_str = ("early() {" +
                    "\n\tfor (s in sim.subpopulations) {" +
                    "\n\t\tfitnesses = c();" +
                    "\n\t\tfor (g in s.individuals.genomes) {" +
                    "\n\t\t\tfitnesses = c(fitnesses, get_fitness(codonsToAminoAcids(nucleotidesToCodons(g.nucleotides()))));\n\t\t}" +
                    "\n\t\ts.setValue(\"mean\", mean(fitnesses));" +
                    "\n\t\ts.setValue(\"sd\", sd(fitnesses));\n\t}\n}\n\n\n")

        self.output_file.write(scaling_str)
        self.output_file.close()
