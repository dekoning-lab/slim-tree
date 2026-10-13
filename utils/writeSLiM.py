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
from utils import nonWFLifeCycle

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
        self.write_fitness(population_parameters)
        
        
        #Write reproduction callback if this is a non-WF model
        if (self.start_params["nonWF"]):
            self.write_reproduction()

        #Make the population and set up fitness effects
        pop_string = ("1 early() {" +
                    "\n\tsetup_fitness();" +
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
        if(population_parameters["jukes_cantor"]):
            initialize_string += "mmJukesCantor(" + str(population_parameters ["mutation_rate"]/3) + ");"
        else:
            initialize_string += population_parameters["mutation_matrix"][1] +";"


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
    def write_fitness(self, population_parameters):
        
        set_up_fitness = "function (void) setup_fitness(void){\n"
    
        if (not self.start_params["neutral_evolution"]):
            #Set up a dictionary in SLiM which takes in amino acids as keys and returns vector of fitnesses
            fitness_profiles = self.start_params["fitness_profiles"]
    
            for key_value in fitness_profiles:
                aa_fitnesses = str(fitness_profiles[key_value])
                aa_fitnesses = "c(" + aa_fitnesses[1:len(aa_fitnesses)-1] + ")" #Remove first and last parts of strings (the brackets) and add parentheses
                set_up_fitness += "\tsim.setValue(\"" + key_value + "\", " + aa_fitnesses + ");\n"
    
    
    
            #Define required constants
            count = 0;
            profile_num = 0;
    
            textfile = open(os.path.join(os.path.dirname(self.start_params["filenames"][1]), "fitness_profile_nums.txt"), "w")
            for element in self.start_params["fitness_profile_nums"]:
                textfile.write(str(element) + "\n")
            textfile.close()
    
            #Write out fitness profiles
            fitness_vector = str(population_parameters["fitness_profile_nums"])
            fitness_vector = "c(" + fitness_vector[1:len(fitness_vector)-1] + ")"
            set_up_fitness += "\n\tsim.setValue(\"fitness_profiles\"," + fitness_vector + ");"
    
    
    
            #List of start and stop codons, remove first and last parts of strings (the brackets) and add parentheses
            start_stop_codons = np.ravel(self.start_params["coding_seqs"]).tolist()
            start_stop_codons = str(start_stop_codons)
            set_up_fitness += ("\n\tdefineConstant(\"start_stop_codon_positions\",c(" + start_stop_codons[1:len(start_stop_codons)-1] +"));\n")
    
            #Set up the starting fitnesses
            set_up_fitness += ("\n\tdefineConstant(\"seq_length\", " + str(self.start_params["genome_length"]*3) + ");" +
                            "\n\tget_fitness();\n")

            #Fitness values are looked up in one flat table: FIT_TABLE[aa_index * N_PROF + profile]. With the site stop model
            #(the default) two extra profiles, the last two columns, score each gene's start codon (ATG optimal) and terminal
            #stop codon (any stop optimal); every other state at those sites costs the stop penalty (stop_penalty_delta)
            aa_letters = list(fitness_profiles.keys())
            columns = {aa: [float(v) for v in fitness_profiles[aa]] for aa in aa_letters}
            if (self.stop_model() == "site"):
                delta = self.stop_penalty_delta()
                for aa in aa_letters:
                    columns[aa] = columns[aa] + [1.0 if aa == "M" else 1.0 - delta, 1.0 if aa == "X" else 1.0 - delta]
            fit_table = [v for aa in aa_letters for v in columns[aa]]
            set_up_fitness += ("\n\tdefineConstant(\"AA_LETTERS\", c(" + ", ".join("\"" + aa + "\"" for aa in aa_letters) + "));" +
                            "\n\tdefineConstant(\"N_PROF\", " + str(len(columns[aa_letters[0]])) + ");" +
                            "\n\tdefineConstant(\"FIT_TABLE\", c(" + ", ".join(repr(v) for v in fit_table) + "));\n")
            set_up_fitness += self.fitness_constants()
            if (self.start_params["nonWF"]):
                set_up_fitness += "\tdefineConstant(\"NONWF_SELECTION_SCALE\", " + repr(nonWFLifeCycle.selection_scale()) + ");\n"



        #Write code to start a fixed state from the starting nucleotide sequence
        set_up_fitness += "\n\tsim.setValue(\"fixations_p1\", sim.chromosome.ancestralNucleotides(format = \"integer\"));"


        #At the start of the sim there are no fixations counted and no non-synonymous or synonymous mutations
        if (self.start_params["count_subs"] | self.start_params["calculate_selection"]):
            set_up_fitness += "\n\tsim.setValue(\"fixations_counted_p1\", 0);"
            
            if(self.start_params["calculate_selection"]):
                set_up_fitness += "\n\tsim.setValue(\"dN_p1\", 0);"
                set_up_fitness += "\n\tsim.setValue(\"dS_p1\", 0);"
                
            set_up_fitness += "\n\tsim.setValue(\"subs_p1\", \"\\n\\nSubstitutions:\");"
        set_up_fitness += "\n}\n\n\n"


        self.output_file.write(set_up_fitness)

        if (not self.start_params["neutral_evolution"]):
            #Defining a function in SLiM which returns the fitness of the ancestral amino acid sequence
            fitness_function_string = ("function (void) get_fitness (void){" +
                                    "\n\tfitnesses = c();" +
                                    "\n\taas = codonsToAminoAcids(sim.chromosome.ancestralNucleotides(format = \"codon\"), " +
                                    "paste = F);" +
                                    "\n\n\tsim.setValue(\"ancestral_aa_seq\",aas); " +
                                    "\n\n\tcount = 0;" +
                                    "\n\tfor (aa in aas){" +
                                    "\n\t\tfitnesses = c(fitnesses, sim.getValue(aa)[sim.getValue(\"fitness_profiles\")[count]]);" +
                                    "\n\t\tcount = count + 1; \n\t\t}\n" +
                                    "\n\tsim.setValue(\"ancestral_fitnesses\", fitnesses);" +
                                    "\n\tsim.setValue(\"ancestral_fitness_value\", sum(fitnesses));"+
                                    "\n\tsim.setValue(\"ancestral_aas\", aas);\n}\n\n\n")
    
            self.output_file.write(fitness_function_string)
    
    
    
            #Defining a function in SLiM which returns the fitness of an individual genome

            #Fitness is multiplicative across sites. The fitness profiles are fitted (utils/calculateFitnesses.py) treating the
            #difference between two amino acids' values as the selection coefficient of that change at a single site, so each
            #codon whose amino acid differs from the ancestor's multiplies a haplosome's selection factor by
            #exp(f_new - f_ancestral). Profile values are therefore log-fitness contributions. get_genome_fitness works on
            #integer codons and amino-acid indices (CODON_TO_AA) rather than amino-acid strings, which is about 4x faster
            #when few codons differ from the ancestor.
            #Site stop model (default): a stop codon is one more state at each site, scored from the table like any amino
            #acid (its value is set in findFitness.apply_stop_model); get_genome_fitness returns the selection factor.
            #Gene stop model: in addition, a gene is a null allele in a haplosome if its start or terminal stop codon has
            #changed or it has a premature stop, and each gene multiplies fitness by NULL_FACTORS[number of functional
            #copies]; get_genome_fitness returns c(selection factor, functional flag of gene 0, gene 1, ...). A null copy's
            #codons stay in the selection factor, so a knockout can never raise fitness
            gene_model = self.stop_model() == "gene"
            genome_fitness_function_string = ("function (float) get_genome_fitness (object nucs){" +
                                    "\n\tcodons = nucs.nucleotides(format = \"codon\");" +
                                    "\n\tposes = which(codons != ANC_CODONS);" +
                                    "\n\tif (size(poses) == 0) return " + ("c(1.0, ALL_FUNCTIONAL)" if gene_model else "1.0") + ";" +
                                    "\n\taa = CODON_TO_AA[codons[poses]];" +
                                    "\n\tchanged = aa != ANC_AA[poses];" +
                                    "\n\tif (!any(changed)) return " + ("c(1.0, ALL_FUNCTIONAL)" if gene_model else "1.0") + ";" +
                                    "\n\tposes = poses[changed];" +
                                    "\n\taa = aa[changed];")
            if (gene_model):
                genome_fitness_function_string += ("\n\tfunctional = ALL_FUNCTIONAL;" +
                                    "\n\tnull_genes = c(NULL_IF_CHANGED[poses], NULL_IF_STOP[poses[aa == X_INDEX]]);" +
                                    "\n\tnull_genes = null_genes[null_genes >= 0];" +
                                    "\n\tif (size(null_genes) > 0) functional[null_genes] = 0.0;" +
                                    "\n\treturn c(exp(sum(FIT_TABLE[aa * N_PROF + PROFILE_OF_CODON[poses]] - ANC_FIT[poses])), functional);\n}\n\n\n")
            else:
                genome_fitness_function_string += "\n\treturn exp(sum(FIT_TABLE[aa * N_PROF + PROFILE_OF_CODON[poses]] - ANC_FIT[poses]));\n}\n\n\n"

            self.output_file.write(genome_fitness_function_string)

            if (not self.start_params["nonWF"]):
                if (gene_model):
                    self.output_file.write("fitnessEffect() {" +
                                           "\n\tw1 = get_genome_fitness(individual.haploidGenome1);" +
                                           "\n\tw2 = get_genome_fitness(individual.haploidGenome2);" +
                                           "\n\treturn w1[0] * w2[0] * product(NULL_FACTORS[asInteger(w1[1:N_GENES] + w2[1:N_GENES])]);\n}\n\n\n")
                else:
                    self.output_file.write("fitnessEffect() {" +
                                           "\n\treturn get_genome_fitness(individual.haploidGenome1) * get_genome_fitness(individual.haploidGenome2);\n}\n\n\n")
            else:
                self.output_file.write(self.nonwf_fitness_function())


    #Stop model and penalty (see utils/findFitness.apply_stop_model and README "Stop codons")
    def stop_model(self):
        return self.start_params.get("stop_model", "site")

    #Per-copy log-fitness gap of the site stop model: S_stop / (4 N_ref - 2), with N_ref the population size the
    #profiles were fitted at (the same convention as calculateFitnesses.ratio_func)
    def stop_penalty_delta(self):
        return float(self.start_params.get("stop_penalty", 100.0)) / (4.0 * self.start_params["population_size"] - 2.0)

    #Eidos constants for fitness (written inside setup_fitness). Per gene: GENE_STARTS/GENE_ENDS (codon indices of its
    #start and terminal stop codon). Site stop model: those codons use the two boundary profiles (the last two columns of
    #FIT_TABLE). Gene stop model: per codon (-1 = none) GENE_OF_CODON, NULL_IF_CHANGED (the gene whose start or stop codon
    #this is) and NULL_IF_STOP (the gene that a stop codon here knocks out: any codon strictly between its start and stop
    #codons), and NULL_FACTORS (fitness factor for 0, 1 or 2 functional copies). Then the integer lookups used by
    #get_genome_fitness
    def fitness_constants(self):
        coding = np.asarray(self.start_params["coding_seqs"]).reshape(-1, 2)
        starts = [int(s) for s in coding[:, 0]]
        ends = [int(e) for e in coding[:, 1]]

        def eidos_vector(values):
            return "c(" + ", ".join(repr(v) for v in values) + ")"

        out = ("\n\tdefineConstant(\"N_GENES\", " + str(len(starts)) + ");" +
               "\n\tdefineConstant(\"GENE_STARTS\", " + eidos_vector(starts) + ");" +
               "\n\tdefineConstant(\"GENE_ENDS\", " + eidos_vector(ends) + ");")
        if (self.stop_model() == "gene"):
            h = float(self.start_params.get("null_dominance", 0.5))
            s_hom = float(self.start_params.get("null_cost", 1.0))
            out += ("\n\tdefineConstant(\"ALL_FUNCTIONAL\", rep(1.0, N_GENES));" +
                    "\n\tdefineConstant(\"NULL_FACTORS\", " + eidos_vector([1.0 - s_hom, 1.0 - h * s_hom, 1.0]) + ");" +
                    "\n\tgene_of_codon = rep(-1, " + str(self.start_params["genome_length"]) + ");" +
                    "\n\tnull_if_changed = gene_of_codon;" +
                    "\n\tnull_if_stop = gene_of_codon;" +
                    "\n\tfor (gene in 0:(N_GENES - 1)){" +
                    "\n\t\tgene_of_codon[GENE_STARTS[gene]:GENE_ENDS[gene]] = gene;" +
                    "\n\t\tnull_if_changed[c(GENE_STARTS[gene], GENE_ENDS[gene])] = gene;" +
                    "\n\t\tif (GENE_ENDS[gene] - GENE_STARTS[gene] > 1) null_if_stop[(GENE_STARTS[gene] + 1):(GENE_ENDS[gene] - 1)] = gene;" +
                    "\n\t}" +
                    "\n\tdefineConstant(\"GENE_OF_CODON\", gene_of_codon);" +
                    "\n\tdefineConstant(\"NULL_IF_CHANGED\", null_if_changed);" +
                    "\n\tdefineConstant(\"NULL_IF_STOP\", null_if_stop);")
        out += ("\n\tprofile_of_codon = sim.getValue(\"fitness_profiles\");" +
                ("\n\tprofile_of_codon[GENE_STARTS] = N_PROF - 2;\n\tprofile_of_codon[GENE_ENDS] = N_PROF - 1;"
                 if self.stop_model() == "site" else "") +
                #Integer lookups for get_genome_fitness: codon (0-63) -> index into AA_LETTERS, and ancestral state per codon
                "\n\tdefineConstant(\"PROFILE_OF_CODON\", profile_of_codon);" +
                "\n\tdefineConstant(\"CODON_TO_AA\", match(codonsToAminoAcids(0:63, paste = F), AA_LETTERS));" +
                "\n\tdefineConstant(\"X_INDEX\", match(\"X\", AA_LETTERS));" +
                "\n\tdefineConstant(\"ANC_CODONS\", sim.chromosome.ancestralNucleotides(format = \"codon\"));" +
                "\n\tdefineConstant(\"ANC_AA\", CODON_TO_AA[ANC_CODONS]);" +
                "\n\tdefineConstant(\"ANC_FIT\", FIT_TABLE[ANC_AA * N_PROF + PROFILE_OF_CODON]);\n")
        return out


    #Non-WF models raise each individual's selection factor to the power NONWF_SELECTION_SCALE (gene-model null-allele
    #factors are not scaled). The value, K / Ne per tick, depends on the life cycle written by write_reproduction and the early()
    #events: utils/nonWFLifeCycle.py derives it and states the assumptions (11/16 with one offspring per tick)

    #Non-WF models compute fitness in each subpopulation's early() event instead of a fitnessEffect() callback, so that
    #survival can be normalized by total fitness: an individual with fitness w survives with probability w * K / sum(w),
    #which keeps the expected population size at K whatever the mean fitness (soft selection)
    def nonwf_fitness_function(self):
        return ("function (void) set_fitness_nonwf (object<Subpopulation>$ subpop, numeric$ K){" +
                "\n\tinds = subpop.individuals;" +
                "\n\tif (size(inds) == 0) return;" +
                ("\n\tw = sapply(inds, \"w1 = get_genome_fitness(applyValue.haploidGenome1); w2 = get_genome_fitness(applyValue.haploidGenome2); " +
                 "(w1[0] * w2[0]) ^ NONWF_SELECTION_SCALE * product(NULL_FACTORS[asInteger(w1[1:N_GENES] + w2[1:N_GENES])]);\");"
                 if self.stop_model() == "gene" else
                 "\n\tw = sapply(inds, \"(get_genome_fitness(applyValue.haploidGenome1) * get_genome_fitness(applyValue.haploidGenome2)) " +
                 "^ NONWF_SELECTION_SCALE;\");") +
                "\n\tinds.fitnessScaling = w;" +
                "\n\ttotal = sum(w);" +
                "\n\tif (total > 0) subpop.fitnessScaling = K / total; else subpop.fitnessScaling = 0.0;\n}\n\n\n")

    #Eidos for the density regulation in a non-WF subpopulation's early() event
    def nonwf_early_commands(self, pop_name, population_size):
        if (self.start_params["neutral_evolution"]):
            return pop_name + ".fitnessScaling = " + str(int(population_size)) + "/" + pop_name + ".individualCount;"
        return "set_fitness_nonwf(" + pop_name + ", " + str(int(population_size)) + ");"


    #Write the reproduction callback for non-Wright-Fisher models: each individual is the first parent of
    #nonWFLifeCycle.OFFSPRING_PER_TICK offspring, each with a random second parent. The non-WF selection scale is derived
    #for this life cycle (utils/nonWFLifeCycle.py), so change both together.
    def write_reproduction(self):
        b = nonWFLifeCycle.OFFSPRING_PER_TICK
        cross = "subpop.addCrossed(individual, subpop.sampleIndividuals(1));"
        reproduction_string = ("reproduction(){\n\t" + (cross if b == 1 else "for (i in seqLen(" + str(b) + ")) " + cross) + "\n}\n\n\n")

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
                        ".haplosomes.nucleotides(NULL, NULL, \"integer\");"+
                        "\n\t\tmuts_mat = matrix(muts_mat, nrow = row_num, byrow = T);" +
                        "\n\t\tcompare_seq = c(muts_mat[0,]);"+
                        "\n\n\t\tfixed_nucs = apply(muts_mat, 1, " +
                        "\"max(applyValue) - min (applyValue) == 0;\");" +
                        "\n\n\t\tdifferent_muts = (ancestral_genome != compare_seq);" +
                        "\n\t\tnew_fixations = different_muts & fixed_nucs;" +
                        "\n\n\t\tif(any(new_fixations)){" +
                        "\n\t\t\tnew_fixed = ancestral_genome;"+
                        "\n\t\t\tnew_fixed[new_fixations] = compare_seq[new_fixations];")
                            
                            
            #If there is a flag to also calculate selection - ie. count fixed dN/dS, figure out if dN or dS
            if(population_parameters["calculate_selection"]):
                repeated_commands_string += ("\n\t\t\tnew_fixations_space = which(new_fixations);" +
                                "\n\n\t\t\tdN_name = \"dN_" + pop_name + "\";" +
                                "\n\t\t\tdS_name = \"dS_" + pop_name + "\";" +
                                "\n\t\t\tsubs_name = \"subs_" + pop_name + "\";"+
                                "\n\t\t\tfor(fix in new_fixations_space){" +
                                "\n\t\t\t\tfix_pos = (fix + 1) % 3;" +
                                "\n\t\t\t\tif (fix_pos == 0) {" +
                                "\n\t\t\t\t\told_nucs = ancestral_genome[(fix-2):fix];" +
                                "\n\t\t\t\t\tnew_nucs = new_fixed[(fix-2):fix];" +
                                "\n\t\t\t\t} else if (fix_pos == 1) {" +
                                "\n\t\t\t\t\told_nucs = ancestral_genome[fix:(fix+2)];" +
                                "\n\t\t\t\t\tnew_nucs = new_fixed[fix:(fix+2)];" +
                                "\n\t\t\t\t} else if (fix_pos == 2) {" +
                                "\n\t\t\t\t\told_nucs = ancestral_genome[(fix-1):(fix+1)];" +
                                "\n\t\t\t\t\tnew_nucs = new_fixed[(fix-1):(fix+1)];" +
                                "\n\t\t\t\t}"+
                                "\n\n\t\t\t\told_AA = codonsToAminoAcids(nucleotidesToCodons(old_nucs));" +
                                "\n\t\t\t\tnew_AA = codonsToAminoAcids(nucleotidesToCodons(new_nucs));" +
                                "\n\t\t\t\tif (old_AA !=  \"X\" & new_AA != \"X\"){" +
                                "\n\t\t\t\t\tif (old_AA == new_AA & new_AA != \"X\"){" +
                                "\n\t\t\t\t\t\tsim.setValue(dS_name, sim.getValue(dS_name) + 1);" + 
                                "\n\t\t\t\t\t} else {" +
                                "\n\t\t\t\t\t\tsim.setValue(dN_name, sim.getValue(dN_name) + 1);" +
                                "\n\t\t\t\t\t}" +
                                "\n\t\t\t\t\tsim.setValue(subs_name, sim.getValue(subs_name)+" +
                                "paste0(\"\\n\",asInteger(floor(fix/3)), old_AA, \":\", new_AA));" +
                                "\n\t\t\t\t}\n\t\t\t\t}")
            
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
                    "\n\tsim.setValue(\"dS_"+ population_parameters["pop_name"]+"\", 0);" +
                    "\n\tsim.setValue(\"subs_"+ population_parameters["pop_name"]+"\", \"\\n\\nSubstitutions:\");")

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
                                    "\n\tsim.setValue(\"fixations_counted_"+ pop_name+"\", 0);" +
                                    "\n\tsim.setValue(\"dN_"+ pop_name+"\", 0);" +
                                    "\n\tsim.setValue(\"dS_"+ pop_name+"\", 0);" +
                                    "\n\tsim.setValue(\"subs_"+ pop_name+"\", \"\\n\\nSubstitutions:\");")


            define_population_string += "\n}\n\n\n"

            self.output_file.write(define_population_string)

        #Write the early commands - this may need tweaking w/ the fitness algorithm
        early_event = (str(int(population_parameters["dist_from_start"]) + 1) + ":" + str(int(population_parameters["end_dist"])) +
                        " early(){\n\t" + self.nonwf_early_commands(pop_name, population_parameters["population_size"]))

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
            end_population_string += ("\n\twriteFile(\"" + self.start_params["filenames"][5] +"/"+ population_parameters["clade_name"] + 
                "_num_subs.txt\"," + "asString(sim.getValue(\"fixations_counted_" + population_parameters["pop_name"] + "\")));")

        #Write file with the number of synonymous and synonymous substitutions
        if(population_parameters["calculate_selection"]):
            end_population_string += ("\n\twriteFile(\"" + self.start_params["filenames"][4] +"/"+ population_parameters["clade_name"] + "_dNdS.txt\"," +
                "paste0(\"dN: \", sim.getValue(\"dN_" + population_parameters["pop_name"] + "\"), \" / " + str(self.start_params["dn_denom"]) + " = \", " +
                "sim.getValue(\"dN_" + population_parameters["pop_name"] + "\") / " + str(self.start_params["dn_denom"]) + ", " +
                "\"\\ndS: \", sim.getValue(\"dS_" + population_parameters["pop_name"] + "\"), \" / " + str(self.start_params["ds_denom"]) + " = \", " +
                "sim.getValue(\"dS_" + population_parameters["pop_name"] + "\") / " + str(self.start_params["ds_denom"]) + ", " +
                "sim.getValue(\"subs_" + population_parameters["pop_name"] + "\")));" )

        #Write files containing polymorphisms in each population and relative proportions
        if(population_parameters["polymorphisms"]):
            end_population_string += ("\n\tpop_seq = sample("+ population_parameters["pop_name"] +".individuals.haplosomes, 1).nucleotides();" +
                                        "\n\tpop_seq = strsplit(codonsToAminoAcids(nucleotidesToCodons(pop_seq)), sep = \"\");" +
                                        "\n\tpolymorph_str = c();\n\tfixed_str=c();" +
                                        "\n\tfor (a in 0:(length(pop_seq)-1)) " +
                                        "{\n\t\tdiffs = c();" + 
                                        "\n\t\tfor (g in " + population_parameters["pop_name"] + ".individuals.haplosomes.nucleotides()){" +
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
                                        "\n\t}\n\twriteFile(\"" + self.start_params["filenames"][6]  + "/" + population_parameters["clade_name"] + "_polymorphisms.txt\", " +
                                        "paste(polymorph_str, sep = \"\"));" +
                                        "\n\twriteFile(\"" + self.start_params["filenames"][6] +"/" + population_parameters["clade_name"] + 
                                        "_fixed_sites.txt\", paste(fixed_str, sep = \"\"));")


        if(population_parameters["terminal_clade"] and self.start_params["nonWF"]):
            end_population_string += "\n\t" + population_parameters["pop_name"] + ".removeSubpopulation();"

        end_population_string += "\n}\n\n\n"

        self.output_file.write(end_population_string)






    #Write code to write the output for terminal populations after they have reached their population
    def write_terminal_output(self, population_parameters, pop = "p1"):


        #Set up sampling of the population
        clade_name = population_parameters["clade_name"]
        pop_size = population_parameters["population_size"]
        samp_size = population_parameters["sample_size"]



        #Formulate and output consensus sequence
        if (samp_size == "consensus"):
            terminal_output_string = ("\n\n\tconsensus = \"\";" +
                                    "\n\tfor (i in 0:" + str(self.start_params["genome_length"] * 3 - 1) + "){" +
                                    "\n\t\tconsensus = consensus+ c(\"A\", \"C\", \"G\", \"T\")[whichMax(nucleotideCounts(paste0(matrix(sapply(" + pop +
                                    ".haplosomes.nucleotides(), \"strsplit(applyValue, sep = '');\"), ncol = " + str(self.start_params["genome_length"] * 3) + 
                                    ", byrow = T)[,i])))];\n\t}" +
                                    "\n\n\tfasta_string_nuc = paste0(\">" + clade_name + ": \\n\", consensus);" +
                                    "\n\twriteFile(\"" + self.start_params["filenames"][7] + "/nuc.fasta" + "\", fasta_string_nuc,append = T);" +
                                    "\n\n\tfasta_string_prot = paste0(\">" + clade_name + ": \\n\", codonsToAminoAcids(nucleotidesToCodons(consensus)));" +
                                        "\n\twriteFile(\"" + self.start_params["filenames"][8] + "/AA.fasta" + "\", fasta_string_prot,append = T);")

        else:
            if(samp_size == "all"):
                terminal_output_string = "\n\tgenomes = " + pop + ".haplosomes;"
            else:
                terminal_output_string = ("\n\tgenomes = sample(" + pop + ".haplosomes, min(" + str(int(samp_size)) +
                                    ", 2*" + pop + ".individualCount), replace=F);")


            
            #Iterate through each random sample to write script to output samples of amino acids and nucleotides to fasta files
            terminal_output_string += ("\n\n\tfor (g in genomes){" +
                                        "\n\t\tfasta_string_nuc = \"\";" +
                                        "\n\t\tfasta_string_nuc = paste0(\">\", g.individual, \", " + clade_name + ": \\n\", g.nucleotides());" +
                                        "\n\t\twriteFile(\"" + self.start_params["filenames"][7] + "/" + clade_name + "_nuc.fasta" + "\", fasta_string_nuc,append = T);" +
                                        "\n\t\tfasta_string_prot = paste0(\">\", g.individual, \", " + clade_name +
                                        ": \\n\", codonsToAminoAcids(nucleotidesToCodons(g.nucleotides())));" +
                                        "\n\t\twriteFile(\"" + self.start_params["filenames"][8] + "/" + clade_name + "_aa.fasta" + "\", fasta_string_prot,append = T);}" )            


        return terminal_output_string





    #Closes the file that has been appended to
    def close_file(self):
        self.output_file.close()
