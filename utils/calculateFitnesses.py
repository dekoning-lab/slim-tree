import numpy as np
import scipy.optimize as sp
import csv

#parallel library
import multiprocessing

#USED IN TESTING!
import time


def KLD(truth, est):
    '''
    :param truth: Array of Float, true stationary distribution of codons.
    :param est: Array of Float, estimated stationary distribution of codons.
    :return: Kullback Leibler Divergence between the two distributions.
    '''
    kld = 0.0
    for i in range(len(truth)):
        if(est[i] == 0.0 and est[i] != truth[i]):
            kld += (truth[i] * np.log2(truth[i]/np.finfo(float).eps))
        elif(est[i] == truth[i]):
            continue
        else:
            kld += (truth[i] * np.log2(truth[i]/est[i]))
    return(kld)

def Get_nucleotide_from_idx(codon_idx):
    '''
    :param codon_idx: Int, the index of the codon that we want to get the nucleotides of.
    :return: String, the nucleotides associated with the codon idx.
    '''
    # codon map
    codonMap = {11:"TTT", 10:"TTC", 27:"TTA", 28:"TTG", 47:"TCT", 45:"TCC",
                44:"TCA", 46:"TCG", 59:"TAT", 60:"TAC", 5:"TGT",  4:"TGC",
                58:"TGG", 26:"CTT", 24:"CTC", 23:"CTA", 25:"CTG", 35:"CCT",
                33:"CCC", 32:"CCA", 34:"CCG", 17:"CAT", 16:"CAC", 36:"CAA",
                37:"CAG", 43:"CGT", 41:"CGC", 40:"CGA", 42:"CGG", 20:"ATT",
                19:"ATC", 18:"ATA", 29:"ATG", 53:"ACT", 51:"ACC", 50:"ACA",
                52:"ACG", 31:"AAT", 30:"AAC", 21:"AAA", 22:"AAG", 48:"AGT",
                49:"AGC", 38:"AGA", 39:"AGG", 57:"GTT", 55:"GTC", 54:"GTA",
                56:"GTG", 3:"GCT",   1:"GCC",  0:"GCA",  2:"GCG",  7:"GAT",
                6:"GAC",  8:"GAA",  9:"GAG",  15:"GGT", 13:"GGC", 12:"GGA",
                14:"GGG"}
    return(codonMap[codon_idx])

def Get_idx_from_nucleotide(codon_nucleotieds):
    '''
    :param codon_idx: String, Coded triplet for a codon.
    :return: Int, the codon index associated with nucleotides.
    '''
    # codon map
    codonMap = {"TTT":11, "TTC":10, "TTA":27, "TTG":28, "TCT":47, "TCC":45,
                "TCA":44, "TCG":46, "TAT":59, "TAC":60, "TGT":5 , "TGC":4,
                "TGG":58, "CTT":26, "CTC":24, "CTA":23, "CTG":25, "CCT":35,
                "CCC":33, "CCA":32, "CCG":34, "CAT":17, "CAC":16, "CAA":36,
                "CAG":37, "CGT":43, "CGC":41, "CGA":40, "CGG":42, "ATT":20,
                "ATC":19, "ATA":18, "ATG":29, "ACT":53, "ACC":51, "ACA":50,
                "ACG":52, "AAT":31, "AAC":30, "AAA":21, "AAG":22, "AGT":48,
                "AGC":49, "AGA":38, "AGG":39, "GTT":57, "GTC":55, "GTA":54,
                "GTG":56, "GCT":3 , "GCC":1 , "GCA":0 , "GCG":2 , "GAT":7 ,
                "GAC":6 , "GAA":8 , "GAG":9 , "GGT":15, "GGC":13, "GGA":12,
                "GGG":14}
    return(codonMap[codon_nucleotieds])

def Get_amino_acid(codon_idx):
    '''
    :param codon_idx: Int, Codon idx used in MutSel model
    :return: Int,
    '''
    aminoAcidMap = {11:4,  10:4,  27:9,  28:9,  47:15, 45:15,
                    44:15, 46:15, 59:19, 60:19, 5:1,   4:1,
                    58:18, 26:9,  24:9,  23:9,  25:9,  35:12,
                    33:12, 32:12, 34:12, 17:6,  16:6,  36:13,
                    37:13, 43:14, 41:14, 40:14, 42:14, 20:7,
                    19:7,  18:7,  29:10, 53:16, 51:16, 50:16,
                    52:16, 31:11, 30:11, 21:8,  22:8,  48:15,
                    49:15, 38:14, 39:14, 57:17, 55:17, 54:17,
                    56:17, 3:0,   1:0,   0:0,   2:0,   7:2,
                    6:2,   8:3,   9:3,   15:5,  13:5,  12:5,
                    14:5}
    return(aminoAcidMap[codon_idx])

def Get_codons(amino_acid_idx):
    '''
    :param codon_idx: Int, Codon idx used in MutSel model
    :return: Int,
    '''
    codonMap = {0:[0,1,2,3],
                1:[4,5],
                2:[6,7],
                3:[8,9],
                4:[10,11],
                5:[12,13,14,15],
                6:[16,17],
                7:[18,19,20],
                8:[21,22],
                9:[23,24,25,26,27,28],
                10:[29],
                11:[30,31],
                12:[32,33,34,35],
                13:[36,37],
                14:[38,39,40,41,42,43],
                15:[44,45,46,47,48,49],
                16:[50,51,52,53],
                17:[54,55,56,57],
                18:[58],
                19:[59,60]}

    return(codonMap[amino_acid_idx])

def Read_frequency_file(file_name):
    '''
    :param file_name: String, the name of the csv file containing codon frequencies.
    :return: List of (Lists of floats), a list of the codon frequencies in the correct order.
    '''
    with open(file_name) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        freqs = []
        for row in spamreader:
            #Initalize data structure
            if(len(freqs) == 0):
                for i in range(len(row) - 1):
                    freqs.append([0]*61)

            #reorder if needed
            codon_idx = Get_idx_from_nucleotide(row[0])
            for i in range(1, len(row)):
                freqs[i-1][codon_idx] = float(row[i])
    return(freqs)
def Calc_stationary_KLD(A):
    '''
    :param A: numpy array, an instantaneous rates matrix.
    :return: List of Floats, the stationary distribution of the instantaneous rates matrix.
    '''
    size = len(A[0,])
    A_prime = np.vstack([A.transpose(), np.ones(size)])
    b = np.zeros(size+1) #target solution
    b[size] = 1.0        #make sure freqs sum to 1
    pi = np.linalg.lstsq(A_prime, b)[0]
    #Square the matrix using QR decomp
    #Q ,R = sp.linalg.qr(A_prime)
    #pi = np.linalg.lstsq(R ,np.dot(Q.transpose(), b)) #(Numpy only takes square martices DUMB)
    for idx in range(len(pi)):
        if(pi[idx] < 0.0 and -1e-12 < pi[idx]):
            pi[idx] = np.finfo(float).eps
        elif(pi[idx] < -1e-5):
            raise Exception("Large negative frequency detected" + str(pi[idx]) +"\n" + ",".join(str(y) for y in pi))
    return(pi/sum(pi))

def Compute_matrix(fitness, pop, mu = 1.0):
    '''
    NO MUTATION RATE NEEDED FOR THIS APPLICATION!!!
    :param fitness: List of Float, the expected number of offspring for an individual with a particular amino acid.
                                    (Amino Acids in alphabetic order)
    :param pop: Float, the population size (0 <= pop).
    :return: Numpy array, instantaneous rates matrix for a mutation selection model.
    '''
    rates_matrix = np.zeros((61,61))
    for i in range(61):
        for j in range(61):
            #get Hamming distance
            nucs_i = Get_nucleotide_from_idx(i)
            nucs_j = Get_nucleotide_from_idx(j)
            hamm = 0
            for bp  in range(3):
                if(nucs_i[bp] != nucs_j[bp]):
                    hamm += 1

            fit_i = fitness[Get_amino_acid(i)]
            fit_j = fitness[Get_amino_acid(j)]

            if(hamm == 1):
                #syn
                if(abs(fit_j - fit_i) <= 1e-30):
                    rates_matrix[i,j] = mu
                else:
                    rates_matrix[i,j] = mu * 2 * pop * ((1 - np.exp(-2 * (fit_j - fit_i))) / (1 - np.exp(-4 *pop * (fit_j - fit_i))))
        rates_matrix[i,i] = -1 * rates_matrix[i,].sum()
    return(rates_matrix)




def Compute_matrix_w_HKY(fitness, pop, nuc_freqs, kappa,  mu = 1.0):
    '''
    NO MUTATION RATE NEEDED FOR THIS APPLICATION!!!
    :param fitness: List of Float, the expected number of offspring for an individual with a particular amino acid.
                                    (Amino Acids in alphabetic order)
    :param pop: Float, the population size (0 <= pop).
    :return: Numpy array, instantaneous rates matrix for a mutation selection model.
    '''

    nuc_to_idx = {'A':0, 'C':1, 'G':2, 'T':3}
    idx_to_nuc = {0:'A', 1:'C', 2:'G', 3:'T'}
    transversion_dic = {(0,2):False, (2,0):False, (1,3):False, (3,1):False,
                        (0,1): True, (1,0): True, (0,3): True, (3,0): True,
                        (2,1): True, (1,2): True, (2,3): True, (3,2): True}
    hky_matrix = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if(i != j):
                if(transversion_dic[(i,j)]):
                    hky_matrix[i,j] = nuc_freqs[j] * kappa * mu
                else:
                    hky_matrix[i,j] = nuc_freqs[j] * mu
        hky_matrix[i,i] = -1 * hky_matrix[i,].sum()

    rates_matrix = np.zeros((61,61))
    for i in range(61):
        for j in range(61):
            #get Hamming distance
            nucs_i = Get_nucleotide_from_idx(i)
            nucs_j = Get_nucleotide_from_idx(j)
            hamm = 0
            dif_idx = -1
            for bp  in range(3):
                if(nucs_i[bp] != nucs_j[bp]):
                    hamm += 1
                    dif_idx = bp

            fit_i = fitness[Get_amino_acid(i)]
            fit_j = fitness[Get_amino_acid(j)]

            if(hamm == 1):
                x = nuc_to_idx[nucs_i[dif_idx]]
                y = nuc_to_idx[nucs_j[dif_idx]]

                #syn
                if(abs(fit_j - fit_i) <= 1e-30):
                    rates_matrix[i,j] = hky_matrix[x,y]
                else:
                    rates_matrix[i,j] = hky_matrix[x,y] * 2 * pop * ((1 - np.exp(-2 * (fit_j - fit_i))) / (1 - np.exp(-4 *pop * (fit_j - fit_i))))
        rates_matrix[i,i] = -1 * rates_matrix[i,].sum()
    return(rates_matrix)



def opt_ratio(selection, pop, ratio):
    '''
    This is the objective function of the optimization algorithm.
    :param selection: selection coefficient between amino acids
    :param pop: The number of individuals in the population
    :param ratio: Ratio pre-computed in the prior function
    :return:
    '''
    numerator = ((1 - np.exp(2 * selection)) / (1 - np.exp(4 * pop * selection)))
    denominator = ((1 - np.exp(-2 * selection)) / (1 - np.exp(-4 * pop * selection)))

    if(denominator <= 0):
        denominator = np.finfo(float).eps
    if(numerator <= 0):
        numerator = np.finfo(float).eps
    return(abs((numerator/denominator) - ratio))

def ratio_func(pi_i, pi_j, c_aa_i, c_aa_j, pop,  mu_ij = 0.0, mu_ji = 0.0):
    '''
    :param pi_i: proportion stationary frequency of codon i wrt to amino acid AA(i) under a neutral model(0 < pi_i < 1.0)
    :param pi_j:  proportion stationary frequency of codon j wrt to amino acid AA(j) under a neutral model (0 < pi_j < 1.0)
    :param c_aa_i: correction factor for frequency of codon i based on amino acid frequency that i maps to
    :param c_aa_j: correction factor for frequency of codon j based on amino acid frequency that j maps to
    :param pop: the population size (0 <= pop)
    :param mu_ij: mutation rate from i to j under nucleotide mode (default assumed mu_ij = mu_ji)
    :param mu_ji: mutation rate from j to i under nucleotide mode (default assumed mu_ij = mu_ji)
    :return: selection coefficient between i and j, s_ij
    '''
    ratio_const = 0.0
    if(mu_ij == 0.0 or mu_ji == 0.0):
        ratio_const = (pi_i * c_aa_i) / (pi_j * c_aa_j)
    else:
        ratio_const = (pi_i * c_aa_i * mu_ij) / (pi_j * c_aa_j * mu_ji)

    #catch cases that are neutral
    if(abs(ratio_const - 1.0 ) < 1e-5):
        return(0.0)

    if(ratio_const < 1.0):
        invert = True
        ratio_const = 1.0 / ratio_const
    else:
        invert = False

    #Solve for the (negative) selection coefficient in closed form. With Kimura's fixation probability
    #u(s) = (1 - exp(-2s)) / (1 - exp(-4 pop s)), the ratio u(-s) / u(s) that opt_ratio matches equals
    #exp(-(4 pop - 2) s) exactly, so s = -log(ratio) / (4 pop - 2). The numerical optimizer this replaces
    #(L-BFGS-B from -1e-7 with tol 1e-7) lost accuracy as pop grew: errors in 4 pop s reached ~0.1 at pop = 1000
    #and ~10 at pop = 10000. The bounds of the old search, [-1, -1e-7], are kept as a clamp.
    sel = -np.log(ratio_const) / (4.0 * pop - 2.0)
    sel = min(-1e-7, max(-1.0, sel))

    if(invert):
        sel = sel * -1.0
    return(sel)

def Find_nuc_freqs(true_distribution, pop):
    """
    :param true_distribution: Codon distribution we want to compute the HKY parameters for. Kappa is NON-IDENTIFYABLE!
    :return: Neutral stationary frequencies, the mutation matrix
    """
    nuc_dic = {'A':0, 'C':1, 'G':2, 'T':3}
    #Check if amino acids have non uniform frequecies for codon rep
    same_freq = True
    aa_freqs = [[] for i in range(20)]

    for i in range(61):
        aa_freqs[Get_amino_acid(i)].append(true_distribution[i])
    for aa in aa_freqs:
        for codon in aa:
            if(abs(codon - (sum(aa)/len(aa))) > 1e-5):
                same_freq = False

    if(same_freq):
        return(([1.0/61.0] * 61, np.ones((4,4)))) #NO BIAS JK69
    else:
        #find syn codons to find mutation/representation bias
        ratios = []
        nuc_tally = [0,0,0,0]
        for i in range(61):
            for j in range(i,61):
                if(i != j):
                    #get Hamming distance
                    nucs_i = Get_nucleotide_from_idx(i)
                    nucs_j = Get_nucleotide_from_idx(j)
                    hamm = 0
                    dif_idx = -1
                    for bp  in range(3):
                        if(nucs_i[bp] != nucs_j[bp]):
                            hamm += 1
                            dif_idx = bp
                    aa_i = Get_amino_acid(i)
                    aa_j = Get_amino_acid(j)

                    if(hamm == 1) and (aa_i == aa_j):
                        if(abs(true_distribution[i] - true_distribution[j]) > 1e-5):
                            if(true_distribution[i] < true_distribution[j]):
                                ratios.append((true_distribution[j]/ true_distribution[i], nuc_dic[nucs_i[dif_idx]], nuc_dic[nucs_j[dif_idx]]))
                                nuc_tally[nuc_dic[nucs_j[dif_idx]]] += 1
                            elif(true_distribution[i] > true_distribution[j]):
                                ratios.append((true_distribution[i]/ true_distribution[j], nuc_dic[nucs_j[dif_idx]], nuc_dic[nucs_i[dif_idx]]))
                                nuc_tally[nuc_dic[nucs_i[dif_idx]]] += 1

        nuc_freqs = [1.0, 1.0, 1.0, 1.0]
        nucs_infered = [0,0,0,0] #counts the most common nuc
        #fix the most frequent nuc
        nucs_infered[max(range(len(nuc_tally)), key=nuc_tally.__getitem__)] = 1
        for frac in ratios:
            (ratio, i, j) = frac
            if(nucs_infered[j] == 1 ):
                nuc_freqs[i] = nuc_freqs[j] / ratio
                nucs_infered[i] = 1
            if(sum(nucs_infered) == 4):
                break
        nuc_freqs /= sum(nuc_freqs)

    #make HKY85 NO KAPPA :(
    hky_matrix = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if(i != j):
                hky_matrix[i,j] = nuc_freqs[j]
        hky_matrix[i,i] = -1 * hky_matrix[i,].sum()

    stationary = Calc_stationary_KLD(Compute_matrix_w_HKY([1.0]*20, pop, nuc_freqs, 1.0, 1.0))
    return((stationary,hky_matrix))

def Find_fits(true_distribution, pop, eps=1e-5):
    '''
    NUMERICALLY STABLE FOR FREQUENCIES DOWN TO 1e-10
    :param true_distribution: the frequency of each codon
    :param pop: the population size (0 <= pop)
    :param mu: the mutation rate (default we ignore mutational bias)
    :param eps: the epsilon value (default 1e-5)
    :return:
    '''

    #get amino acid frequencies
    codon_pi = true_distribution
    epsilon = 1e-10
    codon_pi = [epsilon if x < epsilon else x for x in codon_pi]
    codon_pi = [x/sum(codon_pi) for x in codon_pi] #ideal stationary distribution!

    aa_pi = np.zeros(20)
    aa_pi_neutral = np.zeros(20)
    aa_fitness_key = np.zeros(20)

    for i in range(61):
        aa_pi[Get_amino_acid(i)] += codon_pi[i]
    aa_pi /= sum(aa_pi)

    #fix the most frequent amino acid
    max_aa_idx = np.argmax(aa_pi)
    aa_fitness_key[max_aa_idx] = 1.0

    #find the neutral stationary distribution
    neutral_freqs, mu_matrix = Find_nuc_freqs(true_distribution, pop)
    nuc_dic = {'A':0, 'C':1, 'G':2, 'T':3}

    #this is used to account for how much each codon contributes to amino acid frequency under a neutral model
    for i in range(61):
        aa_pi_neutral[Get_amino_acid(i)] += neutral_freqs[i]
    aa_pi_neutral /= sum(aa_pi_neutral)

    #Find codons that have different amino acids
    for i in Get_codons(max_aa_idx):
        for j in range(61):
            #get Hamming distance
            nucs_i = Get_nucleotide_from_idx(i)
            nucs_j = Get_nucleotide_from_idx(j)
            hamm = 0
            diff_idx = -1
            for bp  in range(3):
                if(nucs_i[bp] != nucs_j[bp]):
                    hamm += 1
                    diff_idx = bp
            aa_i = Get_amino_acid(i)
            aa_j = Get_amino_acid(j)

            #make sure the codons differ by one BP
            if(hamm == 1):
                if(aa_i == max_aa_idx and aa_j != max_aa_idx):
                    if(aa_fitness_key[aa_j] == 0):

                        c_aa_i =  aa_pi[aa_i] #desired frequency of amino acid associated with codon i
                        c_aa_j =  aa_pi[aa_j] #desired frequency of amino acid associated with codon j
                        proption_i = neutral_freqs[i] / aa_pi_neutral[aa_i] #proportion of amino acid represented by codon i when selection is not a factor
                        proption_j = neutral_freqs[j] / aa_pi_neutral[aa_j] #proportion of amino acid represented by codon j when selection is not a factor
                        mu_ij = mu_matrix[nuc_dic[nucs_i[diff_idx]],nuc_dic[nucs_j[diff_idx]]] #probability of mutation from i to j
                        mu_ji = mu_matrix[nuc_dic[nucs_j[diff_idx]],nuc_dic[nucs_i[diff_idx]]] #probability of mutation from j to i
                        aa_fitness_key[aa_j] = 1.0 + ratio_func(proption_i,proption_j, c_aa_i, c_aa_j, pop, mu_ij, mu_ji)

    #find missing fitness codons
    missing_fits = []
    for aa_idx in range(20):
        if(aa_fitness_key[aa_idx] == 0):
            missing_fits.extend(Get_codons(aa_idx))

    #Keep going until all fits are found
    while len(missing_fits) > 0:
        for i in missing_fits:
            aa_i = Get_amino_acid(i)
            if(aa_fitness_key[aa_i] != 0.0):
                continue
            for j in range(61):
                aa_j = Get_amino_acid(j)
                if(aa_fitness_key[aa_j] == 0.0):
                    continue

                #get Hamming distance
                nucs_i = Get_nucleotide_from_idx(i)
                nucs_j = Get_nucleotide_from_idx(j)
                hamm = 0
                diff_idx = -1
                for bp  in range(3):
                    if(nucs_i[bp] != nucs_j[bp]):
                        hamm += 1
                        diff_idx = bp

                if(hamm == 1):
                    c_aa_i =  aa_pi[aa_i] #desired frequency of amino acid associated with codon i
                    c_aa_j =  aa_pi[aa_j] #desired frequency of amino acid associated with codon j
                    proption_i = neutral_freqs[i] / aa_pi_neutral[aa_i] #proportion of amino acid represented by codon i when selection is not a factor
                    proption_j = neutral_freqs[j] / aa_pi_neutral[aa_j] #proportion of amino acid represented by codon j when selection is not a factor
                    mu_ij = mu_matrix[nuc_dic[nucs_i[diff_idx]],nuc_dic[nucs_j[diff_idx]]] #probability of mutation from i to j
                    mu_ji = mu_matrix[nuc_dic[nucs_j[diff_idx]],nuc_dic[nucs_i[diff_idx]]] #probability of mutation from j to i
                    #run it backwards since we know the fitness of j not i
                    aa_fitness_key[aa_i] = aa_fitness_key[aa_j] + ratio_func(proption_j,proption_i, c_aa_j, c_aa_i, pop, mu_ji, mu_ij)
                    #remove the codons for the fitness we found
                    for x in Get_codons(aa_i):
                        missing_fits.remove(x)  #allow for stopping condition
                    break
    return(aa_fitness_key)
def Find_fits_parallel_boilerplate(params):
    '''
    :param params: List of Tuple, (Codon distribution, population size)
    :return: list of fitness profiles.
    '''
    return(Find_fits(params[0], params[1]))

def Batch_find_fits(true_distributions, pops, thread_count):
    '''
    :param true_distributions: List of (List of Float), A list of the TRUE stationary distribution of codons.
    :param pops: List of Float, A list of the Population sizes (0 <= pops).
    :param thread_count: Int, the number of threads to use 0 < thread_count.
    :return: List of (List of Float), fitness values from optimized solutions.
    '''
    # multiprocessing pool object
    pool = multiprocessing.Pool()

    # pool object with number of element
    pool = multiprocessing.Pool(processes=thread_count)

    # input list
    inputs = []
    for i in range(len(pops)):
        inputs.append((true_distributions[i], pops[i]))

    # function and input list as arguments
    outputs = pool.map(Find_fits_parallel_boilerplate, inputs)
    return(outputs)

def Find_fitnesses(filename, population_size, is_parallel, thread_count=10):
    '''
    :param filename: String, the name of the csv file to be read with stationary distributions.
    :param population_size: Float, the population size (0 <= population_size).
    :param is_parallel: Bool, True=HPC mode, False=Serial mode.
    :param thread_count: Int, the number of threads to use for HPC MODE 0 < thread_count.
    :return: List of fitness profiles found.
    '''
    freqs = Read_frequency_file(filename)
    pops = [population_size] * len(freqs)
    if(is_parallel):
        if(len(freqs) < thread_count):
            thread_count = len(freqs)
        fits = Batch_find_fits(freqs, pops, thread_count)
    else:
        fits = []
        for idx in range(len(freqs)):
            fits.append(Find_fits(freqs[idx], pops[idx]))
    return(fits)












#MULTI THREAD GAURD! HARDCODED SECURITY IN MULTITHRED LIB TO STOP FORK BOMB ACCIDENTS!
if __name__ == '__main__':


    #TEST CASES BELOW
    case1 = False   #Basic
    case2 = False   #Basic
    case3 = False   #Basic
    case4 = False   #Basic
    case5 = False   #Contrived Basic
    case6 = False   #Contrived Basic
    case7 = False   #Contrived Basic
    case8 = False   #Contrived Basic
    case9 = False   #Real Case
    case10 = False  #Real Case
    case11 = False  #HKY TEST CASE
    case12 = False  #Serial Time Case
    case13 = True  #Parallel Time Case
    case14 = False  #Real Challenging Case


    #TEST CASE1
    if(case1):
        true_fits = [1.0, 1.0, 1.0, .99, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0]
        Pop = 100

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))
        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #BASIC TEST
    if(case2):
        true_fits = [.99, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0]
        Pop = 100

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #CORNER CASE1
    if(case3):
        true_fits = [1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, .99]
        Pop = 100

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #CORNER CASE2
    if(case4):
        true_fits = [1.0, .99, .99, .99, .99,
                     .99, .99, .99, .99, .99,
                     .99, .99, .99, .99, .99,
                     .99, .99, .99, .99, .99]
        Pop = 100

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #CORNER CASE3
    if(case5):
        true_fits = [1.0, .99, .99, .99, .99,
                     .99, .99, .99, .99, .99,
                     .99, .99, .99, .99, .99,
                     .99, .99, .99, .99, .99]
        Pop = 10

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #VARITY CASE
    if(case6):
        true_fits = [1.0, .99, .98, .97, .96,
                     .95, .94, .93, .92, .91,
                     .90, .99, .98, .97, .96,
                     .95, .94, .93, .92, .91]
        Pop = 10

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #MAGNITUTE CASE
    if(case7):
        true_fits = [1.0, .9, 1.0, 1.0, 1.0,
                     1.0, .9, 1.0, 1.0, 1.0,
                     1.0, .8, 1.0, 1.0, 1.0,
                     1.0, .8, 1.0, 1.0, 1.0]
        Pop = 5

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #MAGNITUTE CASE2
    if(case8):
        true_fits = [1.0, .99, 1.0, 1.0, 1.0,
                     1.0, .99, 1.0, 1.0, 1.0,
                     1.0, .98, 1.0, 1.0, 1.0,
                     1.0, .98, 1.0, 1.0, 1.0]
        Pop = 100

        freqs = Calc_stationary_KLD(Compute_matrix(true_fits, Pop))

        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(true_fits, "\n",est_fits)
        print("KLD: ", KLD(freqs, est_freqs))

    #REAL CASE
    if(case9):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        Pop = 100
        freqs = freqs[0]
        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(freqs, "\n",est_freqs)
        print("KLD: ", KLD(freqs, est_freqs))

    #REAL CASE2
    if(case10):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        Pop = 100

        freqs = freqs[1]
        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(freqs, "\n",est_freqs)
        print("KLD: ", KLD(freqs, est_freqs))

    #HKY TEST
    if(case11):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        Pop = 100
        freqs = freqs[2]
        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))

        est_freqs2 = Calc_stationary_KLD(Compute_matrix_w_HKY(est_fits, Pop, [0.1,0.5,0.15,0.25], 10.0))
        est_fits2 = Find_fits(est_freqs2, Pop)

        print(est_fits,"\n", est_fits2 )


    #Sequential TIME Test
    if(case12):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        pops = [100]*50
        Total_time = 0.0
        for idx in range(len(freqs)):
            start = time.time()
            est_fits = Find_fits(freqs[idx], pops[idx])
            end = time.time()
            est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, pops[idx]))
            print("Frequency_distribution_",str(idx),"KLD:\t {res:.3E}".format(res=KLD(freqs[idx], est_freqs)), "\tTime: {res:.1f}".format(res=end - start),"Sec")
            Total_time += (end - start)
        print("TOTAL TIME: ", Total_time)

    #Parallel TIME Test
    if(case13):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        pops = [100]*50
        start = time.time()
        est_fits = Batch_find_fits(freqs, pops, 10)
        end = time.time()
        print("TOTAL TIME: ", (end - start))

    #Challenging real case
    if(case14):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        pop = 100
        freqs = freqs[24]
        est_fits = Find_fits(freqs, pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, pop))
        print("KLD: ", KLD(freqs, est_freqs))