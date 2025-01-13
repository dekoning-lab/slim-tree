import numpy as np
import scipy as sp
import csv
import warnings

#parallel library
import multiprocessing

#USED IN TESTING!
import time

warnings.filterwarnings("ignore")

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

def Compute_matrix(fitness, pop):
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
                    rates_matrix[i,j] = 1.0
                else:
                    rates_matrix[i,j] = 1.0 * 2 * pop * ((1 - np.exp(-2 * (fit_j - fit_i))) / (1 - np.exp(-4 *pop * (fit_j - fit_i))))
        rates_matrix[i,i] = -1 * rates_matrix[i,].sum()
    return(rates_matrix)


def Frequency_distance(fits, true_freqs, pop, fixed_idx):
    '''
    This is the function called by the optimizer
    :param fits: List of Float, the expected number of offspring for an individual with a particular amino acid.
                                (Amino Acids in alphabetic order, 0.0 <= fitness <= 2.0 ) ***ONE MISSING***
    :param true_freqs: List of Float, the TRUE stationary distribution of codons
    :param pop: Float, the population size (0 <= pop).
    :param fixed_idx: Int, value of the index of the amino acid with a fixed fitness value.
    :return: Float, KLD between the true stationary distribution and the MutSel stationary distribution.
    '''
    fits = np.insert(fits, fixed_idx, 1.0)
    est_freqs = Calc_stationary_KLD(Compute_matrix(fits, pop))
    return(KLD(true_freqs, est_freqs))
def Inital_fits(freqs, pop):
    '''
    :param freqs: List of Float, the TRUE stationary distribution of codons
    :param pop: Float, the population size (0 <= pop).
    :return: List of Float, the initial fitness values from analytical solution.
    '''
    pi = freqs
    epsilon = np.finfo(float).eps
    Ne = 2 * pop
    pi = [epsilon if x < epsilon else x for x in pi]
    pi = [x/sum(pi) for x in pi]
    pi = [(np.log(x) - np.log(max(pi)) + Ne) / Ne for x in pi]
    fits = np.zeros(20)
    map_count = [0] * 20
    for i in range(61):
        fits[Get_amino_acid(i)] += pi[i]
        map_count[Get_amino_acid(i)] += 1
    for i in range(20):
        fits[i] /= map_count[i]
    return(fits)

def Find_max_freq_AA(codon_distribution):
    '''
    :param codon_distribution: List of Float, the TRUE stationary distribution of codons
    :return: Int, index of the maximum frequency amino acid.
    '''
    AA_freq = [0.0] * 20
    AA_reps = [0] * 20
    for idx in range(0,61):
        AA_freq[Get_amino_acid(idx)] += codon_distribution[idx]
        AA_reps[Get_amino_acid(idx)] += 1

    #norm amino acid frequency based of number of codons
    for idx in range(20):
        AA_freq[idx] /= AA_reps[idx]

    max_val = max(AA_freq)
    max_idx = AA_freq.index(max_val)
    return(max_idx)
def Find_fits(true_distribution, pop, eps=1e-5):
    '''
    :param true_distribution: List of Float, the TRUE stationary distribution of codons
    :param pop: Float, the population size (0 <= pop).
    :param eps: Float, the allowable error in the KLD of the model.
    :return: List of Float, the fitness values from optimized solution.
    '''
    #find the fitness value held constant
    max_freq = Find_max_freq_AA(true_distribution)
    fitness_values = Inital_fits(true_distribution, pop)
    fitness_values = np.delete(fitness_values, max_freq)

    bounds = [(np.finfo(float).eps, 2.0) for n in range(19)]
    fits = sp.optimize.minimize(Frequency_distance, fitness_values, args=(true_distribution, pop, max_freq), bounds=bounds, method='L-BFGS-B', tol = 10e-7)
    fits = np.insert(fits.x, max_freq, 1.0)

    #check solution is accurate enough
    KLDerr = KLD(true_distribution, Calc_stationary_KLD(Compute_matrix(fits, pop)))
    if(KLDerr <= eps):
        return(fits)
    else:
        print("Changing method, checking if a better solution can be found.")
        fits2 = sp.optimize.minimize(Frequency_distance, fitness_values, args=(true_distribution, pop, max_freq), bounds=bounds, method='TNC', tol = 10e-7)
        fits2 = np.insert(fits2.x, max_freq, 1.0)
        KLDerr2 = KLD(true_distribution, Calc_stationary_KLD(Compute_matrix(fits2, pop)))
        if(KLDerr2 <= KLDerr):
            print("Better solution found!")
            return(fits2)
    return(fits)

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

    #Speed of transition does not effect stationary distribtuion
    #the relative rates are all we are concerned with
    #THUS NO MUTATION RATE NEEDED

    #BASIC EXAMPLE
    pop = 1000
    print(Find_fitnesses('table_stationary_distributions2.csv', pop, False))

    #PARALLEL EXAMPLE
    print(Find_fitnesses('table_stationary_distributions3.csv', pop, True))



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
    case11 = False  #Serial Time Case
    case12 = False  #Parallel Time Case
    case13 = False  #Real Challenging Case


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

        est_fits = Find_fits(freqs[0], Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(freqs[0], "\n",est_freqs)
        print("KLD: ", KLD(freqs[0], est_freqs))

    #REAL CASE2
    if(case10):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        Pop = 100

        freqs = freqs[2]
        est_fits = Find_fits(freqs, Pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, Pop))
        print(freqs, "\n",est_freqs)
        print("KLD: ", KLD(freqs, est_freqs))

    #Sequential TIME Test
    if(case11):
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
    if(case12):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        pops = [100]*50
        start = time.time()
        est_fits = Batch_find_fits(freqs, pops, 10)
        end = time.time()
        print("TOTAL TIME: ", (end - start))

    #Challenging real case
    if(case13):
        freqs = Read_frequency_file('table_stationary_distributions.csv')
        pop = 100
        freqs = freqs[24]
        est_fits = Find_fits(freqs, pop)
        est_freqs = Calc_stationary_KLD(Compute_matrix(est_fits, pop))
        print("KLD: ", KLD(freqs, est_freqs))