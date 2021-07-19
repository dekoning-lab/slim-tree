package simulate;

import java.util.Hashtable;
import java.util.Vector;

/**
 *
 * @author Richard Goldstein
 */
public class Protein {
    
    char[] aaSequence = new char[Params.PROTEIN_SIZE];          // amino acid sequence, char[] format
    char[][] dnaSequence = new char[Params.PROTEIN_SIZE][3];    // dna sequence, in array of codons
    String proteinSequence = "";                                // aa sequence, String format

    double nGen = 0.0;          // number of generations
    int nSubstitutions = 0;     // number of aa substitutions
    double time = 0.0;          // simulation time in units of dna substitutions per site under neutral conditions

    double nativeDeltaG = 0.0;
    double nativePFolded = 0.0;
    
    double probTransition = Params.TRANSTRANSRATIO / (2.0 + Params.TRANSTRANSRATIO);    // Probability of transition
    double probTransversion = 1.0 / (2.0 + Params.TRANSTRANSRATIO);                     // Probability of transversion
    double[] pNeutral = new double[3];                                                  // Probability of three types of mutations

    Structure nativeStructure;
    Vector<Structure> dummyStructureVector;     // Vector of structures representing unfolded conformationsn
    Hashtable<Integer, Character> fixedResidues = null;     // List of residues fixed during the simulation

    Hashtable<String, Double> previousSequenceHash = new Hashtable<String, Double>();   // Storage for stabilities of recent sequences

    /**
     * Initialise Protein when an initial (DNA) sequence is provided
     * 
     */
    Protein(Structure nativeStructure, Vector<Structure> dummyStructureVector, String startSequence,
            Hashtable<Integer, Character> fixedResidues) {
        this.nativeStructure = nativeStructure;
        this.dummyStructureVector = dummyStructureVector;
        this.fixedResidues = fixedResidues;
        pNeutral[0] = probTransition;
        pNeutral[1] = probTransversion;
        pNeutral[2] = probTransversion;
        for (int iCodon = 0; iCodon < Params.PROTEIN_SIZE; iCodon++) {
            dnaSequence[iCodon][0] = startSequence.charAt(3 * iCodon);
            dnaSequence[iCodon][1] = startSequence.charAt(3 * iCodon + 1);
            dnaSequence[iCodon][2] = startSequence.charAt(3 * iCodon + 2);
        }
        for (int iRes = 0; iRes < Params.PROTEIN_SIZE; iRes++) {
            aaSequence[iRes] = HashTables.translationHash.get(new String(dnaSequence[iRes]));
        }
        for (int iRes : Params.FIXED_SITE_HASH.keySet()) {
            aaSequence[iRes] = Params.FIXED_SITE_HASH.get(iRes);
        }
        proteinSequence = new String(aaSequence);;
		
		// if(Params.VERBOSE){
			// System.out.println("hhh\t" + proteinSequence);
			// nativeDeltaG = computeDeltaG(aaSequence);
			// nativePFolded = computePFolded(nativeDeltaG);
			// System.out.printf("sss\t%6.3f\t%6.3f\t%10.5g\t\n", nativeDeltaG, nativePFolded, (1.0 - nativePFolded));
                
			// System.out.print("jjj\t");
			// for (int iRes = 0; iRes < Params.PROTEIN_SIZE; iRes++) {
				// System.out.print(new String(dnaSequence[iRes]));
			// }
			// System.out.println();
		// }
    }
    
    /**
     * Initialise Protein when the number of a (DNA) sequence is provided, corresponding to
     * a sequence from the file 'Sequences'
     * 
     */
    Protein(Structure nativeStructure, Vector<Structure> dummyStructureVector, int initialSequence,
            Hashtable<Integer, Character> fixedResidues) {
        this(nativeStructure, dummyStructureVector, Sequences.getSequence(initialSequence), fixedResidues);
    }
    
    /**
     * If no sequence string or integer reference to Sequences is provided, 
     * initialise Protein with a random sequence
     * 
     */
    Protein(Structure nativeStructure, Vector<Structure> dummyStructureVector, 
            Hashtable<Integer, Character> fixedResidues) {
        this(nativeStructure, dummyStructureVector, Sequences.getRandomSequence(), fixedResidues);
    }
    
    /**
     * Print current sequences and (if PRINT_ALL_MUTATIONS) effects of all mutations
     * 
     */
    void summarise() {
        // System.out.println("hhh\t" + proteinSequence);
        // System.out.print("jjj\t");
        for (int iRes = 0; iRes < Params.PROTEIN_SIZE; iRes++) {
            System.out.print(new String(dnaSequence[iRes]));
        }
        // System.out.println();
        // if (Params.PRINT_ALL_MUTATIONS) {
            // for (int iFocal = 0; iFocal < 300; iFocal++) {
                // char oldAA = aaSequence[iFocal];
                // System.out.print("kkk\t" + iFocal + "\t" + oldAA + "\t");
                // for (int iAA = 0; iAA < 21; iAA++) {
                    // aaSequence[iFocal] = HashTables.aaIntToChar[iAA]; // aaSequence is now mutated 
                    // double newDeltaG = computeDeltaG(aaSequence);
                    // System.out.format("\t%1$6.3f", newDeltaG);
                // }
                // aaSequence[iFocal] = oldAA;
                // System.out.println();
            // }
        // }
    }

    double getTime() {
        return time;
    }

    /**
     * Evaluate all possible mutations, choose one based on probability of its occurrence,
     * and advance clock appropriately
     * 
     */
    void mutate(int iGen) {
        // Calculate substitution relative probability for all possible mutations
        double[][][] substitutionProb = new double[Params.PROTEIN_SIZE][3][3];
        double substitutionProbTotal = 0.0;
        for (int iRes = 0; iRes < Params.PROTEIN_SIZE; iRes++) {        // Consider all sites
            if (!Params.FIXED_SITE_HASH.containsKey(iRes)) {                     // Do not change fixed residue
                char oldAA = aaSequence[iRes];
                char[] oldCodon = dnaSequence[iRes];
                char[] newCodon = new char[3];
                newCodon[0] = oldCodon[0];
                newCodon[1] = oldCodon[1];
                newCodon[2] = oldCodon[2];
                for (int iBase = 0; iBase < 3; iBase++) {               // Consider all bases
                    for (int iChoose = 0; iChoose < 3; iChoose++) {     // Consider all possible changes
                        substitutionProb[iRes][iBase][iChoose] = pNeutral[iChoose];
                        newCodon[iBase] = HashTables.mutationHash.get(oldCodon[iBase])[iChoose];
                        char newAA = HashTables.translationHash.get(new String(newCodon));
                        aaSequence[iRes] = newAA;
                        String mutantProteinSequence = new String(aaSequence);      // New mutated sequence
                        if (newAA == 'X') {
                            substitutionProb[iRes][iBase][iChoose] = 0.0;               // Prob of acceptance for stop codons = 0
                        } else if (newAA != oldAA) {                                // Non synonymous mutation
                            double mutantPFolded = computePFolded(computeDeltaG(aaSequence));
                            substitutionProb[iRes][iBase][iChoose] *= computePFix(mutantPFolded);
                        }
                        substitutionProbTotal += substitutionProb[iRes][iBase][iChoose];
                        newCodon[iBase] = oldCodon[iBase];
                        aaSequence[iRes] = oldAA;
                    }
                }
            }
        }

        // Advance clock by appropriate amount and choose substitution to accept
        int iRes = 0;
        int iBase = 0;
        int iChoose = 0;
        time += -Math.log(1.0 - Simulate.random.nextDouble()) / substitutionProbTotal;
        double randomNumber = Simulate.random.nextDouble();
        boolean chosen = false;
        while (!chosen) {
            substitutionProb[iRes][iBase][iChoose] /= substitutionProbTotal;
            if (substitutionProb[iRes][iBase][iChoose] >= randomNumber) {
                chosen = true;
            } else {
                randomNumber -= substitutionProb[iRes][iBase][iChoose];
                iChoose++;
                if (iChoose == 3) {
                    iChoose = 0;
                    iBase++;
                }
                if (iBase == 3) {
                    iBase = 0;
                    iRes++;
                }
                if (iRes == Params.PROTEIN_SIZE) {
                    iChoose = 0;
                    iBase = 0;
                    iRes = 0;
                    randomNumber = Simulate.random.nextDouble();
                }
            }
        }

        // Update sequence
        char oldNuc = dnaSequence[iRes][iBase];
        char newNuc = HashTables.mutationHash.get(dnaSequence[iRes][iBase])[iChoose];
        dnaSequence[iRes][iBase] = newNuc;
        char oldAA = aaSequence[iRes];
        char newAA = HashTables.translationHash.get(new String(dnaSequence[iRes]));

        // Update and print event
        boolean synonymous = true;
        if (newAA != oldAA) {
            synonymous = false;
            aaSequence[iRes] = newAA;
            proteinSequence = new String(aaSequence);
            nSubstitutions++;
            double mutantDeltaG = computeDeltaG(aaSequence);
            double mutantPFolded = computePFolded(mutantDeltaG);
			if (Params.VERBOSE) {
				System.out.format("xxx\t%1$d\t%2$d\t%3$10.4f\t%4$c%5$03d%6$c\t%7$c%8$1d%9$c\t%10$6.3f\t%11$6.3f\n",
                        iGen, nSubstitutions, time, oldAA, iRes, newAA, oldNuc, iBase, newNuc, nativeDeltaG, mutantDeltaG);
			}
            nativeDeltaG = mutantDeltaG;
            nativePFolded = mutantPFolded;
            // summarise();
        } else if (Params.VERBOSE){
            System.out.format("yyy\t%1$d\t%2$d\t%3$10.4f\t%4$c%5$03d%6$c\t%7$c%8$1d%9$c\t%10$6.3f\t%11$6.3f\n",
                    iGen, nSubstitutions, time, oldAA, iRes, newAA, oldNuc, iBase, newNuc, nativeDeltaG, nativeDeltaG);
        }
    }

    
    /**
     * Compute probability protein sequence would be folded at equilibrium
     * 
     */
    double computePFolded(double deltaG) {
        double factor = Math.exp(-deltaG/Params.TEMPERATURE);
        return (factor / (1 + factor));
    }

    /**
     * Compute probability of fixation of a mutation relative to neutral following Kimura 1969
     * 
     */
    double computePFix(double newPFolded) {
        double s = (newPFolded - nativePFolded) / nativePFolded;
        double prob = 1.0;
        if ((Math.abs(s) * Params.N_EFF) > 0.001) {
            prob = 2.0 * Params.N_EFF * (1.0 - Math.exp(-2.0 * s)) / (1.0 - Math.exp(-4.0 * s * Params.N_EFF));
        }
        return prob;
    }

    /**
     * Compute stability of protein
     * 
     */
    double computeDeltaG(char[] currentAaSequence) {
        String currentProteinSequence = new String(currentAaSequence);
        if (previousSequenceHash.containsKey(currentProteinSequence)) {
            return previousSequenceHash.get(currentProteinSequence);
        } else {
            double nativeEnergy = nativeStructure.getEnergy(currentAaSequence);
            double avgDummyEnergy = 0.0;
            double avgDummyEnergy2 = 0.0;
            for (Structure dummyStructure : dummyStructureVector) {
                double dummyEnergy = dummyStructure.getEnergy(currentAaSequence);
				avgDummyEnergy += dummyEnergy;
                avgDummyEnergy2 += (dummyEnergy * dummyEnergy);
            }
            avgDummyEnergy /= dummyStructureVector.size();
            avgDummyEnergy2 /= dummyStructureVector.size();
            double sigma2 = avgDummyEnergy2 - (avgDummyEnergy * avgDummyEnergy);
			double deltaG = nativeEnergy - avgDummyEnergy + (Params.TEMPERATURE * Params.LOG_UNFOLDED_STATES) 
                    + (sigma2 / (2.0 * Params.TEMPERATURE));
            if (previousSequenceHash.size() > 100000) {
                previousSequenceHash.clear();
            }
            previousSequenceHash.put(currentProteinSequence, deltaG);
            return deltaG;
        }

    }

}
