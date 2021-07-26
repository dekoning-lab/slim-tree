package simulate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

/**
 * Reads in set of sequences from file called 'Sequences'
 * that are acceptable as starting sequences
 * This sequences have evolved to a reasonable degree of stability
 * Also contains code for creating a random dna sequence without stop codons
 * 
 * 
 * @author Richard Goldstein
 */
public class Sequences {
    static public Vector<String> sequenceVector = new Vector<String>();

    Sequences() {
        File sequenceFile = new File("Sequences");
        if (sequenceFile.exists()) {
            try {
                FileReader file = new FileReader(sequenceFile);
                BufferedReader buff = new BufferedReader(file);
                boolean eof = false;
                while (!eof) {
                    String line = buff.readLine();
                    if (line == null) {
                        eof = true;
                    } else {
                        sequenceVector.add(line);
                    }
                }
            } catch (IOException ioe) {
                System.out.println("Error -- " + ioe.toString());
                System.exit(1);
            }
        }
    }

    /**
     * Get specific sequence from list of available sequences
     * 
     */
    static public String getSequence(int iSeq) {
        if (iSeq > sequenceVector.size()-1) {
            System.out.println("Request of sequence that does not exist");
            System.exit(1);
        }
        return sequenceVector.get(iSeq);
    }
    
    static public int getSequenceListSize() {
        return sequenceVector.size();
    }

    /**
     * Generate random sequence with no stop codons
     * 
     */    
    static public String getRandomSequence() {
        char[][] dnaSequence = new char[Params.PROTEIN_SIZE][3];
        String randomSequence = "";
        for (int iSite = 0; iSite < Params.PROTEIN_SIZE; iSite++) {
            boolean ok = false;
            while (!ok) {
                for (int iBase = 0; iBase < 3; iBase++) {
                    dnaSequence[iSite][iBase] = HashTables.dnaIntToChar[Simulate.random.nextInt(4)];
                }
                String codon = new String(dnaSequence[iSite]);
                if (HashTables.translationHash.get(codon) != 'X') {
                    ok = true;
                }
            }
            randomSequence = randomSequence + new String(dnaSequence[iSite]);
        }
        return randomSequence;
    }

}
