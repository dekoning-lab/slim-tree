package simulate;

import java.io.*;
import java.util.*;

/**
 *
 * @author Richard Goldstein
 */
public class Simulate {
    
    Protein protein;
    Structure nativeStructure;
    Vector<Structure> dummyStructureVector = new Vector<Structure>();
    Vector<Integer> fixedSiteVector = new Vector<Integer>();
    Hashtable<Integer, Character> fixedResidues = new Hashtable<Integer, Character>();
    String startingSequence = "";
    Params params = new Params();
    Sequences sequences = new Sequences();
    HashTables hashTables = new HashTables();
    public static Random random = new Random();
	String[][] dummyProteinList;
	String[] nativeProtein = new String [2];
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Simulate simulate = new Simulate(args);
        simulate.run();
    }

    Simulate(String[] arguments) {
		//Not in original code - get location of proteins and dummy proteins
		String protLoc = arguments[0];
		String dummyProtLoc = arguments[1];
		
		//Not in original code - get native protein structure
		String nativeProt = arguments[2];
		nativeProtein = nativeProt.split(",");

		//Not in original code - convert arguments into dummy protein list
		String dummyProts = arguments[3];
		String [] dummyProteinList1d = dummyProts.split(",");
		
		dummyProteinList = new String[dummyProteinList1d.length/2][2];
		
		for(int pos = 0; pos < dummyProteinList1d.length; pos ++){
			dummyProteinList[pos/2][pos%2] = dummyProteinList1d[pos];
		}
		
		//In original code with addition of new protein locations:
		
        InteractionMatrix interactionMatrix = new InteractionMatrix();
        nativeStructure = new Structure(nativeProtein, interactionMatrix, protLoc);
        for (int iDummy = 0; iDummy < dummyProteinList.length; iDummy++) {
            Structure dummyStructure = new Structure(dummyProteinList[iDummy], interactionMatrix, 
					dummyProtLoc + "/" + dummyProteinList[iDummy][0] + ".pdb");
            dummyStructureVector.add(dummyStructure);
        }
    }
	

    void run() {
        if (Params.CREATE_RANDOM_SEQUENCE) {
            protein = new Protein(nativeStructure, dummyStructureVector, fixedResidues);
        } else if (Params.CHOOSE_RANDOM_SEQUENCE) {
            if (Sequences.getSequenceListSize() == 0) {
                System.out.println("No available sequences: Does file Sequences exist?");
                System.exit(1);
            }
            int initialSequence = new Random().nextInt(Sequences.getSequenceListSize());
            protein = new Protein(nativeStructure, dummyStructureVector, initialSequence, fixedResidues);
        } else if (Params.INITIAL_SEQUENCE >= 0) {
            if (Sequences.getSequenceListSize() < Params.INITIAL_SEQUENCE+1) {
                System.out.println("Sequence not available: Does file Sequences exist?");
                System.exit(1);
            }
            protein = new Protein(nativeStructure, dummyStructureVector, Params.INITIAL_SEQUENCE, fixedResidues);
        } else if (Params.START_SEQUENCE.length() > 0) {
            protein = new Protein(nativeStructure, dummyStructureVector, Params.START_SEQUENCE, fixedResidues);
        }

        boolean finished = false;
        int iGen = 0;
        while (!finished) {
            protein.mutate(iGen);
            iGen++;
            // if ((iGen > Params.N_GEN) || (protein.getTime() > Params.MAX_TIME)) {
                // finished = true;
				// protein.summarise();
            // }
			if(protein.nativePFolded >= 0.99){
				finished = true;
				protein.summarise();
			}
        }
    }

}
