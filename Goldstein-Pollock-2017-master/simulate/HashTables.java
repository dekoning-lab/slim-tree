package simulate;

import java.util.Hashtable;

/**
 * Container for mappings from various representations of sequence information
 * 
 * @author Richard Goldstein
 */
public class HashTables {
    public static Hashtable<Character, Integer> aaCharToIntHash = new Hashtable<Character, Integer>();
    public static Hashtable<String, Character> translationHash = new Hashtable<String, Character>();
    public static Hashtable<Character, char[]> mutationHash = new Hashtable<Character, char[]>();
    public static char[] dnaIntToChar = {'A', 'C', 'G', 'T'};
    public static char[] aaIntToChar = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'};

    HashTables() {
        aaCharToIntHash.put('A', 0);
        aaCharToIntHash.put('R', 1);
        aaCharToIntHash.put('N', 2);
        aaCharToIntHash.put('D', 3);
        aaCharToIntHash.put('C', 4);
        aaCharToIntHash.put('Q', 5);
        aaCharToIntHash.put('E', 6);
        aaCharToIntHash.put('G', 7);
        aaCharToIntHash.put('H', 8);
        aaCharToIntHash.put('I', 9);
        aaCharToIntHash.put('L', 10);
        aaCharToIntHash.put('K', 11);
        aaCharToIntHash.put('M', 12);
        aaCharToIntHash.put('F', 13);
        aaCharToIntHash.put('P', 14);
        aaCharToIntHash.put('S', 15);
        aaCharToIntHash.put('T', 16);
        aaCharToIntHash.put('W', 17);
        aaCharToIntHash.put('Y', 18);
        aaCharToIntHash.put('V', 19);
        aaCharToIntHash.put('X', 20);
        
        translationHash.put("TTT", 'F');
        translationHash.put("TTC", 'F');
        translationHash.put("TTA", 'L');
        translationHash.put("TTG", 'L');
        translationHash.put("TCT", 'S');
        translationHash.put("TCC", 'S');
        translationHash.put("TCA", 'S');
        translationHash.put("TCG", 'S');
        translationHash.put("TAT", 'Y');
        translationHash.put("TAC", 'Y');
        translationHash.put("TAA", 'X');
        translationHash.put("TAG", 'X');
        translationHash.put("TGT", 'C');
        translationHash.put("TGC", 'C');
        translationHash.put("TGA", 'X');
        translationHash.put("TGG", 'W');
        translationHash.put("CTT", 'L');
        translationHash.put("CTC", 'L');
        translationHash.put("CTA", 'L');
        translationHash.put("CTG", 'L');
        translationHash.put("CCT", 'P');
        translationHash.put("CCC", 'P');
        translationHash.put("CCA", 'P');
        translationHash.put("CCG", 'P');
        translationHash.put("CAT", 'H');
        translationHash.put("CAC", 'H');
        translationHash.put("CAA", 'Q');
        translationHash.put("CAG", 'Q');
        translationHash.put("CGT", 'R');
        translationHash.put("CGC", 'R');
        translationHash.put("CGA", 'R');
        translationHash.put("CGG", 'R');
        translationHash.put("ATT", 'I');
        translationHash.put("ATC", 'I');
        translationHash.put("ATA", 'I');
        translationHash.put("ATG", 'M');
        translationHash.put("ACT", 'T');
        translationHash.put("ACC", 'T');
        translationHash.put("ACA", 'T');
        translationHash.put("ACG", 'T');
        translationHash.put("AAT", 'N');
        translationHash.put("AAC", 'N');
        translationHash.put("AAA", 'K');
        translationHash.put("AAG", 'K');
        translationHash.put("AGT", 'S');
        translationHash.put("AGC", 'S');
        translationHash.put("AGA", 'R');
        translationHash.put("AGG", 'R');
        translationHash.put("GTT", 'V');
        translationHash.put("GTC", 'V');
        translationHash.put("GTA", 'V');
        translationHash.put("GTG", 'V');
        translationHash.put("GCT", 'A');
        translationHash.put("GCC", 'A');
        translationHash.put("GCA", 'A');
        translationHash.put("GCG", 'A');
        translationHash.put("GAT", 'D');
        translationHash.put("GAC", 'D');
        translationHash.put("GAA", 'E');
        translationHash.put("GAG", 'E');
        translationHash.put("GGT", 'G');
        translationHash.put("GGC", 'G');
        translationHash.put("GGA", 'G');
        translationHash.put("GGG", 'G');
        
        // Representation of transitions and two transversions
        mutationHash.put('A', new char[]{'G', 'C', 'T'});
        mutationHash.put('C', new char[]{'T', 'A', 'G'});
        mutationHash.put('G', new char[]{'A', 'T', 'C'});
        mutationHash.put('T', new char[]{'C', 'G', 'A'});
    }
}
