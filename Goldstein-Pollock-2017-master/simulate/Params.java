package simulate;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Properties;

/**
 * Reads in file called "ParameterList" which includes various relevant parameters
 * Note that it must include one method for generating, choosing or specifying initial sequence
 * 
 * 
 * @author Richard Goldstein
 */
public class Params {
    // Simulation parameters
        public static String PDB_DIRECTORY;         // Location of pdb files
        public static int PROTEIN_SIZE;             // Size of protein
        public static double CUTOFF;                // How close sites have to be to be in 'contact'
        public static double TRANSTRANSRATIO;       // Transition transversion ratio
        public static boolean VERBOSE;              // Print more stuff than you want to look at
        public static boolean PRINT_ALL_MUTATIONS;  // Output stability of all possible mutations
        public static double LOG_UNFOLDED_STATES;   // Log of the number of unfolded states
        public static double N_EFF;                 // Effective population size
        public static int N_GEN;                    // Number of generations to simulate
        public static double MAX_TIME;              // Length of time to simulate
        public static double TEMPERATURE;           // Temperature in kcal/mol (room temp = 0.6)
        public static boolean CREATE_RANDOM_SEQUENCE;   // Create random sequence
        public static boolean CHOOSE_RANDOM_SEQUENCE;   // Choose sequence from set of sequences in file "Sequences"
        public static String START_SEQUENCE;            // Starting sequence (String)
        public static int INITIAL_SEQUENCE;             // Initial sequence from set of sequences in file "Sequences"
        public static Hashtable<Integer, Character> FIXED_SITE_HASH
                = new Hashtable<>();      // Sites and occupants of fixed locations
        
    // Evolutionary dynamics parameters
        
/** Read in from "ParameterList" and set various parameter values 
*/
    public Params() {
        Properties properties = new Properties();
        try {
            FileInputStream inputStream = new FileInputStream("ParameterList");
            properties.load(inputStream);
        } catch (FileNotFoundException e) {
            System.out.println("ParameterList not found");
            System.exit(1);
        } catch (IOException e) {
            System.out.println("SecurityException in ParameterList");
            System.exit(1);
        } catch (SecurityException e) {
            System.out.println("SecurityException in ParameterList");
            System.exit(1);
        }
        
        PDB_DIRECTORY = properties.getProperty("PDB_DIRECTORY", "pdbfiles/");
        PROTEIN_SIZE = Integer.parseInt(properties.getProperty("PROTEIN_SIZE", "300"));
        N_GEN = Integer.parseInt(properties.getProperty("N_GEN", "10"));
        CUTOFF = Double.parseDouble(properties.getProperty("CUTOFF", "7.0"));
        LOG_UNFOLDED_STATES = Double.parseDouble(properties.getProperty("LOG_UNFOLDED_STATES", "368.413615")); // log(1.0E160)
        N_EFF = Double.parseDouble(properties.getProperty("N_EFF", "1.0E6"));
        MAX_TIME = Double.parseDouble(properties.getProperty("MAX_TIME", "1.0"));
        TRANSTRANSRATIO = Double.parseDouble(properties.getProperty("TRANSTRANSRATIO", "2.0"));
        TEMPERATURE = Double.parseDouble(properties.getProperty("TEMPERATURE", "0.6"));
        VERBOSE = Boolean.parseBoolean(properties.getProperty("VERBOSE", "true"));
        PRINT_ALL_MUTATIONS = Boolean.parseBoolean(properties.getProperty("PRINT_ALL_MUTATIONS", "true"));
        CREATE_RANDOM_SEQUENCE = Boolean.parseBoolean(properties.getProperty("CREATE_RANDOM_SEQUENCE", "false"));
        CHOOSE_RANDOM_SEQUENCE = Boolean.parseBoolean(properties.getProperty("CHOOSE_RANDOM_SEQUENCE", "false"));
        START_SEQUENCE = properties.getProperty("START_SEQUENCE", "");
        INITIAL_SEQUENCE = Integer.parseInt(properties.getProperty("INITIAL_SEQUENCE", "-999"));
        String fixedSitesString = properties.getProperty("FIXED_SITES", "");
        String fixedSiteOccupantsString = properties.getProperty("FIXED_SITE_OCCUPANTS", "");
        
        
        // Make sure that there is one and only one method of specifying initial sequence
        int nOptions = 0;
        if (CREATE_RANDOM_SEQUENCE) nOptions++;
        if (CHOOSE_RANDOM_SEQUENCE) nOptions++;
        if (START_SEQUENCE.length()> 0) {
            nOptions++;
            if (START_SEQUENCE.length() != (PROTEIN_SIZE*3)) {
                System.out.println("Starting sequence wrong length");
                System.exit(1);
            }
        }
        if (INITIAL_SEQUENCE >= 0) nOptions++;
        if (nOptions == 0) {
            System.out.println("No initial sequence option choosen");
            System.exit(1);
        }
        if (nOptions > 1) {
            System.out.println(nOptions + "initial sequence options choosen: Only one allowed");
            System.exit(1);
        }

        if (fixedSitesString.length() > 0) {
            String[] fixedSites = fixedSitesString.split(",");
            String[] fixedSiteOccupants = fixedSiteOccupantsString.split(",");
            if (fixedSites.length != fixedSiteOccupants.length) {
                System.out.println("Number of fixed sites does not equal number of fixed site occupants");
                System.exit(1);
            }
            for (int iSite = 0; iSite < fixedSites.length; iSite++) {
                int iSiteLoc = Integer.parseInt(fixedSites[iSite]);
                char occupant = fixedSiteOccupants[iSite].charAt(0);
                FIXED_SITE_HASH.put(iSiteLoc, occupant);
            }
        }
        
    }

}
