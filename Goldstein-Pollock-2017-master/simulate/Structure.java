package simulate;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

/**
 * Repository for pdb protein structure information,
 * including checking the structure for obvious problems,
 * making a list of contacts made in the structure,
 * and computing the energy of a sequence in the protein's structure
 * 
 * @author Richard Goldstein
 * 
 * 
 */
public class Structure {

	boolean first = true;
    String proteinName;
    String chainName;
    Vector<double[]> alphaCoordinateVector = new Vector<double[]>();
    Vector<double[]> betaCoordinateVector = new Vector<double[]>();
    InteractionMatrix interactionMatrix;
    Vector<int[]> contactVector = new Vector<int[]>();
    double[][] distanceMatrix = new double[Params.PROTEIN_SIZE][Params.PROTEIN_SIZE];
    int[] sites = {12, 27, 41, 62, 75, 91, 107, 124, 130, 132, 157, 164,
        191, 195, 210, 219, 242, 270, 284, 297};

    Structure(String[] pdbName, InteractionMatrix interactionMatrix, String fileLoc) {
        proteinName = pdbName[0];
        chainName = pdbName[1];
        this.interactionMatrix = interactionMatrix;
        readPDBFile(fileLoc);
        if (!allOK()) {
            System.out.println("Yell!!!");
        }
        makeContactVector();
        if (Params.VERBOSE) {
            System.out.println(proteinName + chainName + " with " + contactVector.size() + " contacts");
        }
    }

    /**
     * Computes energy of an amino acid sequence in the structure of the protein
     */
    double getEnergy(char[] aaSequence) {
        double energy = 0.0;

        for (int[] contact : contactVector) {
            energy += interactionMatrix.getInteraction(aaSequence[contact[0]], aaSequence[contact[1]]);
        }
        return energy;
    }

    /**
     * Find distance between two sites
     */
    double getDistance(int iRes, int jRes) {
        return distanceMatrix[iRes][jRes];
    }

    /**
     * Make a list of all contacts made between sites in the structure
     */
    void makeContactVector() {
        for (int iRes = 0; iRes < Params.PROTEIN_SIZE - 1; iRes++) {
            double[] iCoord = betaCoordinateVector.get(iRes);
            for (int jRes = iRes + 1; jRes < Params.PROTEIN_SIZE; jRes++) {
                double[] jCoord = betaCoordinateVector.get(jRes);
                distanceMatrix[iRes][jRes] = Math.sqrt(((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0])
                                                        + (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1])
                                                        + (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2])));
                distanceMatrix[jRes][iRes] = distanceMatrix[iRes][jRes];
            }
        }
        for (int iRes = 0; iRes < Params.PROTEIN_SIZE - 2; iRes++) {
            for (int jRes = iRes + 2; jRes < Params.PROTEIN_SIZE; jRes++) {
                if (distanceMatrix[iRes][jRes] < Params.CUTOFF) {
                    int[] contact = new int[2];
                    contact[0] = iRes;
                    contact[1] = jRes;
                    contactVector.add(contact);
                }
            }
        }
		
		
		

    }

    /**
     * Check for obvious problems in protein structure
     */
    boolean allOK() {
        boolean ok = true;
        if ((alphaCoordinateVector.size() != betaCoordinateVector.size())
                || (alphaCoordinateVector.size() < Params.PROTEIN_SIZE)
                || (alphaCoordinateVector.size() > Params.PROTEIN_SIZE * 1.5)) {
            System.out.print(proteinName + chainName + "\t");
            System.out.println(alphaCoordinateVector.size() + "\t" + betaCoordinateVector.size());
            ok = false;
        }

        for (int iRes = 0; iRes < alphaCoordinateVector.size() - 1; iRes++) {
            int jRes = iRes + 1;
            double[] iCoord = alphaCoordinateVector.get(iRes);
            double[] jCoord = alphaCoordinateVector.get(jRes);
            double dist2 = ((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0])
                            + (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1])
                            + (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2]));
            if (dist2 > (4.25 * 4.25)) {
                System.out.print("Discontinuity in protein " + proteinName + chainName);
                System.out.println("  Residues " + iRes + " and " + jRes + "\t" + (Math.sqrt(dist2)));
                ok = false;
            }
        }

        double rg = 0.0;
        for (int iRes = 0; iRes < alphaCoordinateVector.size(); iRes++) {
            for (int jRes = 0; jRes < alphaCoordinateVector.size(); jRes++) {
                double[] iCoord = alphaCoordinateVector.get(iRes);
                double[] jCoord = alphaCoordinateVector.get(jRes);
                double dist2 = ((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0])
                                + (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1])
                                + (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2]));
                rg += dist2;
            }
        }
        rg /= (2.0 * alphaCoordinateVector.size() * alphaCoordinateVector.size());
        if (rg > 650) {
            System.out.println("Problem with rg " + proteinName + chainName + "\t" + alphaCoordinateVector.size() + "\t" + rg);
            ok = false;
        }

        double nContacts = 0.0;
        for (int iRes = 0; iRes < betaCoordinateVector.size(); iRes++) {
            for (int jRes = 0; jRes < betaCoordinateVector.size(); jRes++) {
                double[] iCoord = betaCoordinateVector.get(iRes);
                double[] jCoord = betaCoordinateVector.get(jRes);
                double dist2 = ((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0])
                                + (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1])
                                + (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2]));
                if (dist2 < 36.0) {
                    nContacts++;
                }
            }
        }

        nContacts /= betaCoordinateVector.size();
        if ((nContacts < 5) || (nContacts > 7)) {
            System.out.println("Problem with nContacts " + proteinName + chainName + "\t" + alphaCoordinateVector.size() + "\t" + nContacts);
            ok = false;
        }
        return ok;
    }

    /**
     * Read in data from pdb file
     */
    void readPDBFile(String fileLoc) {
        boolean ok = true;
        try {
            FileReader file = new FileReader((fileLoc));
            BufferedReader buff = new BufferedReader(file);
            int prevAlpha = -999;
            int prevBeta = -999;
            boolean eof = false;
            while (!eof) {
                String line = buff.readLine();
                if (line == null) {
                    eof = true;
                } else {
                    if ((line.substring(0, 4).equals("ATOM"))
                            && (line.substring(21, 22).equals(chainName))
                            && (line.substring(16, 17).equals(" ") || line.substring(16, 17).equals("A"))) {
                        String atomType = line.substring(13, 15);
                        String resType = line.substring(17, 20);
                        if (atomType.equals("CA")) {
                            int currentAlpha = Integer.parseInt(line.substring(22, 26).replace(" ", ""));
                            if (prevAlpha < -100) {
                                prevAlpha = currentAlpha - 1;
                            }
                            if (currentAlpha != (prevAlpha + 1)) {
                                System.out.println("Problem with Alpha " + proteinName + chainName + ", residue " + currentAlpha + "	" + prevAlpha);
                                System.out.println(line);
                                ok = false;
                                System.exit(1);
                            }
                            prevAlpha = currentAlpha;
                            double[] coords = new double[3];
                            coords[0] = Double.parseDouble(line.substring(30, 38));
                            coords[1] = Double.parseDouble(line.substring(38, 46));
                            coords[2] = Double.parseDouble(line.substring(46, 54));
                            alphaCoordinateVector.add(coords);
                            if (resType.equals("GLY")) {
                                betaCoordinateVector.add(coords);
                                prevBeta++;
                            }
                        } else if (atomType.equals("CB")) {
                            int currentBeta = Integer.parseInt(line.substring(22, 26).replace(" ", ""));
                            if (prevBeta < -100) {
                                prevBeta = currentBeta - 1;
                            }
                            if (currentBeta != (prevBeta + 1)) {
                                System.out.println("Problem with Beta " + proteinName + chainName + ", residue " + currentBeta + "	" + prevBeta);
                                System.out.println(line);
                                ok = false;
                                System.exit(1);
                            }
                            prevBeta = currentBeta;
                            double[] coords = new double[3];
                            coords[0] = Double.parseDouble(line.substring(30, 38));
                            coords[1] = Double.parseDouble(line.substring(38, 46));
                            coords[2] = Double.parseDouble(line.substring(46, 54));
                            betaCoordinateVector.add(coords);
                        }
                    }
                }
            }
            buff.close();
            file.close();
        } catch (IOException ioe) {
            System.out.println("Error -- " + ioe.toString());
            System.exit(1);
        }
    }
    
}
