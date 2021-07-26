package simulate;

/**
 * Interaction potential from 
 * Miyazawa and Jernigan Macromolecules 1985, 18, 534-552
 * 
 * 
 * @author Richard Goldstein
 */
public class InteractionMatrix {
  
    double getInteraction(char resA, char resB) {
        return (interactionMatrix[HashTables.aaCharToIntHash.get(resA)][HashTables.aaCharToIntHash.get(resB)]);
    }
    
    

    double[][] interactionMatrix =
        {{-0.13,  0.43,  0.28,  0.12,  0.00,  0.08,  0.26, -0.07,  0.34, -0.22, -0.01,  0.14,  0.25,  0.03,  0.10, -0.06, -0.09, -0.09,  0.09, -0.10,  0.00},
         { 0.43,  0.11, -0.14, -0.72,  0.24, -0.52, -0.74, -0.04, -0.12,  0.42,  0.35,  0.75,  0.31,  0.41, -0.38,  0.17, -0.35, -0.16, -0.25,  0.30,  0.00},
         { 0.28, -0.14, -0.53, -0.30,  0.13, -0.25, -0.32, -0.14, -0.24,  0.53,  0.30, -0.33,  0.08,  0.18, -0.18, -0.14, -0.11,  0.06, -0.20,  0.50,  0.00},
         { 0.12, -0.72, -0.30,  0.04,  0.03, -0.17, -0.15, -0.22, -0.39,  0.59,  0.67, -0.76,  0.65,  0.39,  0.04, -0.31, -0.29,  0.24,  0.00,  0.58,  0.00},
         { 0.00,  0.24,  0.13,  0.03, -0.20,  0.05,  0.69, -0.08, -0.19,  0.16, -0.08,  0.71,  0.19, -0.23,  0.00, -0.02,  0.19,  0.08,  0.04,  0.06,  0.00},
         { 0.08, -0.52, -0.25, -0.17,  0.05,  0.29, -0.17, -0.06, -0.02,  0.36,  0.26, -0.38,  0.46,  0.49, -0.42, -0.14, -0.14,  0.08, -0.20,  0.24,  0.00},
         { 0.26, -0.74, -0.32, -0.15,  0.69, -0.17, -0.03,  0.25, -0.45,  0.35,  0.43, -0.97,  0.44,  0.27, -0.10, -0.26,  0.00,  0.29, -0.10,  0.34,  0.00},
         {-0.07, -0.04, -0.14, -0.22, -0.08, -0.06,  0.25, -0.38,  0.20,  0.25,  0.23,  0.11,  0.19,  0.38, -0.11, -0.16, -0.26,  0.18,  0.14,  0.16,  0.00},
         { 0.34, -0.12, -0.24, -0.39, -0.19, -0.02, -0.45,  0.20, -0.29,  0.49,  0.16,  0.22,  0.99, -0.16, -0.21, -0.05, -0.19, -0.12, -0.34,  0.19,  0.00},
         {-0.22,  0.42,  0.53,  0.59,  0.16,  0.36,  0.35,  0.25,  0.49, -0.22, -0.41,  0.36, -0.28, -0.19,  0.25,  0.21,  0.14,  0.02,  0.11, -0.25,  0.00},
         {-0.01,  0.35,  0.30,  0.67, -0.08,  0.26,  0.43,  0.23,  0.16, -0.41, -0.27,  0.19, -0.20, -0.30,  0.42,  0.25,  0.20, -0.09,  0.24, -0.29,  0.00},
         { 0.14,  0.75, -0.33, -0.76,  0.71, -0.38, -0.97,  0.11,  0.22,  0.36,  0.19,  0.25,  0.00,  0.44,  0.11, -0.13, -0.09,  0.22, -0.21,  0.44,  0.00},
         { 0.25,  0.31,  0.08,  0.65,  0.19,  0.46,  0.44,  0.19,  0.99, -0.28, -0.20,  0.00,  0.04, -0.42, -0.34,  0.14,  0.19, -0.67, -0.13, -0.14,  0.00},
         { 0.03,  0.41,  0.18,  0.39, -0.23,  0.49,  0.27,  0.38, -0.16, -0.19, -0.30,  0.44, -0.42, -0.44,  0.20,  0.29,  0.31, -0.16,  0.00, -0.22,  0.00},
         { 0.10, -0.38, -0.18,  0.04,  0.00, -0.42, -0.10, -0.11, -0.21,  0.25,  0.42,  0.11, -0.34,  0.20,  0.26,  0.01, -0.07, -0.28, -0.33,  0.09,  0.00},
         {-0.06,  0.17, -0.14, -0.31, -0.02, -0.14, -0.26, -0.16, -0.05,  0.21,  0.25, -0.13,  0.14,  0.29,  0.01, -0.20, -0.08,  0.34,  0.09,  0.18,  0.00},
         {-0.09, -0.35, -0.11, -0.29,  0.19, -0.14,  0.00, -0.26, -0.19,  0.14,  0.20, -0.09,  0.19,  0.31, -0.07, -0.08,  0.03,  0.22,  0.13,  0.25,  0.00},
         {-0.09, -0.16,  0.06,  0.24,  0.08,  0.08,  0.29,  0.18, -0.12,  0.02, -0.09,  0.22, -0.67, -0.16, -0.28,  0.34,  0.22, -0.12, -0.04, -0.07,  0.00},
         { 0.09, -0.25, -0.20,  0.00,  0.04, -0.20, -0.10,  0.14, -0.34,  0.11,  0.24, -0.21, -0.13,  0.00, -0.33,  0.09,  0.13, -0.04, -0.06,  0.02,  0.00},
         {-0.10,  0.30,  0.50,  0.58,  0.06,  0.24,  0.34,  0.16,  0.19, -0.25, -0.29,  0.44, -0.14, -0.22,  0.09,  0.18,  0.25, -0.07,  0.02, -0.29,  0.00},
         { 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00}};
    
  
}
