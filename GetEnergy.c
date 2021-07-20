#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

/* Program to calculate fitness of a series of individuals from their amino acid sequences
based on the methods by Pollock et al. (2012) */

//Miyazawa-Jernigan matrix in RT units - 0 rows/columns are from letters of the alphabet not represented by amino acids
const double miyazawa_jernigan [26][26] = {
		-0.13,0,0,0.12,0.26,0.03,-0.07,0.34,-0.22,0,0.14,-0.01,0.25,0.28,0,0.1,0.08,0.43,-0.06,-0.09,0,-0.1,-0.09,0,0.09,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0,0,-0.2,0.03,0.69,-0.23,-0.08,-0.19,0.16,0,0.71,-0.08,0.19,0.13,0,0,0.05,0.24,-0.02,0.19,0,0.06,0.08,0,0.04,0,
		0.12,0,0.03,0.04,-0.15,0.39,-0.22,-0.39,0.59,0,-0.76,0.67,0.65,-0.3,0,0.04,-0.17,-0.72,-0.31,-0.29,0,0.58,0.24,0,0,0,
		0.26,0,0.69,-0.15,-0.03,0.27,0.25,-0.45,0.35,0,-0.97,0.43,0.44,-0.32,0,-0.1,-0.17,-0.74,-0.26,0,0,0.34,0.29,0,-0.1,0,
		0.03,0,-0.23,0.39,0.27,-0.44,0.38,-0.16,-0.19,0,0.44,-0.3,-0.42,0.18,0,0.2,0.49,0.41,0.29,0.31,0,-0.22,-0.16,0,0,0,
		-0.07,0,-0.08,-0.22,0.25,0.38,-0.38,0.2,0.25,0,0.11,0.23,0.19,-0.14,0,-0.11,-0.06,-0.04,-0.16,-0.26,0,0.16,0.18,0,0.14,0,
		0.34,0,-0.19,-0.39,-0.45,-0.16,0.2,-0.29,0.49,0,0.22,0.16,0.99,-0.24,0,-0.21,-0.02,-0.12,-0.05,-0.19,0,0.19,-0.12,0,-0.34,0,
		-0.22,0,0.16,0.59,0.35,-0.19,0.25,0.49,-0.22,0,0.36,-0.41,-0.28,0.53,0,0.25,0.36,0.42,0.21,0.14,0,-0.25,0.02,0,0.11,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.14,0,0.71,-0.76,-0.97,0.44,0.11,0.22,0.36,0,0.25,0.19,0,-0.33,0,0.11,-0.38,0.75,-0.13,-0.09,0,0.44,0.22,0,-0.21,0,
		-0.01,0,-0.08,0.67,0.43,-0.3,0.23,0.16,-0.41,0,0.19,-0.27,-0.2,0.3,0,0.42,0.26,0.35,0.25,0.2,0,-0.29,-0.09,0,0.24,0,
		0.25,0,0.19,0.65,0.44,-0.42,0.19,0.99,-0.28,0,0,-0.2,0.04,0.08,0,-0.34,0.46,0.31,0.14,0.19,0,-0.14,-0.67,0,-0.13,0,
		0.28,0,0.13,-0.3,-0.32,0.18,-0.14,-0.24,0.53,0,-0.33,0.3,0.08,-0.53,0,-0.18,-0.25,-0.14,-0.14,-0.11,0,0.5,0.06,0,-0.2,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.1,0,0,0.04,-0.1,0.2,-0.11,-0.21,0.25,0,0.11,0.42,-0.34,-0.18,0,0.26,-0.42,-0.38,0.01,-0.07,0,0.09,-0.28,0,-0.33,0,
		0.08,0,0.05,-0.17,-0.17,0.49,-0.06,-0.02,0.36,0,-0.38,0.26,0.46,-0.25,0,-0.42,0.29,-0.52,-0.14,-0.14,0,0.24,0.08,0,-0.2,0,
		0.43,0,0.24,-0.72,-0.74,0.41,-0.04,-0.12,0.42,0,0.75,0.35,0.31,-0.14,0,-0.38,-0.52,0.11,0.17,-0.35,0,0.3,-0.16,0,-0.25,0,
		-0.06,0,-0.02,-0.31,-0.26,0.29,-0.16,-0.05,0.21,0,-0.13,0.25,0.14,-0.14,0,0.01,-0.14,0.17,-0.2,-0.08,0,0.18,0.34,0,0.09,0,
		-0.09,0,0.19,-0.29,0,0.31,-0.26,-0.19,0.14,0,-0.09,0.2,0.19,-0.11,0,-0.07,-0.14,-0.35,-0.08,0.03,0,0.25,0.22,0,0.13,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		-0.1,0,0.06,0.58,0.34,-0.22,0.16,0.19,-0.25,0,0.44,-0.29,-0.14,0.5,0,0.09,0.24,0.3,0.18,0.25,0,-0.29,-0.07,0,0.02,0,
		-0.09,0,0.08,0.24,0.29,-0.16,0.18,-0.12,0.02,0,0.22,-0.09,-0.67,0.06,0,-0.28,0.08,-0.16,0.34,0.22,0,-0.07,-0.12,0,-0.04,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.09,0,0.04,0,-0.1,0,0.14,-0.34,0.11,0,-0.21,0.24,-0.13,-0.2,0,-0.33,-0.2,-0.25,0.09,0.13,0,0.02,-0.04,0,-0.06,0,
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};





int *sequence_int;
int sequence_length;

int max_contacts;
int max_line_size;
int num_diverse_prots;

int ** diverse_contacts;
double * energies;
double* main_energies;
double main_p_fold;

const float kT = 0.6; //Conversion to account for RT units
long double log_N_U;
bool second = true;


/* Returns the contacts of a particular protein in a 2D array from the line 
read from the csv file of the 189 structurally diverse proteins */
int * get_protein_contacts(char* line)
{
	int* contacts = calloc(max_contacts*2, sizeof(int));
	
    char* split_string;
	int pos = 0;
	int contact;

	//Read each contact from the csv and add to array
    for (split_string = strtok(line, ",");
            split_string && *split_string;
            split_string = strtok(NULL, ",\n"))
    {
    	contact = atoi(split_string);

    	//Only look at positions up to the sequence size
    	if (contact >= sequence_length){
    		break;
    	}

    	contacts[pos] = contact;
		pos++;
    }
    return contacts;
}


/* Reads a csv data file containing the contacts of a set of proteins
and adds the contacts of these proteins to a 3D array with first dimension being each 
individual protein, second dimension being the rows of contacts and 3rd dimension being the 
2 different contacting proteins */
int ** read_dat(char* filename, int num_prots)
{
    FILE* stream = fopen(filename, "r");

    char line[max_line_size]; // Make large enough for all sequences

	int count = 0;


	long toAllocate = num_prots*sizeof(int)*max_contacts;
	int** protein_contact_maps = malloc(toAllocate);


	//Gets the contacts of each protein
    while (fgets(line, max_line_size, stream))
    {
        char* tmp = strdup(line);
		int* contact_array = get_protein_contacts(tmp);
		protein_contact_maps[count] = contact_array;
        free(tmp);
		count ++;
    }

	fclose(stream);
	return protein_contact_maps;
}



/* Finds the value of G_NS by computing the sum of Miyazawa-Jernigan values for all amino acids in contact. 
Second size of array is 2 dimensional as the array will always have only 2 amino acids in contact*/
double find_G_NS(int * contact_matrix, int protein_num)
{


	double contact_energies = 0;
	double energy;
	int contact_num = 0;


	int contact_one = contact_matrix[contact_num];
	int contact_two = contact_matrix[contact_num+1];

	do {
		// Find contact energy and add to total, as well as list of all energies
		energy = miyazawa_jernigan[sequence_int[contact_one]][sequence_int[contact_two]];
		contact_energies += energy;

		energies[(contact_num/2)+protein_num*max_contacts] = energy;

		contact_num += 2;

		contact_one = contact_matrix[contact_num];
		contact_two = contact_matrix[contact_num+1];


	} while (contact_one != 0 || contact_two != 0); //Both having zero marks end of array

	main_energies[protein_num] = contact_energies;
	return contact_energies;
	
}


/* Finds the value of G_NS by computing the sum of Miyazawa-Jernigan values for amino acids that have
 * changed sine the first position*/
double find_G_NS_not_first(int * mutations, int * contact_matrix, int protein_num)
{
	double contact_energies = main_energies[protein_num];
	int contact_num = 0;


	int contact_one = contact_matrix[contact_num];
	int contact_two = contact_matrix[contact_num+1];

	int dif_contact_pointer_1 = -2;
	int dif_contact_pointer_2 = -2;
	int new_contact_1 = -1;
	int new_contact_2 = -1;

	float change_in_energy;
	do {
		//Check to see if contact one has changed
		while ((new_contact_1 < contact_one)){
			dif_contact_pointer_1 += 2;
			new_contact_1 = mutations[dif_contact_pointer_1];
			if((new_contact_2 == 0) & (dif_contact_pointer_2 != 0)){
					break;
			}
		} 
		//Check to see if contact 2 has changed
		while ((new_contact_2 != contact_two)){
			dif_contact_pointer_2 += 2;
			new_contact_2 = mutations[dif_contact_pointer_2];
			if ((new_contact_2 == 0) & (dif_contact_pointer_2 != 0)){
				new_contact_2 = -1;
				break;
			}
		}
		
		// If either of contact 1 or contact 2 have changed - update with correct energy
		if ((new_contact_1 == contact_one) & (new_contact_2 == contact_two)){
			change_in_energy = (miyazawa_jernigan[mutations[dif_contact_pointer_1+1]][mutations[dif_contact_pointer_2+1]]-
					energies[(contact_num/2)+protein_num*max_contacts]);
			contact_energies += change_in_energy;
		}
		else if (new_contact_1 == contact_one){
			change_in_energy = (miyazawa_jernigan[mutations[dif_contact_pointer_1+1]][sequence_int[contact_two]]-
												energies[(contact_num/2)+protein_num*max_contacts]);
			contact_energies += change_in_energy;

		}
		else if (new_contact_2 == contact_two){
			change_in_energy = (miyazawa_jernigan[sequence_int[contact_one]][mutations[dif_contact_pointer_2+1]]-
									energies[(contact_num/2)+protein_num*max_contacts]);
			contact_energies += change_in_energy;

		}
		
		//Advance contact num and reset new contact_2 as this might be 0 later on
		contact_num += 2;

		contact_one = contact_matrix[contact_num];
		contact_two = contact_matrix[contact_num+1];

		dif_contact_pointer_2 = -2;
		new_contact_2 = -1;
	} while (contact_one != 0 || contact_two != 0);

	return contact_energies;

}



/*Finds the mean and standard deviation of the distribution of sequences based on the
protein sequence given to the program */
double * find_distribution_stats ( int** protein_contact_maps, bool first, int* mutations)
{

	double mean = 0.0;
	double mean2 = 0.0;
	double contact_energy;

	//Find the mean of the diverse proteins based on the sequence
	for (int protein_num = 0; protein_num < num_diverse_prots; protein_num ++){


		if (first){
			contact_energy = find_G_NS(protein_contact_maps[protein_num], protein_num);
		}
		else{
			contact_energy = find_G_NS_not_first(mutations, protein_contact_maps[protein_num],
					protein_num);
		}
		mean += contact_energy;
		mean2 += (contact_energy*contact_energy);
	}

	mean /= num_diverse_prots;
	mean2 /= num_diverse_prots;
	double sigma = mean2 - (mean*mean);


	double * summary_stats = malloc(2*sizeof(double));
	summary_stats[0] = mean;
	summary_stats[1] = sigma;
	
	return summary_stats;
	
	
}





/* Converts a sequence of characters to integers to be fed into M-J matrix */
void find_integer_sequence(char* sequence){

	//Allocate memory for the sequence in integer form
	sequence_int = malloc (sequence_length * sizeof(int));
	
	for (int seq_pos = 0; seq_pos < sequence_length; seq_pos ++){
		sequence_int[seq_pos] = sequence[seq_pos] - 65; //Converts ASCII character to alphabet character starting at 0

	}
}


/* Function made entirely for testing purposes - prints out an array*/
void print_array (int* array_to_print){
	for (int i =0; i < 600; i++){
		printf("%d, ", array_to_print[i]);
	}
	printf("\n");
}



/* Finds positions that are different from the consensus sequence */
int* find_dif_poses(char* sequence){
	//Allocate memory for the sequence in integer form
	int* dif_seqs = calloc (2*sequence_length, sizeof(int));
	int int_pos_val;
	int dif_seqs_pos = 0;


	for (int seq_pos = 0; seq_pos < sequence_length; seq_pos ++){
		int_pos_val = sequence[seq_pos] - 65; //Converts ASCII character to alphabet character starting at 0
		//If the sequence is different, record position of difference and what the new value is
		if(sequence_int[seq_pos] != int_pos_val){
			dif_seqs[dif_seqs_pos] = seq_pos;
			dif_seqs[dif_seqs_pos+1] = int_pos_val;
			dif_seqs_pos += 2;
		}
	}
	
	if(dif_seqs_pos == 0){
		return NULL;
	}
	
	return dif_seqs;
}



/*Computes the probability of folding for a certain protein using the formula by Pollock et al. (2012) */
void find_prob_fold(bool first, int* mutations, int* main_contact_mat)
{
	float G_NS;
	float sigma;
	float G_bar;
	double change_G;
	double p_fold;
	
	//Converts amino acid sequence to integer form and calculates the distribution for sequence

	double* summary_stats = find_distribution_stats(diverse_contacts, first, mutations);

	G_bar = summary_stats[0];
	sigma = summary_stats[1];

	free(summary_stats);

	//Finds the contact energy of the contact matrix provided by the user
	if(first){
		G_NS = find_G_NS(main_contact_mat, num_diverse_prots);
	}
	else{
		G_NS = find_G_NS_not_first(mutations, main_contact_mat, num_diverse_prots);
	}

	change_G = G_NS - G_bar + (kT * log_N_U) + (sigma/(2.0* kT));

	double factor = exp(-change_G/kT);
	p_fold = (factor/(1+factor));
	
	if(first){
		main_p_fold = p_fold;
	} else if (second){
		printf("%f", p_fold);
		second = false;
	} else {
		printf(",%f", p_fold);
	}
}






/* Main function to compute fitness based on set of sequences and contact map given in arguments */
int main(int argc, char *argv[])
{

	// clock_t begin = clock();

	sequence_length = atoi(argv[1]); //length of genome given by the user
	num_diverse_prots = atoi(argv[2]);
	log_N_U = log(ceil(pow(3.4, sequence_length))); //From Pollock et al.


	max_contacts = atoi(argv[6]); //argv[6] is the maximum number of contacts that a protein has
	max_line_size = atoi(argv[7]); //argv[7] is the maximum size of a line of contacts
	
	energies = calloc(max_contacts*(num_diverse_prots + 1), sizeof(double));
	main_energies = malloc((num_diverse_prots + 1)* sizeof(double));

	diverse_contacts = read_dat(argv[3], num_diverse_prots); //argv[3] is the file location of diverse contacts maps
	
	int** main_contact_mat = read_dat(argv[4], 1);//argv[4] is the file location of the main contact map



	FILE* stream = fopen(argv[5], "r"); //argv[5] is the file containing all of the sequences
	
	char seq[sequence_length + 1]; // Read each sequence of length sequence_length
	int len;

	bool first = true;
	char* sequence;
	int* dif_seqs;
	
	
	// Recurse through each sequence and find the probability of folding ie. fitness of the protein
   while (fgets(seq, sequence_length + 1, stream))
   {

		// Skip over newlines
		len = strlen(seq);
		if( seq[len-1] == '\n')
			continue;

		sequence = strdup(seq);
		
		// Put in commas after probability of folding type to form array of probabilities
		if(first){
			find_integer_sequence(sequence);
		} else {
			dif_seqs = find_dif_poses(sequence);
		}
		
		// Here the sequence is the same as the consensus so we just print out main_p_fold
		if(dif_seqs == NULL & !first){
			if (second){
				printf("%f", main_p_fold);
				second = false;
			} else {
				printf(",%f", main_p_fold);
			}
		} else {
			find_prob_fold(first, dif_seqs, main_contact_mat[0]);
		}
		
		
		if(!first){
			free(dif_seqs);
		}
		
		first = false;
		free(sequence);
   } 
	
	//Free everything that is open
	fclose(stream);
	
	free(energies);
	free(main_energies);


	free(*diverse_contacts);
	free(diverse_contacts);
	free(*main_contact_mat);
	free(main_contact_mat); 

	free(sequence_int);
	


	// clock_t end = clock();
	// double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	// time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	// printf("\n %f \n", time_spent);
	

	return 0;
	
}



