#Program to take output of SLiMTree simulation and create a consensus sequence
#and polymorphic sites document for further analysis in RaxML

#Required packages:
#BioPython
#random


import argparse, random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo

class process_data:
    #Controls the running of each of the commands required to create the consensus sequence data
    def __init__(self, nuc_fasta_file):

        self.parse_nucleotide_data(nuc_fasta_file)
        self.find_consensus_sequences()
        self.output_data(nuc_fasta_file)



    #Parses through the nucleotide fasta file and creates a dictionary of sequences in each population
    def parse_nucleotide_data(self, nuc_fasta_file):
        self.nuc_dict = {}
        
        for record in SeqIO.parse(nuc_fasta_file, "fasta"):
            population = record.id.split("_")[0]

            if (population in self.nuc_dict.keys()):
                self.nuc_dict[population].append(SeqRecord(record.seq, id = record.id))
            else:
                self.nuc_dict[population]=[SeqRecord(record.seq, id = record.id)]



    #Takes in a dictionary of lists of sequences and finds the consensus sequence for the sequences
    def find_consensus_sequences(self):
        self.consensus_sequences = []
        self.scoring_matrices = {}
        for population in self.nuc_dict:
            sequences = self.nuc_dict[population]
            align = MultipleSeqAlignment(sequences)
            align_info = AlignInfo.SummaryInfo(align)

            consensus_sequence = align_info.dumb_consensus(threshold=0.25)
            mutable_consensus_seq = consensus_sequence.tomutable()
            consensus_scoring_matrix = align_info.pos_specific_score_matrix()

            for i in range(len(mutable_consensus_seq)):
                if (mutable_consensus_seq[i] == 'X'):
                    possible_nucleotides = ['A', 'C', 'G', 'T']

                    max_score = 0
                    max_nucs = []

                    for nuc in possible_nucleotides:
                        score = consensus_scoring_matrix[i][nuc]

                        if (score > max_score):
                            max_nucs = [nuc]
                            max_score = score
                        elif (score == max_score):
                            max_nucs.append(nuc)
                    mutable_consensus_seq[i] = max_nucs[random.randint(0,len(max_nucs)-1)]
            print("Population "+ population + " complete")
            self.consensus_sequences.append(SeqRecord(mutable_consensus_seq, id = population, description = "" ))
            self.scoring_matrices[population] = str(consensus_scoring_matrix)



    #Outputs the consensus sequences to a fasta file and the scoring matrices to a tsv file
    def output_data(self, file_name):
        file_name_fasta = file_name.split(".fasta")[0]+"_consensus.fasta"
        
        SeqIO.write(self.consensus_sequences, file_name_fasta, "fasta")

        general_file_name = file_name.split(".fasta")[0]+"_consensus_"
        
        for i in self.scoring_matrices:
            writer = open(general_file_name + i + ".tsv", "w")
            writer.write( self.scoring_matrices[i])
            writer.close()


        

            


def main():
    #Parse for the fasta data files from the input given by the user
    parser = argparse.ArgumentParser(description='A program that reads output of SLiMTree simulations' +
                'and processes the data.')
    parser.add_argument('-i','--nucleotide_fasta_file', nargs = 1, required = True, type = str,
            help = 'Nucleotide fasta file output by the SLiMTree simulation')



    arguments = parser.parse_args()

    #Post process the data using commands in the process data class
    data_processor = process_data(arguments.nucleotide_fasta_file[0])




if __name__ == '__main__':
    main()
        
