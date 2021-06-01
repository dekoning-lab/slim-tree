#Program to take file of trees and make a csv of the branch lengths 

#Required packages:
#Pandas
#Argparse


import pandas as pd
import argparse

class find_branch_lengths:
    #Controls the running of each of the commands required to make the csv file
    def __init__(self, tree_file):
        self.make_branch_dictionary(tree_file)


    #Creates a dictionary of branch lengths from a tree file
    def make_branch_dictionary(self, tree_file):
        line = tree_file.readline()
        branch_length_dict = {}
        populations_dict = {'p6':'p4', 'p7': 'p3', 'p8': 'p2'}
        count = 1;
        while (line != ''):
            #Separate each tree into a list of strings containing each of the branches
            line = line.replace('(', ',')
            line = line.replace(')', ',')
            line = line.replace(';\n', ',')
            branches = line.split(',')
            branches = list(filter(None, branches))

            previous_population = ''
            for leaf in branches:
                leaf = leaf.split(':')
                population = leaf[0]
                branch_length = leaf[1]
                if (population == ''):
                    population = populations_dict[previous_population]

                if(population in branch_length_dict):
                    branch_length_dict[population].append(branch_length)
                else:
                    branch_length_dict[population] = [branch_length]

                previous_population = population


            # for populations in branch_length_dict:
                # if(len(branch_length_dict[populations]) < count):
                    # branch_length_dict[populations].append("NA")

            count += 1

            line = tree_file.readline()
            
                

            
        pd.DataFrame(branch_length_dict).to_csv('branch_lengths.csv', header = True)
        tree_file.close()
        

            


def main():
    #Parse for the data file of trees given by the user
    parser = argparse.ArgumentParser(description='A program which takes input of trees file and' +
                                     'constructs a csv file of branch lengths')
    parser.add_argument('-i','--trees_file', nargs = 1, required = True, type = argparse.FileType('r'),
            help = 'File containing trees to find branch lengths from')

    

    arguments = parser.parse_args()

    #Find branch lengths from the data file
    data_processor = find_branch_lengths(arguments.trees_file[0])




if __name__ == '__main__':
    main()
        
