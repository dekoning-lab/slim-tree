#Program to take in a DNA model name and mutation rate and return a matrix of mutation
#rates from each mutation to each mutation to be input into SLiM

#Required packages:
#string
#sys

import string, sys

class getDNAModel:

    #Main method to set up mutation model 
    def __init__(self, model_name, mutation_rate):
        self.model_name = model_name.translate(str.maketrans('', '', string.punctuation)).lower()
        self.mutation_rate = mutation_rate
        
        self.make_model_hash()
        

    
    #Creates a dictionary of acceptable mutation models and links to their respective functions
    def make_model_hash(self):
        self.mutation_model_dict = {
            "jukescantor" : self.jukes_cantor(),
            "jc" : self.jukes_cantor()
        }
    
    
    #Creates a matrix of mutation rates according to Kimura's (1980) model
    def jukes_cantor(self):
        model_mutation_rate = self.mutation_rate/3
        mutation_matrix = ([0] + [model_mutation_rate]* 3 + 
                [model_mutation_rate] + [0] + [model_mutation_rate] * 2 +
                [model_mutation_rate] * 2 + [0] + [model_mutation_rate] +
                [model_mutation_rate]* 3 + [0])
        return str(mutation_matrix)[1:-1]
        
    #Creates a matrix of mutation rates according to Jukes and Cantor's (1969) model    
    def kimura_two_parameter
                
               
    #Returns the model given by the user    
    def get_model(self):
        if(self.model_name in self.mutation_model_dict.keys()):
            return self.mutation_model_dict[self.model_name]
        else:
            print("At this time we only accept " +
                "Jukes-Cantor (JC), Kimura-Two-Parameter (K2P) " +
                "models of DNA evolution. If your desired DNA model is not one of these models, " +
                "please make an issue request or, even better, add the model yourself to the 'getDNAModel.py' " +
                "script and send us a pull request! If your desired DNA model is one of these models, make sure" + 
                " you are writing the model as specified above")
            sys.exit(0)
               
               
if (__name__ == '__main__'):
    mutation_model = input("Enter mutation model")
    mutation_rate = float(input("Enter mutation rate"))
    transversion_rate = float(input("Enter transversion rate"))
    model = getDNAModel(mutation_model, mutation_rate)
    print(model.get_model())