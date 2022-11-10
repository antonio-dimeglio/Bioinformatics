from typing import Set
from random import randint
from math import log2
from sys import argv
from termcolor import colored

def get_coordinates(path:str = 'model-based-cpg-islands-hg19.txt') -> Set:
    '''
    This function returns a set containing tuples of coordinates
    '''

    coordinates = set()

    try:
        with open(path, 'r') as fstream:
            for line in fstream:
                lines = line.strip().split()
                if lines[0] == 'chr22':
                    coordinates.add((int(lines[1]), int(lines[2])))

    except FileNotFoundError:
        print(f"Could not find file {path}.")

    return coordinates

def get_genome(path:str = 'chr22.fa'):
    '''
    This function returns the content of a FASTA file
    '''
    genome = ""

    try:
        with open(path, 'r') as fstream:
            for line in fstream:
                if line[0] != '>':
                    genome += line.strip().upper()  
    
    except FileNotFoundError:
        print(f"Could not fine file {path}")
         
    return genome

def compute_model(model_type:str = "inside"):
    alphabet = ['A', 'C', 'G', 'T']    
    probabilities = {
        'AA':0, 'AC':0, 'AG':0, 'AT':0,
        'CA':0, 'CC':0, 'CG':0, 'CT':0,
        'GA':0, 'GC':0, 'GG':0, 'GT':0,
        'TA':0, 'TC':0, 'TG':0, 'TT':0,
    }

    coordinates = get_coordinates()
    genome = get_genome()

    for c in coordinates:
        if model_type == 'inside':
            island = genome[c[0]+1:c[1]+1]

        elif model_type == 'outside':
            start = randint(0, len(genome)-1) + 1
            end = start + (c[1] - c[0]) + 1
            island = genome[start:end]
        
        for i in range(1, len(island)):
            if island[i] != 'N' and island[i-1] != 'N':
                probabilities[island[i-1]+island[i]] += 1


    for c in alphabet:
        total = 0
        
        for key in probabilities.keys():
            if key[0] == c:
                total += probabilities[key]

        for key in probabilities.keys():
            if key[0] == c:
                probabilities[key] = round(probabilities[key] / total, 3)
    
    return probabilities

def pretty_print(d:dict):
    alphabet = ['A', 'C', 'G', 'T']   
    for i in alphabet:
        for j in alphabet:
            dimer = i + j
            print(f"{dimer}: {d[dimer]}\t", end='')
        print()
    print('-'*70)

def compute_query(inside_model:dict, outside_model:dict, query:str):
    outside = inside = log2(0.25)
     

    for i in range(1, len(query)):
        dimer = query[i-1]+query[i]
        inside += log2(inside_model[dimer])
        outside += log2(outside_model[dimer])
    
    return inside-outside

def main():
    if len(argv) == 1:
        print(colored("[USAGE] markov.py <query>\n", "green"),
        colored("<query>: DNA query string\n"), "yellow")
    if len(argv) == 2:    
        inside_model = compute_model("inside")
        outside_model = compute_model("outside")

        print("\t\t\tInside model")
        pretty_print(inside_model)
        print("\t\t\tOutside model")
        pretty_print(outside_model)
        print(f'Likelihood of being a CpG island: {compute_query(inside_model, outside_model, argv[1])}')

    else:
        print(colored("Too many arguments!", "red"))
        
if __name__ == '__main__':
    main()