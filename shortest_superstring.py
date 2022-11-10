from typing import List
from itertools import permutations
import random
import numpy as np
from termcolor import colored
from sys import argv


#Helper functions
def get_genome(length):
    genome = "".join(random.choice("AGCT") for i in range(length))
    return genome

def get_substrings(seq, length):
    s = []
    for i in range(len(seq) - length+1):
        s.append(seq[i:i+length])
    return s

def find_overlap(suffix:str, prefix:str) -> int:
    ol = 0 
    for i in range(len(suffix)):
        if suffix[-i:] == prefix[:i]: ol = i
    return ol

def traverse_overlap_graph(reads:List[str]):
    #To build the overlap graph we create a matrix
    #For each read we'll have a row and a column 
    #The overlap between each read is computed 
    #(with the exception being that same read)
    overlap_matrix = np.zeros(shape=(len(reads), len(reads)), dtype=np.int32)
    for i in range(len(reads)):
        for j in range(len(reads)):
            overlap_matrix[i, j] = -1 if i == j else find_overlap(reads[i], reads[j])

    
    i, j = np.unravel_index(overlap_matrix.argmax(), overlap_matrix.shape)
    return (i, j, overlap_matrix[i, j])


def greedy_scs(reads:List[str]) -> str:

    while (len(reads) > 1):
        i, j, overlap = traverse_overlap_graph(reads)
        reads[i]+= reads[j][overlap:]
        reads.pop(j) 

    return reads[0]


def bruteforce_scs(reads:List[str]) -> str:
    perms = permutations(reads)
    best_length = float('inf')
    best = ""
    t = 0
    for perm in perms:
        t +=1 
        current = greedy_scs(list(perm))
        if len(current) < best_length:
            best_length = len(current)
            best = current
    return best

def main():
    if len(argv) == 1 or len(argv) != 4:
        print(colored('[USAGE]: shortest_superstring.py <flag> <genome_length> <kmer_length>\n', 'white'),\
            colored('<flag>: -bf for bruteforce -gr for greedy\n', 'yellow'), \
            colored('<genome_length>: length of the genome to use for benchmarking\n', 'yellow'), \
            colored('<kmer_length>: length to use for kmerization \n', 'yellow'), \
            colored("NOTE: Avoid using kmer_length <= 3", 'red'),
            sep=''
        )
    else:
        genome = get_genome(int(argv[2]))
        kmers = get_substrings(genome, int(argv[3]))

        if argv[1] == '-bf':
            result = bruteforce_scs(kmers)
        else:
            result = greedy_scs(kmers)
        
        print(colored("Initial genome:\n", 'red'), f"{genome}\n", colored("Results:\n", 'green'), f"{result}", sep='')



if __name__ == '__main__':
    main()
