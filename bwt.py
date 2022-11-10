from sys import argv
from termcolor import colored

def bwt(T:str):
    #All rotations are generated
    T += '$'
    rotations = []

    for i in range(len(T)):
        rotations.append(T)
        T = T[-1] + T[0:len(T) - 1]
    

    rotations = sorted(rotations)
    
    return ''.join(map(lambda x : x[-1], rotations))

def count_others(values, count, c):
    total = 0 
    for i in range(len(values)):
        if values[i] == c:
            return total
        else:
            total += count[i]

def count_occurances(T):
    s = sorted(list(set(T)))
    l = [0] * len(s)

    for c in T:
        l[s.index(c)] += 1 
    return (s, l)

def reverse_bwt(T:str):

    starting_str = '$'
    values, counts = count_occurances(T)

    i = 0 

    while T[i] != '$':
        starting_str = T[i] + starting_str
        i = count_others(values, counts, T[i]) + (counts[values.index(T[i])] - T[i:].count(T[i]))
        

    return starting_str

def match_input(T:str, S:str):

    values, counts = count_occurances(T)

    possible_hits = [] 
    offset = 1

    #we first find all position in which the last character in S is found in T

    for i in range(len(T)):
        if T[i] == S[-1]:
            possible_hits.append(i)

    #then we check all possible hits, given that any are found
    #if an hit is not valid we remove it from the list
    while offset <= len(S) and possible_hits:
        i = 0
        while i < len(possible_hits):
            j = 0
            current_str = ""
            current_pos = possible_hits[i]
            while j < offset:
                current_str = T[current_pos] + current_str
                current_pos = count_others(values, counts, T[current_pos]) + counts[values.index(T[current_pos])] - T[current_pos:].count(T[current_pos])
                j+=1
            
            if current_str != S[-offset:]:
                possible_hits.remove(possible_hits[i])
                if i > 0: i -= 1 
            else: 
                i += 1
        
        offset += 1

    for j in range(len(possible_hits)):
        distance = 0
        i = possible_hits[j]
        while T[i] != '$':
            i = count_others(values, counts, T[i]) + (counts[values.index(T[i])] - T[i:].count(T[i]))
            distance += 1
        possible_hits[j] = distance - len(S)
    
    return possible_hits



def main():
    if len(argv) == 1:
        print(
            colored("[USAGE] btw.py <reverse> <input> <match> <input> \n", "green"),
            colored("<input>: input string (must be inline)\n", "yellow"),
            colored("<reverse>: -r to reverse a string already compressed, if this flag is not found the input string will be compressed\n", "yellow"),
            colored("<match>: -m to match the string put after this flag to the input \n", "yellow"),
            colored("WARNING: when reversing a compressed string remember to put a \\ before the dollar sign (tested on ubuntu only)\n Otherwise the $ sign wont be picked up as it is a reserved character in unix systems\n", "red"),
            colored("WARNING: when matching for a string the script assumes that the input string has not yet been compressed\n", "red")
        )
    elif len(argv) == 2:
        print(bwt(argv[1]))
    elif len(argv) == 3 and "-r" in argv:
        print(reverse_bwt(argv[2]))
    elif len(argv) == 4 and "-m" in argv:
        results = match_input(argv[1], argv[3])
        print(f"{len(results)} matches have been found for the query in the input string")
        if len(results) > 0:
            results.sort()
            print(f"Match found at indexes:{results}")
    else:
        print(colored("Something went wrong while parsing the argument passed to the script.","red"))
        print(argv)

if __name__ == "__main__":
    main()
    