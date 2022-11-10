import numpy as np
import sys


#CONSTANTS

#COLOR CONSTANTS
OK = '\033[92m' #GREEN
WARNING = '\033[93m' #YELLOW
FAIL = '\033[91m' #RED
RESET = '\033[0m' #RESET COLOR

#ALIGNMENT CONSTANTS


def global_alignment(first_sequence:str, second_sequence:str, match = 1, mismatch = -1, gap = -2):
    align = lambda c1, c2: match if c1 == c2 else mismatch

    score_matrix = np.zeros((len(first_sequence)+1, len(second_sequence)+1), dtype=np.int32)
    
    for i in range(1, len(second_sequence) + 1):
        score_matrix[0, i] = i*gap

    for i in range(1, len(first_sequence)+1):
        score_matrix[i, 0] = i*gap

    for i in range(1, len(first_sequence) + 1):
        for j in range(1, len(second_sequence) + 1):
            score_matrix[i, j] = max(
                score_matrix[i-1, j-1] + align(first_sequence[i-1], second_sequence[j-1]),
                score_matrix[i-1, j] + gap,
                score_matrix[i, j-1] + gap
            )

    first_aligned = ""
    second_aligned = ""
    i = len(first_sequence)
    j = len(second_sequence)

    while(i != 0 or j != 0):
        if i == 0:
            first_aligned = "-" + first_aligned
            second_aligned = second_aligned[j] + second_aligned
            j-=1 
        elif j == 0:
            second_aligned = "-" + second_aligned
            first_aligned = first_aligned[j] + first_aligned
            i-=1
        else:
            l = [
                score_matrix[i-1, j-1],
                score_matrix[i-1, j],
                score_matrix[i, j-1]
            ]
            local_max = l.index(max(l))
            if first_sequence[i-1] == second_sequence[j-1] or local_max == 0:
                first_aligned = first_sequence[i-1] + first_aligned
                second_aligned = second_sequence[j-1] + second_aligned
                i -= 1
                j -= 1
            elif local_max == 1:
                second_aligned = "-" + second_aligned
                first_aligned = first_sequence[i-1] + first_aligned
                i -= 1
            else:
                first_aligned = "-" + first_aligned
                second_aligned = second_sequence[j-1] + second_aligned
                j -= 1
    return [first_aligned, second_aligned]


def local_alignment(first_sequence:str, second_sequence:str, match=1, mismatch=-1, gap=-2):
    to_zero = lambda x: 0 if x < 0 else x
    align = lambda c1, c2: match if c1 == c2 else mismatch

    score_matrix = np.zeros((len(first_sequence)+1, len(second_sequence)+1), dtype=np.int32)

    for i in range(1, len(first_sequence)+1):
        for j in range(1, len(second_sequence)+1):
            score_matrix[i, j] = max(
                to_zero(score_matrix[i-1, j-1] + align(first_sequence[i-1], second_sequence[j-1])),
                to_zero(score_matrix[i-1, j] + gap),
                to_zero(score_matrix[i, j-1] + gap)
            )


    i, j = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
    first_aligned = first_sequence[i-1]
    second_aligned = second_sequence[j-1]

    while(score_matrix[i, j] != 0):
        if i == 0:
            first_aligned = "-" + first_aligned
            second_aligned = second_aligned[j] + second_aligned
            j-=1 
        elif j == 0:
            second_aligned = "-" + second_aligned
            first_aligned = first_aligned[j] + first_aligned
            i-=1
        else:
            l = [
                score_matrix[i-1, j-1],
                score_matrix[i-1, j],
                score_matrix[i, j-1]
            ]
            local_max = l.index(max(l))
            if first_sequence[i-1] == second_sequence[j-1] or local_max == 0:
                first_aligned = first_sequence[i-1] + first_aligned
                second_aligned = second_sequence[j-1] + second_aligned
                i -= 1
                j -= 1
            elif local_max == 1:
                second_aligned = "-" + second_aligned
                first_aligned = first_sequence[i-1] + first_aligned
                i -= 1
            else:
                first_aligned = "-" + first_aligned
                second_aligned = second_sequence[j-1] + second_aligned
                j -= 1
    
    return [first_aligned, second_aligned]


def parse_fasta(path:str) -> str:
    seq = ""
    with open(path, 'r') as f:
        for line in f:
            if line[0] != '>':
                seq += line.strip()
    return seq
    
def askparams():
    params = [1, -1, -2]
    yn = input(
        "Would you like to set the default value for the score function to something different?\n"
        "They are currently set as\n" + WARNING + 
        f"MATCH: {params[0]}\nMISMATCH: {params[1]}\nGAP: {params[2]}\n" + RESET + 
        "Type y if you want to change them, otherwise just press enter "
    )
    if yn in ['y', 'Y']:
        params[0] = int(input("Match value: "))
        params[1] = int(input("Misatch value: "))
        params[2] = int(input("Gap value: "))
    return params

def print_usage():
    print(
        OK +"[USAGE]: alignment.py <flag> <alignment_type> <sequence> <sequence>\n" +
        RESET + "<flag>: -l to load sequence from fasta files, this requires two paths to be given in input\n"
        "<flag>: -i to input the sequences directly from the command line\n"
        "<alignment_type>: either " + FAIL + "global " + RESET + "or "+ FAIL + "local\n" + RESET +
        "<sequence>: either a raw sequence or the path to the file\n" +
        WARNING +"WARNING: when trying to load a sequence the script expects the files to be in the cwd (i.e the folder from which you launch the script)"
    )

def main():
    if len(sys.argv) != 5: print_usage()
    else:
        flag = sys.argv[1]
        alignment_type = sys.argv[2]
        if (flag not in ['-i', '-l']):
            print(FAIL + "Wrong flag in input!")
            exit(-1)
        if (alignment_type not in ["local", "global"]):
            print(FAIL + "Wrong alignment type in input!")
            exit(-1)
        if flag == '-l':
            first_sequence = parse_fasta(sys.argv[3])
            second_sequence = parse_fasta(sys.argv[4])
        if flag == '-i':
            first_sequence = sys.argv[3]
            second_sequence = sys.argv[4]
        match, mismatch, gap = askparams()
        if alignment_type == "global":
            first_aligned, second_aligned = global_alignment(first_sequence, second_sequence, match, mismatch, gap)
        if alignment_type == "local":
            first_aligned, second_aligned = local_alignment(first_sequence, second_sequence, match, mismatch, gap)
        print(f"Results:\n{first_aligned}\n{second_aligned}")
        
    

if __name__ == '__main__':
    main()