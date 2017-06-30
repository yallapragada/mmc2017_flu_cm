from Bio import AlignIO
from numpy import array, transpose
from scipy.stats import mode
import csv,sys
import operator


#usage: python create_01_file.py [name of aligned fasta file]
#converts an aligned fasta file to 0s and 1s
#0 if residue = consensus residue, 1 otherwise
#writes output 0,1 as input_file.csv


def get_best_mutation(residues, site):
    residue_count = {}
    for residue in residues:
        if residue_count.get(residue) is None:
            residue_count[residue] = 1
        else:
            residue_count[residue] += 1
    sorted_residue_count = sorted(residue_count.items(), key=operator.itemgetter(1), reverse=True)

    if len(sorted_residue_count) > 1:
        consensus_count = int(sorted_residue_count[0][1])
        top_mutated_count = int(sorted_residue_count[1][1])

        if top_mutated_count > consensus_count/10:
            return sorted_residue_count[0][0] + site + sorted_residue_count[1][0]
        else:
            return sorted_residue_count[0][0] + site + sorted_residue_count[0][0]
    else:
        return sorted_residue_count[0][0] + site + sorted_residue_count[0][0]



def write_list_to_file(mutations, p_sequences):
    with open('test.csv', 'w', encoding='utf8', newline='') as test_file:
        writer = csv.writer(test_file)
        writer.writerow(mutations)
        writer.writerows(p_sequences)
    test_file.close()


def perform_transpose(sequences):
    transposed_list = transpose(array(sequences)).tolist()
    return transposed_list


def modify_sequences(sequences):
    sequencesT = perform_transpose(sequences=sequences)
    new_sequencesT = [check_gaps(sequenceT) for sequenceT in sequencesT]
    mutations = [get_best_mutation(new_sequenceT, str(counter+1)) for counter, new_sequenceT in enumerate(new_sequencesT)]
    new_sequences = perform_transpose(sequences=new_sequencesT)
    return new_sequences, mutations


# discard sequence if more than 20% gaps
# we do this by converting the complete sequence to a gap and hence it will imply static residue across strains
def check_gaps(sequence):
    num_gaps = str(sequence).count('-')
    len_sequence = len(sequence)
    if num_gaps > len(sequence)/5:
        sequence = list('-' * len_sequence)
    return sequence


# input = list of sequence strings (raw sequences)
def get_consensus_sequence_new(sequences):
    mod = [mode(y)[0][0] for y in transpose(array([list(z) for z in sequences]))]
    consensus_sequence = ''.join(mod)
    print(consensus_sequence)
    return consensus_sequence


def write_list_to_file(mutations, p_sequences, p):
    with open(p+'_01.csv', 'w', encoding='utf8', newline='') as test_file:
        writer = csv.writer(test_file)
        writer.writerow(mutations)
        writer.writerows(p_sequences)
    test_file.close()


def create_01_sequences_gaps(file):
    p = file.split('.')[0]
    sequences = AlignIO.read(file, 'fasta')

    p_sequences = []  # list of p1 sequences

    for sequence in sequences:
        p_sequences.append(list(sequence.seq))

    p_new_sequences, mutations = modify_sequences(p_sequences)
    consensus_sequence = get_consensus_sequence_new(p_new_sequences)
    p_sequences_01 = [convert_01(p_new_sequence, consensus_sequence) for p_new_sequence in p_new_sequences]

    write_list_to_file(mutations, p_sequences_01, p)


def convert_01(sequence, consensus_sequence):
    return [0 if x==consensus_sequence[i] else 1 for i,x in enumerate(sequence)]


def run(file):
    create_01_sequences_gaps(file)


if __name__ == '__main__':
    file = sys.argv[1]
    run(file=file)