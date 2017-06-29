import os
from Bio import AlignIO
from numpy import transpose, array
from scipy.stats import mode
import sys
from minepy import MINE
import csv

def perform_transpose(sequences):
    transposed_list = transpose(array(sequences)).tolist()
    return transposed_list


# Given a SeqRecord, return strain_name
def get_strain_name(record):
    strain_full_str = record.description.split('|')[3]
    if (strain_full_str.startswith('Strain')):
        strain = strain_full_str.split(':')[1]
    else:
        strain = record.description.split('|')[4].split(':')[1]
    return strain

def get_matching_sequence(records, strain_name):
    for record in records:
        if (get_strain_name(record) == strain_name):
            return record
    return None


def modify_sequences(sequences):
    sequencesT = perform_transpose(sequences=sequences)
    new_sequencesT = [check_gaps(sequenceT) for sequenceT in sequencesT]
    new_sequences = perform_transpose(sequences=new_sequencesT)
    return new_sequences


# discard sequence if more than 10% gaps
# we do this by converting the complete sequence to a gap and hence it will imply static residue across strains
def check_gaps(sequence):
    num_gaps = str(sequence).count('-')
    len_sequence = len(sequence)
    if num_gaps > len(sequence)/10:
        sequence = list('-' * len_sequence)
    return sequence


# input = list of sequence strings (raw sequences)
def get_consensus_sequence_new(sequences):
    mod = [mode(y)[0][0] for y in transpose(array([list(z) for z in sequences]))]
    consensus_sequence = ''.join(mod)
    return consensus_sequence


def create_01_sequences_gaps(file1, file2):
    sequences1 = AlignIO.read(file1, 'fasta')
    sequences2 = AlignIO.read(file2, 'fasta')

    p1_sequences = []  # list of p1 sequences
    p2_sequences = []  # list of p2 sequences

    for sequence1 in sequences1:
        strain_name = get_strain_name(sequence1)
        sequence2 = get_matching_sequence(sequences2, strain_name=strain_name)
        if (sequence2):
            p1_sequences.append(list(sequence1.seq))
            p2_sequences.append(list(sequence2.seq))

    p1_new_sequences = modify_sequences(p1_sequences)
    p2_new_sequences = modify_sequences(p2_sequences)

    consensus_sequence1 = get_consensus_sequence_new(p1_new_sequences)
    consensus_sequence2 = get_consensus_sequence_new(p2_new_sequences)

    p1_sequences_01 = [convert_01(p1_new_sequence, consensus_sequence1) for p1_new_sequence in p1_new_sequences]
    p2_sequences_01 = [convert_01(p2_new_sequence, consensus_sequence2) for p2_new_sequence in p2_new_sequences]

    return p1_sequences_01, p2_sequences_01


def convert_01(sequence, consensus_sequence):
    return [1 if x==consensus_sequence[i] else 0 for i,x in enumerate(sequence)]

def write_mics_to_csv(mics, p1, p2, cutoff):
    keys = ['p1','p2','x','y','weight']
    mics_file_name = p1 + '_' + p2 + '_' + cutoff + '.csv'
    with open(mics_file_name,'w',encoding='utf8', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(mics)
    output_file.close()

def perform_mic_2p(p1_sequences, p2_sequences, p1, p2, cutoff=0.5):
    mic_scores = []
    p1_sequences_t = transpose(array([list(z) for z in p1_sequences])).tolist()
    p2_sequences_t = transpose(array([list(z) for z in p2_sequences])).tolist()

    for idx1, record1 in enumerate(p1_sequences_t):
        for idx2, record2 in enumerate(p2_sequences_t):
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(record1, record2)
            if (mine.mic() > float(cutoff)):
                mic_score = {}
                mic_score['x'] = p1+'_'+str(idx1+1)
                mic_score['y'] = p2+'_'+str(idx2+1)
                mic_score['p1'] = p1
                mic_score['p2'] = p2
                mic_score['weight'] = format(mine.mic(), '.3f')
                mic_scores.append(mic_score)

    write_mics_to_csv(mics=mic_scores, p1=p1, p2=p2, cutoff=cutoff)
    return mic_scores


def run(file1, file2, cutoff):
    p1 = file1.split('.')[0]
    p2 = file2.split('.')[0]

    p1_sequences_01, p2_sequences_01 = create_01_sequences_gaps(file1, file2)
    mic_scores = perform_mic_2p(p1_sequences_01, p2_sequences_01, p1, p2, cutoff=cutoff)

    print('done with run_01 for ', file1, file2, cutoff)


if __name__ == '__main__':
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    cutoff = sys.argv[3]
    run(file1=file1, file2=file2, cutoff=cutoff)

