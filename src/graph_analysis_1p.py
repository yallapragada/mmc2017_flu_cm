import sys
from Bio import AlignIO
from minepy import MINE
from numpy import array, transpose
from scipy.stats import mode
import csv


def write_mics_to_csv(mics, p, cutoff):
    keys = ['x','y','weight']
    mics_file_name = p + '_' + cutoff + '.csv'
    with open(mics_file_name,'w',encoding='utf8', newline='') as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(mics)
    output_file.close()


def performMIC(transposed_list, cutoff):
    mic_scores=[]
    for counter1 in range(0, len(transposed_list)-1):
        for counter2 in range(counter1+1, len(transposed_list)):
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(transposed_list[counter1], transposed_list[counter2])
            if (mine.mic() > float(cutoff)):
                mic_score={}
                mic_score['x']=counter1
                mic_score['y']=counter2
                mic_score['weight']=mine.mic()
                mic_scores.append(mic_score)
    return mic_scores


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


def create_01_sequences_gaps(file):
    sequences = AlignIO.read(file, 'fasta')

    p_sequences = []  # list of p1 sequences

    for sequence in sequences:
        p_sequences.append(list(sequence.seq))

    p_new_sequences = modify_sequences(p_sequences)
    consensus_sequence = get_consensus_sequence_new(p_new_sequences)
    p_sequences_01 = [convert_01(p_new_sequence, consensus_sequence) for p_new_sequence in p_new_sequences]

    return p_sequences_01


def convert_01(sequence, consensus_sequence):
    return [1 if x==consensus_sequence[i] else 0 for i,x in enumerate(sequence)]


def run(file, cutoff):
    p = file.split('.')[0]

    p_sequences_01 = create_01_sequences_gaps(file)
    p_sequences_01_T = perform_transpose(p_sequences_01)
    mic_scores = performMIC(p_sequences_01_T, cutoff=cutoff)
    write_mics_to_csv(mics=mic_scores, p=p, cutoff=cutoff)

    print('done with run_01 for ', file, cutoff)


if __name__ == '__main__':
    file = sys.argv[1]
    cutoff = sys.argv[2]
    run(file=file, cutoff=cutoff)
