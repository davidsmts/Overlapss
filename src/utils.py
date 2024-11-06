import numpy as np
from Bio import SeqIO

def get_precision_recall(estimated_positions, true_positions):
    TP, FP, FN = 0, 0, 0
    for pos in estimated_positions:
        if pos in true_positions:
            TP += 1
        else:
            FP += 1

    FN = len(true_positions) - len(estimated_positions)
    if FN < 0: FN = 0
    if TP + FN == 0: R = 0
    else: R = TP / (TP + FN)

    if TP + FP == 0: P = 0
    else: P = TP / (TP + FP)
    return P, R


def get_sens_spec(estimated_positions, true_positions, all_pos):
    sens = 0
    for pos in true_positions:
        if pos in estimated_positions:
            sens += 1
    sens = sens / len(true_positions)

    # specificity
    FP = 0
    for pos in estimated_positions:
        if pos not in true_positions:
            FP += 1
    spec = (len(all_pos) - len(true_positions) - FP) / (len(all_pos) - len(true_positions))

    fp_rate = FP / (len(all_pos) - len(true_positions))
    
    return sens, spec, fp_rate


def get_reads_from(seq, read_len=200, coverage=5):
    L = len(seq)
    read_amt = int(L/read_len * coverage)
    reads = []
    for _ in range(read_amt):
        start = np.random.randint(0,L-read_len)
        read = seq[start:start+read_len]
        reads.append(read)
    return reads


def parse_nucleotides(sequence):
    new_seq = []
    map_to_vals = {"a": 1, "c": 2, "g": 3, "t":4}
    for symbol in sequence:
        new_seq.append(map_to_vals[symbol])
    
    return new_seq


def getskmer(snippet, profile):
    spaced_kmer = snippet * profile
    spaced_kmer = spaced_kmer[spaced_kmer != 0]
    s = ''.join(str(x) for x in spaced_kmer)
    return s


def load_fasta(filename):
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    seqs = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seqs.append(sequence)
    return seqs