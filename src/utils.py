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


def get_even_reads_from_fasta(file, read_len=250, amt=6):
    seq = load_fasta(file)
    step_size = int(len(seq)/(2*read_len))
    reads = []
    reads.append(seq)
    i=0
    while i <= len(seq)-read_len:
        reads.append(seq[i:i+read_len])
        randomness = np.random.randint(-2,2)
        i += int(read_len / amt) + randomness
    reads.append(seq[len(seq)-read_len:len(seq)])
    return reads


def get_reads_from_fasta(file, read_len=200, coverage=5):
    seq = load_fasta(file)
    L = len(seq)
    print(L)
    read_amt = int(L/read_len * coverage)
    reads = []
    for _ in range(read_amt):
        start = np.random.randint(0,L-read_len)
        read = seq[start:start+read_len]
        reads.append(read)
    return reads


def write_reads(reads, filename):
    with open(filename, "w") as txt_file:
        for i in range(len(reads)-1):
            line = reads[i]
            #txt_file.write(">"+str(i)+"\n")
            txt_file.write("".join(line) + "\n")
        print(line)
        line = reads[len(reads)-1]
        txt_file.write("".join(line))


def parse_nucleotides(sequence, caps=True):
    new_seq = []
    if caps:
        map_to_vals = {"A": 1, "C": 2, "G": 3, "T":4}
    else:
        map_to_vals = {"a": 1, "c": 2, "g": 3, "t":4}
    for symbol in sequence:
        new_seq.append(map_to_vals[symbol])
    
    return new_seq

def parse_num_to_nuc(sequence, caps=True):
    new_seq = []
    if caps:
        map_to_vals = {1: "A", 2: "C", 3: "G", 4 : "T"}
    else: 
        map_to_vals = {1: "a", 2: "c", 3: "g", 4 : "t"}
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
    #print(len(sequence))
    for seq in seqs:
        print(len(seq))
    return seqs[0]

def load_fasta_seqs(filename):
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    seqs = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seqs.append(sequence)
    #print(len(sequence))
    return seqs


def write_fasta(seqs, filename, nums=False):
    with open(filename, "w") as file:
        for i in range(len(seqs)-1):
            if nums:
                line = parse_num_to_nuc(seqs[i], caps=True)
            else: 
                line = seqs[i]
            strline = "".join(line) + "\n"
            file.write(">" + str(i) + "\n" + strline)
        
        #print(line)
        file.write(">last\n")
        if nums:
            parse_num_to_nuc(seqs[len(seqs)-1], caps=True)
        else: 
            line = seqs[len(seqs)-1]
        file.write("".join(line))

