import numpy as np
import matplotlib.pyplot as plt
import utils

def get_maxcount(pos, seqs, seqs_kmers, spaced_kmer_profile, seq_to_investigate=0):
    target = seqs[seq_to_investigate]
    f = len(spaced_kmer_profile)
    counts_i = []
    indexes = []
    len_loop = min(f-1, pos)
    start = pos-len_loop
    end = pos+1
    for i in range(start, end):
        if spaced_kmer_profile[pos-i] != 1:
            counts_i.append(0)
            indexes.append(-1)
            continue
        
        # Extract k mers starting at selected position
        spaced_kmer = target[i : i+f] * spaced_kmer_profile
        spaced_kmer = spaced_kmer[spaced_kmer != 0]
        s = ''.join(str(x) for x in spaced_kmer)
        counts_i.append(seqs_kmers[s])
        indexes.append(i)
    
    return max(counts_i), indexes[np.argmax(counts_i)]



#
# CORRECT COUNTS - AN ATTEMPT AT FIXING THE ECLIPSE ERROR
#
def correct_counts(maxed_counts, maxed_count_indices, target_sequence, start_end_posis, seqs, diff_profile, kmer_profile):
    f = len(kmer_profile)
    correction_artifact = []
    last_correction = 0 
    # When corrections are necessary they are written into the maxed_counts array
    for i in range(1,len(maxed_count_indices)):
        # Make sure we have a change position on our hands.
        if maxed_counts[i] != maxed_counts[i-1]: 
            # Correct for reads that start in between the two maxed_count indices:
            start_sum_at = last_correction
            end_sum_at = maxed_count_indices[i]+1
            #end_sum_at = i+1
            
            if start_sum_at <= end_sum_at:
                correction = sum(start_end_posis[start_sum_at:end_sum_at])
            else:
                correction = -sum(start_end_posis[end_sum_at:start_sum_at])
            
            last_correction = maxed_count_indices[i] + 1
            diff_profile[i-1] += correction
            correction_artifact.append([correction, (start_sum_at, end_sum_at), start_end_posis[start_sum_at:end_sum_at]])

    return diff_profile, correction_artifact



def get_correction_profile(target, seqs, overlap_size, start_dict, end_dict):
    corr_profile = [0 for i in range(len(target)-overlap_size+1)]
    #corr_profile[0] = -1
    for i in range(len(target)-overlap_size):
        start_snippet = ''.join(str(x) for x in target[i : i+overlap_size])      
        if start_snippet in start_dict:

            corr_profile[i] -= start_dict[start_snippet]

    for i in range(overlap_size, len(target)+1):
        end_snippet = ''.join(str(x) for x in target[i-overlap_size : i])
        if end_snippet in end_dict:
            corr_profile[i-overlap_size+1] += end_dict[end_snippet]
    
        
    corr_profile[0]=0
    return corr_profile



def correct_diff_profile(filename, str_profile, seq_to_investigate, data=[]):
    seqs = []
    if data:
        for read in data:
            sequence_chars = [val for val in read]
            sequence = utils.parse_nucleotides(sequence_chars)
            seqs.append(np.array(sequence))
    else:
        with open(filename) as file_in:
            for line in file_in:
                newline = line.rstrip('\n')
                sequence_chars = [char for char in newline]
                sequence = utils.parse_nucleotides(sequence_chars)
                seqs.append(np.array(sequence))
    
    profile = [int(character) for character in str_profile]
    k = sum(profile)
    f = len(profile)
    
    # Turn into np arrays for componentwise multiplication
    profile = np.array(profile)
    
    # Count occurence of spaced k-mers
    starts, ends = {}, {}
    seqs_kmers = {}
    for sequence in seqs:
        for i in range(len(sequence) - f):
            spaced_kmer = sequence[i:i+f] * profile
            spaced_kmer = spaced_kmer[spaced_kmer != 0]
            s = ''.join(str(x) for x in spaced_kmer)
            seqs_kmers[s] = seqs_kmers.get(s, 0) + 1

            if i == 0 or i == len(sequence)-f-1:
                solid_kmer = sequence[i:i+f]
                s2 = ''.join(str(x) for x in solid_kmer)
                addon = 0
                if i == 0:
                    starts[s2] = starts.get(s2, 0) + 1
                else: 
                    ends[s2] = ends.get(s2, 0) + 1

    
    # Get maxcounts from counts
    target = seqs[seq_to_investigate]
    xpoints = np.array([i for i in range(len(target) - f)])
    max_counts = []
    max_count_indices = []
    for i in range(len(xpoints)):
        maxp, argmaxp = get_maxcount(i, seqs, seqs_kmers, profile, seq_to_investigate=seq_to_investigate)
        max_counts.append(maxp)
        max_count_indices.append(argmaxp)

    # Get correction profile:
    start_end_posis = get_correction_profile(target, seqs, f, starts, ends)
    
    # Get diff profile:
    pre_corr_diff_profile = [max_counts[j] - max_counts[j-1] for j in range(1,len(max_counts))]
    
    # Apply correction strategy
    ypoints, correction_artifact = correct_counts(max_counts, max_count_indices, target, start_end_posis, seqs, pre_corr_diff_profile.copy(), profile)
    return xpoints[1:], ypoints, max_counts, start_end_posis, pre_corr_diff_profile, correction_artifact, max_count_indices
