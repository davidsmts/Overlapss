import numpy as np
import matplotlib.pyplot as plt
import utils

def get_maxcount(pos, seqs_kmers, spaced_kmer_profile, target):
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
    corr_profile = [0 for i in range(len(target)-overlap_size+2)]
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
        maxp, argmaxp = get_maxcount(i, seqs_kmers, profile, target)
        max_counts.append(maxp)
        max_count_indices.append(argmaxp)

    # Get correction profile:
    start_end_posis = get_correction_profile(target, seqs, f, starts, ends)
    
    # Get diff profile:
    pre_corr_diff_profile = [max_counts[j] - max_counts[j-1] for j in range(1,len(max_counts))]
    
    # Apply correction strategy
    ypoints, correction_artifact = correct_counts(max_counts, max_count_indices, target, start_end_posis, seqs, pre_corr_diff_profile.copy(), profile)
    return xpoints[1:], ypoints, max_counts, start_end_posis, pre_corr_diff_profile, correction_artifact, max_count_indices


def detect_repeats_greedy(ypoints):
    spikes = [pos for pos in range(len(ypoints)) if ypoints[pos] > 0]
    dips = [pos for pos in range(len(ypoints)) if ypoints[pos] < 0]
    has_spike = []
    repeat_regions = []
    # Map dips to spikes
    for pos in range(len(ypoints)):
        if pos in dips:
            if len(has_spike) > 0:
                furthest_spike = has_spike[0]
                has_spike.pop(0)
                repeat_regions.append((furthest_spike, pos))

        if pos in spikes:
            has_spike.append(pos)

    # Detect repeat at the beginning if the first irregularity is a dip
    for i in range(len(dips)):
        if dips[i] <= spikes[0]:
            repeat_regions.append((0,dips[i]))

    # Detect repeat at the end if the last irregularity is a spike
    for i in range(len(spikes)):
        if dips[len(dips)-1] <= spikes[len(spikes)-i-1]:
            repeat_regions.append((spikes[len(spikes)-i-1],len(ypoints)))


    # Remove enclosures
    for reg in repeat_regions:
        for reg2 in repeat_regions:
            if (reg[0] <= reg2[0] and reg[1] >reg2[1]) or (reg[0] < reg2[0] and reg[1] >= reg2[1]):
                repeat_regions.remove(reg2)
    
    # Calculate space between repeats and remove if too small
    for reg in repeat_regions:
        min_dist = len(ypoints)
        min_reg = (0,len(ypoints))
        for reg2 in repeat_regions:
            if reg2[0] > reg[1] and reg2[0] - reg[1] <= min_dist:
                min_dist = reg2[0] - reg[1]
                min_reg = reg2
        
        if min_dist < 20:
            repeat_regions.remove(reg2)
            repeat_regions.remove(reg)
            repeat_regions.insert(0, (reg[0], reg2[1]))

    return repeat_regions

def detect_repeats_2(ypoints):
    repeat_region = (0,0)
    spikes = [pos for pos in range(len(ypoints)) if ypoints[pos] > 0]
    dips = [pos for pos in range(len(ypoints)) if ypoints[pos] < 0]
    if len(spikes) > len(dips):
        repeat_region = (spikes[0], len(ypoints))
    elif len(dips) > len(spikes):
        repeat_region = (0, dips[len(dips)-1])
    
    return repeat_region
    


def remove_repeats(filename, str_profile, data=[], txt=True):
    seqs = []
    if data:
        for read in data:
            sequence_chars = [val for val in read]
            sequence = utils.parse_nucleotides(sequence_chars)
            seqs.append(np.array(sequence))
    else:
        if txt:
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
    print("building hashmap")
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

    errors = 0
    corrections = 0
    no_irr = 0
    new_seqs = []
    for target_index in range(1, len(seqs)):
        # Get maxcounts from counts
        target = seqs[target_index]
        print(str(target_index) + "/" + str(len(seqs)))
        xpoints = np.array([i for i in range(len(target) - f)])
        max_counts = []
        max_count_indices = []
        for i in range(len(xpoints)):
            maxp, argmaxp = get_maxcount(i, seqs_kmers, profile, target)
            max_counts.append(maxp)
            max_count_indices.append(argmaxp)

        # Get correction profile:
        start_end_posis = get_correction_profile(target, seqs, f, starts, ends)
        
        # Get diff profile:
        pre_corr_diff_profile = [max_counts[j] - max_counts[j-1] for j in range(1,len(max_counts))]
        
        # Apply correction strategy
        ypoints, correction_artifact = correct_counts(max_counts, max_count_indices, target, start_end_posis, seqs, pre_corr_diff_profile.copy(), profile)

        # detect repeats 
        #irregularities = [pos for pos in range(len(ypoints)) if np.abs(ypoints[pos]) > 1]
        repeat_region = detect_repeats_2(ypoints)
        if repeat_region != (0,0):
            start = 0
            end = 0
            if repeat_region[0] > 0:
                start = 0
                end = repeat_region[0]
            else:
                start = repeat_region[1]
                end = len(ypoints)

            read_wout_rep = target[start:end]
            new_seqs.append(read_wout_rep)
        
    print("corrections: " + str(corrections))
    print("errors: " + str(errors))
    print("wout irregularities: " + str(no_irr))
    return seqs
