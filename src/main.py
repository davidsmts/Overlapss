import argparse
import overlapss
import matplotlib.pyplot as plt
import utils

parser = argparse.ArgumentParser(description="Overlapss: Overlap based metagenome assembler. For now just difference profile correction")

parser.add_argument("file", type=str, nargs="?", help="All commands require data.")
#parser.add_argument("targetfile", type=str, nargs="?", help="Most commands require a target file")
parser.add_argument("targetseq", type=int, nargs="?", help="Give index of read for which to display the count difference profile.")
parser.add_argument("--remove_repeats", action="store_true", help="Remove repeats from all reads in file.")
parser.add_argument("--opt_kmer_profile", type=str, help="You can also specify a custom spaced k-mer profile to run this on.")
parser.add_argument("--opt_get_reads", type=str, help="If active, produces reads from the genome in the given file.")
parser.add_argument("--opt_reads_into", type=str, help="If active, writes read of genome into this file")
parser.add_argument("--opt_read_len", type=int, help="Length of reads we want to sample")
parser.add_argument("--opt_read_cov", type=int, help="Cov of reads we want to sample")

args = parser.parse_args()

profile = "1111110110110101110101011101011111101111"
#111111111111111111111111111111
if args.opt_kmer_profile:
    profile = args.opt_kmer_profile

if args.opt_get_reads and args.opt_reads_into:
    print("creating reads from " + args.opt_get_reads)
    reads = utils.get_even_reads_from_fasta(args.opt_get_reads, read_len=args.opt_read_len, amt=args.opt_read_cov)
    utils.write_fasta(reads[0], args.opt_reads_into, nums=False)

elif args.remove_repeats:
    data = utils.load_fasta_seqs(args.file)
    correct_reads, chopped_reads = overlapss.remove_repeats(args.file, profile, data=data)
    utils.write_fasta(correct_reads, "../data/wout_reps/long_correct.fasta", nums=True)
    utils.write_fasta(chopped_reads, "../data/wout_reps/long_chopped.fasta", nums=True)
    #plt.plot(xpoints, ypoints)
    #plt.show()

else:
    print("calling overlapss correction profile")
    xpoints, ypoints, mc, _, _, _, _ = overlapss.correct_diff_profile(args.file, profile, args.targetseq)
    plt.plot([i for i in range(len(mc))], mc)
    plt.show()
