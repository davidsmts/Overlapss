import argparse
import overlapss
import matplotlib.pyplot as plt
import utils

parser = argparse.ArgumentParser(description="Overlapss: Overlap based metagenome assembler. For now just difference profile correction")

parser.add_argument("file", type=str, nargs="?", help="Program requires data.")
parser.add_argument("target", type=int, nargs="?", help="Give index of read for which to display the count difference profile.")
parser.add_argument("--remove_repeats", action="store_true", help="Remove repeats from all reads in file.")
parser.add_argument("--opt_kmer_profile", type=str, help="You can also specify a custom spaced k-mer profile to run this on.")
parser.add_argument("--opt_get_reads", type=str, help="If active, produces reads from the genome in the given file.")
parser.add_argument("--opt_reads_into", type=str, help="If active, writes read of genome into this file")
parser.add_argument("--opt_read_len", type=int, help="Length of reads we want to sample")
parser.add_argument("--opt_read_cov", type=int, help="Cov of reads we want to sample")

args = parser.parse_args()

profile = "1111110110110101110101011101011111101111"
if args.opt_kmer_profile:
    profile = args.opt_kmer_profile

if args.opt_get_reads and args.opt_reads_into:
    print("creating reads from " + args.opt_get_reads)
    reads = utils.get_even_reads_from_fasta(args.opt_get_reads, read_len=args.opt_read_len, amt=args.opt_read_cov)
    utils.write_reads(reads, args.opt_reads_into)

elif args.remove_repeats:
    seqs = overlapss.remove_repeats(args.file, profile)
    #plt.plot(xpoints, ypoints)
    #plt.show()

else:
    print("calling overlapss correction profile")
    xpoints, ypoints, _, _, _, _, _ = overlapss.correct_diff_profile(args.file, profile, args.target)
    plt.plot(xpoints, ypoints)
    plt.show()
