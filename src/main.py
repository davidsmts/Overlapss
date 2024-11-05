import argparse
import overlapss
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Overlapss: Overlap based metagenome assembler. For now just difference profile correction")

parser.add_argument("file", type=str, help="Program requires data.")
parser.add_argument("target", type=int, help="Give index of read for which to display the count difference profile.")
parser.add_argument("--opt_kmer_profile", type=str, help="You can also specify a custom spaced k-mer profile to run this on.")

args = parser.parse_args()

profile = "1111110110110101110101011101011111101111"
if args.opt_kmer_profile:
    profile = args.opt_kmer_profile

print("calling overlapss correction profile")
xpoints, ypoints, _, _, _, _, _ = overlapss.correct_diff_profile(args.file, profile, args.target)
plt.plot(xpoints, ypoints)
plt.show()
