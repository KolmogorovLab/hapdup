
import sys
import hap_dup.fasta_parser as fp
from collections import defaultdict


def apply_inversions(bed_inversions, input_fasta, output_fasta, haplotype):
    inversions = defaultdict(list)
    for line in open(bed_inversions, "r"):
        line = line.strip()
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        chr_id, start, end, hap = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
        if hap == haplotype:
            inversions[chr_id].append((start, end))

    num_applied = 0
    fasta_dict = fp.read_sequence_dict(input_fasta)
    for seq, inv_list in inversions.items():
        for inv_start, inv_end in inv_list:
            fasta_dict[seq] = fasta_dict[seq][:inv_start] + fp.reverse_complement(fasta_dict[seq][inv_start:inv_end]) + fasta_dict[seq][inv_end:]
            num_applied += 1

    print("Applied", num_applied, "inversions for haplotype", haplotype, file=sys.stderr)
    fp.write_fasta_dict(fasta_dict, output_fasta)


if __name__ == "__main__":
    apply_inversions(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
