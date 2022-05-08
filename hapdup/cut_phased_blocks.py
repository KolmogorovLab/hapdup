from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cut_phased_blocks(blocks_bed, fasta_in, fasta_out):
    block_borders = defaultdict(list)
    out_stream = open(fasta_out, "w")

    for line in open(blocks_bed, "r"):
        fields = line.strip().split()
        #if fields[3] != "ContigEnd":
        block_borders[fields[0]].append((int(fields[1]), int(fields[2])))

    for ctg in block_borders:
        block_borders[ctg].sort(key=lambda p: p[0])

    for seq in SeqIO.parse(fasta_in, "fasta"):
        if seq.id not in block_borders:
            SeqIO.write(seq, out_stream, "fasta")
        else:
            partition = []
            prev_split = 0
            for block_1, block_2 in zip(block_borders[seq.id][:-1], block_borders[seq.id][1:]):
                if block_2[0] > block_1[1]:
                    new_split = (block_2[0] + block_1[1]) // 2
                    if new_split > prev_split:
                        partition.append((prev_split, new_split))
                        prev_split = new_split

            if len(seq.seq) > prev_split:
                partition.append((prev_split, len(seq.seq)))

            for i, (begin, end) in enumerate(partition):
                seq_part = SeqRecord(seq.seq[begin:end], id=seq.id + "_phaseblock_" + str(i), description="")
                SeqIO.write(seq_part, out_stream, "fasta")

