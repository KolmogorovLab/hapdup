#!/usr/bin/env python3

import pysam
import sys


def project(bam_path, ref_seq, ref_pos):
    FLANK = 50
    name, pos, sign = project_flank(bam_path, ref_seq, ref_pos, 1)
    if name:
        return name, pos, sign
    else:
        return project_flank(bam_path, ref_seq, ref_pos, 50)


def project_flank(bam_path, ref_seq, ref_pos, flank):
    samfile = pysam.AlignmentFile(bam_path, "rb")

    #min_ref_diff = None
    last_aln = None
    last_qry_pos = None

    for pileup_col in samfile.pileup(ref_seq, max(0, ref_pos - flank), ref_pos + flank, truncate=True,
                                     max_depth=5, stepper="samtools"):
        #getting the longest alignment over this coordinate
        selected_pileup_aln = None
        largest_len = 0
        for pileup_read in pileup_col.pileups:
            if not pileup_read.is_del and not pileup_read.is_refskip:
                if pileup_read.alignment.query_alignment_length > largest_len:
                    selected_pileup_aln = pileup_read
                    largest_len = pileup_read.alignment.query_alignment_length

        if not selected_pileup_aln:
            continue

        #computing read position (not that simple heh)
        selected_aln = selected_pileup_aln.alignment
        left_hard = 0
        if selected_aln.cigartuples[0][0] == 5:
            left_hard = selected_aln.cigartuples[0][1]
        right_hard = 0
        if selected_aln.cigartuples[-1][0] == 5:
            right_hard = selected_aln.cigartuples[-1][1]
        query_length = selected_aln.infer_query_length() + left_hard + right_hard

        read_pos = selected_pileup_aln.query_position_or_next + left_hard
        if selected_aln.is_reverse:
            read_pos = query_length - read_pos

        last_aln = selected_pileup_aln.alignment
        last_qry_pos = read_pos
        ###

        #print(last_qry_pos)

        if ref_pos <= pileup_col.pos and last_aln is not None:
            break

    if not last_aln:
        return None, None, None

    return last_aln.query_name, last_qry_pos, 1 if not last_aln.is_reverse else 0


def bed_liftover(bed_file, bam_file, out_stream, output_failed=False):
    for line in open(bed_file, "r"):
        line = line.strip()
        if line.startswith("#"):
            out_stream.write(line + "\n")
            continue

        fields = line.split("\t")
        chr_id, chr_start, chr_end = fields[0], int(fields[1]), int(fields[2])

        proj_start_chr, proj_start_pos, proj_start_sign = project(bam_file, chr_id, chr_start)
        proj_end_chr, proj_end_pos, proj_end_sign = project(bam_file, chr_id, chr_end)

        #proj_start_chr, proj_start_pos, proj_start_sign = project(bam_file, chr_id, chr_start)
        #proj_end_chr, proj_end_pos, proj_end_sign = project(bam_file, chr_id, chr_end)

        if (not proj_start_chr or not proj_end_chr or
                proj_end_chr != proj_start_chr or
                proj_start_sign != proj_end_sign):
            if output_failed:
                print("#Failed:", line)
            continue

        if proj_start_sign < 0:
            proj_start_pos, proj_end_pos = proj_end_pos, proj_start_pos
        #if proj_end_pos < proj_start_pos:
        #    print(proj_start_pos, proj_start_sign, proj_end_pos, proj_end_sign)
        #    raise Exception("Negative length interval")

        fields[0], fields[1], fields[2] = proj_start_chr, str(proj_start_pos), str(proj_end_pos)
        out_stream.write("\t".join(fields) + "\n")


def main():
    if len(sys.argv) != 3:
        print("usage bed_liftover.py bed_file indexed_bam")
        return 1

    bed_file = sys.argv[1]
    bam_file = sys.argv[2]
    bed_liftover(bed_file, bam_file, sys.stdout, output_failed=True)


if __name__ == "__main__":
    main()
