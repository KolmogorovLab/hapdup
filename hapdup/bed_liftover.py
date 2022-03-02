#!/usr/bin/env python3

import pysam
import sys
from multiprocessing import Pool
import random


def project(bam_path, ref_seq, ref_pos):
    FLANK = 50
    name, pos, sign = project_flank(bam_path, ref_seq, ref_pos, 1)
    if name:
        return name, pos, sign
    else:
        return project_flank(bam_path, ref_seq, ref_pos, 50)


def project_flank(bam_file, ref_seq, ref_pos, flank):
    bam_handle = pysam.AlignmentFile(bam_file, "rb")

    #min_ref_diff = None
    last_aln = None
    last_qry_pos = None

    #print("flank", flank)
    for pileup_col in bam_handle.pileup(ref_seq, max(0, ref_pos - flank), ref_pos + flank, truncate=True,
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


def bed_liftover(bed_file, bam_file, output_failed, ctg_id):
    out_strings = []
    #print("Processing", ctg_id)

    for line in open(bed_file, "r"):
        line = line.strip()
        if line.startswith("#"):
            out_stream.write(line + "\n")
            continue

        fields = line.split("\t")
        chr_id, chr_start, chr_end = fields[0], int(fields[1]), int(fields[2])

        if chr_id != ctg_id:
            continue

        MIN_BLOCK = 1000
        if chr_end - chr_start < MIN_BLOCK:
            continue

        proj_start_chr, proj_start_pos, proj_start_sign = project(bam_file, chr_id, chr_start)
        proj_end_chr, proj_end_pos, proj_end_sign = project(bam_file, chr_id, chr_end)

        if (not proj_start_chr or not proj_end_chr or
                proj_end_chr != proj_start_chr or
                proj_start_sign != proj_end_sign):
            if output_failed:
                out_strings.append("#Failed: " + line)
                #print("#Failed:", line)
            continue

        if proj_start_sign < 0:
            proj_start_pos, proj_end_pos = proj_end_pos, proj_start_pos
        #if proj_end_pos < proj_start_pos:
        #    print(proj_start_pos, proj_start_sign, proj_end_pos, proj_end_sign)
        #    raise Exception("Negative length interval")

        fields[0], fields[1], fields[2] = proj_start_chr, str(proj_start_pos), str(proj_end_pos)
        out_strings.append("\t".join(fields))

    return ctg_id, out_strings


def _unpacker(args):
    return bed_liftover(*args)


def liftover_parallel(bed_file, bam_file, out_stream, output_failed, num_threads):
    all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
    random.shuffle(all_reference_ids)
    tasks = [(bed_file, bam_file, output_failed, r)
                for r in all_reference_ids]

    thread_outputs = None
    with Pool(num_threads) as p:
        thread_outputs = p.map(_unpacker, tasks)

    thread_outputs.sort(key=lambda o: o[0])
    thread_outputs = [x[1] for x in thread_outputs if len(x[1]) > 0]
    bed_strings = sum(thread_outputs, [])
    for s in bed_strings:
        out_stream.write(s + "\n")


def main():
    if len(sys.argv) != 3:
        print("usage bed_liftover.py bed_file indexed_bam")
        return 1

    bed_file = sys.argv[1]
    bam_file = sys.argv[2]
    liftover_parallel(bed_file, bam_file, sys.stdout, False, 10)


if __name__ == "__main__":
    main()
