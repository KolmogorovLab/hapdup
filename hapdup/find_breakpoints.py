#!/usr/bin/env python3

"""
This script takes (phased) bam file as input, and outputs coordinates
of tentative 2-breaks, in particular inversion coordinates
"""

import sys
import re
import numpy as np
import math
from copy import copy
from collections import namedtuple, defaultdict
import pysam
from multiprocessing import Pool
import random
import argparse
import os
import logging

logger = logging.getLogger()

ReadSegment = namedtuple("ReadSegment", ["read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                                         "strand", "read_length", "haplotype"])
#ReadConnection = namedtuple("ReadConnection", ["ref_id", "position", "seg_1", "seg_2", "haplotype"])

class ReadConnection(object):
    __slots__ = "ref_id_1", "pos_1", "sign_1", "ref_id_2", "pos_2", "sign_2", "haplotype"

    def __init__(self, ref_id_1, pos_1, sign_1, ref_id_2, pos_2, sign_2, haplotype):
        self.ref_id_1 = ref_id_1
        self.ref_id_2 = ref_id_2
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.sign_1 = sign_1
        self.sign_2 = sign_2
        self.haplotype = haplotype

    def signed_coord_1(self):
        return self.sign_1 * self.pos_1

    def signed_coord_2(self):
        return self.sign_2 * self.pos_2


class Breakpoint(object):
    __slots__ = "ref_id", "position", "spanning_reads", "connections"

    def __init__(self, ref_id, ref_position):
        self.ref_id = ref_id
        self.position = ref_position
        self.spanning_reads = []
        self.connections =[]


class DoubleBreak(object):
    __slots__ = "bp_1", "direction_1", "bp_2", "direction_2", "connections"

    def __init__(self, bp_1, direction_1, bp_2, direction_2):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
        self.direction_2 = direction_2
        self.connections = []

    def directional_coord_1(self):
        return self.direction_1 * self.bp_1.position

    def directional_coord_2(self):
        return self.direction_2 * self.bp_2.position
        #self.spanning_reads = []

cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")
def get_segment(read_id, ref_id, ref_start, strand, cigar, haplotype):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """
    first_clip = False
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    #read_end = 0

    for token in cigar_parser.findall(cigar):
        op = token[-1]
        op_len = int(token[:-1])

        if op == "H" or op == "S":
            if not first_clip:
                first_clip = True
                read_start = op_len
            read_length += op_len

        if op == "M" or op == "=" or op == "X":
            read_aligned += op_len
            ref_aligned += op_len
            read_length += op_len
        if op == "D":
            ref_aligned += op_len
        if op == "I":
            read_aligned += op_len
            read_length += op_len

    ref_end = ref_start + ref_aligned
    read_end = read_start + read_aligned

    if strand == "-":
        read_start, read_end = read_length - read_end, read_length - read_start

    return ReadSegment(read_start, read_end, ref_start, ref_end, read_id,
                       ref_id, strand, read_length, haplotype)


def get_split_reads(bam_file, ref_id, inter_contig):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """
    split_reads = []

    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(ref_id, multiple_iterators=True):
        fields = aln.to_string().split()
        #fields = line.split()
        read_id, flags, chr_id, position = fields[0:4]
        cigar = fields[5]
        ref_id, ref_start = fields[2], int(fields[3])

        is_supplementary = int(flags) & 0x800
        is_secondary = int(flags) & 0x100
        is_unmapped = int(flags) & 0x4
        strand = "-" if int(flags) & 0x10 else "+"

        if is_supplementary or is_secondary or is_unmapped:
            continue

        sa_tag = None
        hp_tag = None
        for tag in fields[11:]:
            if tag.startswith("SA"):
                sa_tag = tag[5:]
            if tag.startswith("HP"):
                hp_tag = tag[5:]
        if not hp_tag:
            haplotype = 0
        else:
            haplotype = int(hp_tag)

        segments = []
        segments.append(get_segment(read_id, ref_id, ref_start, strand, cigar, haplotype))

        if sa_tag:
            for sa_aln in sa_tag.split(";"):
                if sa_aln:
                    sa_fields = sa_aln.split(",")
                    sa_ref, sa_ref_pos, sa_strand, sa_cigar = sa_fields[0], int(sa_fields[1]), sa_fields[2], sa_fields[3]

                    if sa_ref == ref_id or inter_contig:
                        segments.append(get_segment(read_id, sa_ref, sa_ref_pos, sa_strand, sa_cigar, haplotype))

        segments.sort(key=lambda s: s.read_start)
        split_reads.append(segments)

    return split_reads

def _unpacker(args):
    return get_split_reads(*args)

def get_all_reads_parallel(bam_file, num_threads, inter_contig):
    all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
    random.shuffle(all_reference_ids)
    tasks = [(bam_file, r, inter_contig) for r in all_reference_ids]

    fetched_reads = None
    with Pool(num_threads) as p:
        fetched_reads = p.map(_unpacker, tasks)

    all_reads = sum(fetched_reads, [])
    return all_reads


def resolve_overlaps(split_reads, min_ovlp_len):
    """
    Some supplementary alignments may be overlapping (e.g. in case of inversions with flanking repeat).
    This function checks if the overlap has ok structe, trims and outputs non-overlapping alignments
    """
    def _get_ovlp(seg_1, seg_2):
        max_ovlp_len = min(seg_1.read_end - seg_1.read_start, seg_2.read_end - seg_2.read_start) // 2
        if (seg_1.read_end - seg_2.read_start > min_ovlp_len and
            seg_1.read_end - seg_2.read_start < max_ovlp_len and
            seg_2.read_end > seg_1.read_end and
            seg_2.read_start > seg_1.read_start):

            return seg_1.read_end - seg_2.read_start
        else:
            return 0

    new_reads = []
    for read_segments in split_reads:
        upd_segments = []
        for i in range(len(read_segments)):
            left_ovlp = 0
            right_ovlp = 0
            if i > 0 and read_segments[i - 1].ref_id == read_segments[i].ref_id:
                left_ovlp = _get_ovlp(read_segments[i - 1], read_segments[i])
            if i < len(read_segments) - 1 and read_segments[i].ref_id == read_segments[i + 1].ref_id:
                right_ovlp = _get_ovlp(read_segments[i], read_segments[i + 1])

            left_ovlp = left_ovlp // 2
            right_ovlp = right_ovlp // 2
            seg = read_segments[i]

            if left_ovlp > 0:
                if seg.strand == "+":
                    seg = seg._replace(read_start = seg.read_start + left_ovlp,
                                       ref_start = seg.ref_start + left_ovlp)
                else:
                    seg = seg._replace(read_start = seg.read_start + left_ovlp,
                                       ref_end = seg.ref_end - left_ovlp)
            if right_ovlp > 0:
                if seg.strand == "+":
                    seg = seg._replace(read_end = seg.read_end - right_ovlp,
                                       ref_end = seg.ref_end - right_ovlp)
                else:
                    seg = seg._replace(read_end = seg.read_end - right_ovlp,
                                       ref_start = seg.ref_start + right_ovlp)

            upd_segments.append(seg)

        new_reads.append(upd_segments)

    return new_reads


def get_breakpoints(all_reads, split_reads, clust_len, min_reads, min_ref_flank, ref_lengths):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """
    seq_breakpoints = defaultdict(list)
    for read_segments in split_reads:
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            if abs(s1.read_end - s2.read_start) < clust_len:
                #if s1.ref_id != s2.ref_id:
                #    raise Exception("Inter-contig connection!")

                ref_bp_1 = s1.ref_end if s1.strand == "+" else s1.ref_start
                ref_bp_2 = s2.ref_start if s2.strand == "+" else s2.ref_end

                sign_1 = 1 if s1.strand == "+" else -1
                sign_2 = -1 if s2.strand == "+" else 1
                #conn_1 = ref_bp_1 if s1.strand == "+" else -ref_bp_1
                #conn_2 = -ref_bp_2 if s2.strand == "+" else ref_bp_2

                #seq_breakpoints[s1.ref_id].append(ReadConnection(s1.ref_id, ref_bp_1, conn_1, conn_2, s1.haplotype))
                #seq_breakpoints[s2.ref_id].append(ReadConnection(s2.ref_id, ref_bp_2, conn_2, conn_1, s1.haplotype))

                seq_breakpoints[s1.ref_id].append(ReadConnection(s1.ref_id, ref_bp_1, sign_1, s2.ref_id, ref_bp_2, sign_2, s1.haplotype))
                seq_breakpoints[s2.ref_id].append(ReadConnection(s2.ref_id, ref_bp_2, sign_2, s1.ref_id, ref_bp_1, sign_1, s2.haplotype))

    bp_clusters = defaultdict(list)
    for seq, bp_pos in seq_breakpoints.items():
        clusters = []
        cur_cluster = []
        bp_pos.sort(key=lambda bp: bp.pos_1)

        for rc in bp_pos:
            if cur_cluster and rc.pos_1 - cur_cluster[-1].pos_1 > clust_len:
                clusters.append(cur_cluster)
                cur_cluster = [rc]
            else:
                cur_cluster.append(rc)
        if cur_cluster:
            clusters.append(cur_cluster)

        for cl in clusters:
            if len(cl) >= min_reads:
                position = int(np.median([x.pos_1 for x in cl]))
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank:
                    bp_cluster = Breakpoint(seq, position)
                    bp_cluster.connections = cl
                    bp_clusters[seq].append(bp_cluster)

    #find reads that span putative breakpoints
    for read_segments in all_reads:
        for seg in read_segments:
            for bp in bp_clusters[seg.ref_id]:
                if bp.position - seg.ref_start > clust_len and seg.ref_end - bp.position > clust_len:
                    bp.spanning_reads.append(seg)

    return bp_clusters


def get_2_breaks(bp_clusters, clust_len, min_connections):
    """
    Matches one-sided breakpoints into 2-breaks, complementery pairs of one-sided breakpoints
    """
    def _normalize_coord(coord, clusters):
        for cl in clusters:
            if abs(abs(coord) - cl.position) < clust_len:
                return int(math.copysign(1, coord)), cl
        return None, None

    #Breakpoints are clustered by a single reference coordinate.
    #Our goal here is to translate it into connections (defined for two directional reference coordinates)
    #and then find balanced breakes (2-breaks)

    double_breaks = []
    double_connections = defaultdict(list)
    for seq in bp_clusters:
        for cl in bp_clusters[seq]:
            for conn in cl.connections:
                #normalizing the coordinates of the original read, wrt to clustered breakpoint coordinates
                dir_1, bp_1 = _normalize_coord(conn.signed_coord_1(), bp_clusters[conn.ref_id_1])
                dir_2, bp_2 = _normalize_coord(conn.signed_coord_2(), bp_clusters[conn.ref_id_2])
                if None in [bp_1, bp_2]:
                    continue

                #same breakpoint, just connected from two different sides
                if bp_1 == bp_2:
                    continue

                if bp_2.position < bp_1.position:
                    bp_1, bp_2 = bp_2, bp_1
                    dir_1, dir_2 = dir_2, dir_1
                double_connections[(bp_1, dir_1, bp_2, dir_2)].append(conn)

    for (bp_1, dir_1, bp_2, dir_2), conn_list in double_connections.items():
        if len(conn_list) // 2 >= min_connections:
            double_breaks.append(DoubleBreak(bp_1, dir_1, bp_2, dir_2))
            double_breaks[-1].connections = conn_list

    double_breaks.sort(key=lambda b: (b.bp_1.ref_id, b.bp_1.position, b.direction_1))

    #find breakend pairs
    all_breakends = set()
    balanced_breaks = []
    for db in double_breaks:
        all_breakends.add(db.directional_coord_1())
        all_breakends.add(db.directional_coord_2())

    for db in double_breaks:
        if -db.directional_coord_1() in all_breakends and -db.directional_coord_2() in all_breakends:
            balanced_breaks.append(db)

    return double_breaks, balanced_breaks


def output_breaks(breaks, out_stream):
    header = "#seq_id_1\tpos_1\tdirection_1\tseq_id_2\tposition_2\tdirection_2\thaplotype\tsupporting_reads\topposing_reads"
    out_stream.write(header + "\n")
    for br in breaks:
        by_hp = defaultdict(int)
        for conn in br.connections:
            by_hp[conn.haplotype] += 1
        max_hp = max(by_hp, key=by_hp.get)

        num_opposing = 0
        for r in br.bp_1.spanning_reads + br.bp_2.spanning_reads:
            if r.haplotype == max_hp:
                num_opposing += 1

        out_stream.write("\t".join([br.bp_1.ref_id, str(br.bp_1.position), "+" if br.direction_1 > 0 else "-",
                                   br.bp_2.ref_id, str(br.bp_2.position), "+" if br.direction_2 > 0 else "-",
                                   str(max_hp), str(by_hp[max_hp] // 2), str(num_opposing)]) + "\n")


def output_inversions(breaks, out_stream):
    out_stream.write("#chrom\tchrmStart\tchrmEnd\thaplotype\tsupportReads\topposingReads\n")
    for br in breaks:
        if br.bp_1.ref_id == br.bp_2.ref_id and br.direction_1 > 0 and br.direction_2 > 0:
            matching = False
            for other_br in breaks:
                if (other_br.bp_1.ref_id == br.bp_1.ref_id and other_br.bp_2.ref_id == br.bp_2.ref_id and
                        other_br.bp_1.position == br.bp_1.position and other_br.bp_2.position == br.bp_2.position and
                        other_br.direction_1 < 0 and other_br.direction_2 < 0):
                    matching = True
                    break

            if matching:
                by_hp = defaultdict(int)
                for conn in br.connections:
                    by_hp[conn.haplotype] += 1
                max_hp = max(by_hp, key=by_hp.get)

                num_opposing = 0
                for r in br.bp_1.spanning_reads + br.bp_2.spanning_reads:
                    if r.haplotype == max_hp:
                        num_opposing += 1
                out_stream.write("\t".join([br.bp_1.ref_id, str(br.bp_1.position), str(br.bp_2.position),
                                            str(max_hp), str(by_hp[max_hp] // 2), str(num_opposing)]) + "\n")


def _run_pipeline(arguments):
    BP_CLUSTER_SIZE = 100
    MIN_BREAKPOINT_READS = 10
    MIN_DOUBLE_BP_READS = 5
    MIN_REF_FLANK = 5000

    parser = argparse.ArgumentParser \
        (description="Find breakpoints given a bam file")

    parser.add_argument("-b", "--bam", dest="bam_path",
                        metavar="path", required=True,
                        help="path to assembly bam file (must be indexed)")
    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("-t", "--threads", dest="threads",
                        default=8, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument("--cluster-size", dest="cluster_size",
                        default=BP_CLUSTER_SIZE, metavar="int", type=int, help="size of breakpoint cluster in bp [100]")
    parser.add_argument("--breakpoint-min-reads", dest="bp_min_reads",
                        default=MIN_BREAKPOINT_READS, metavar="int", type=int, help="minimum reads in breakpoint [10]")
    parser.add_argument("--double-breakpoint-min-reads", dest="double_bp_min_reads",
                        default=MIN_DOUBLE_BP_READS, metavar="int", type=int, help="minimum reads in double breakpoint [5]")
    parser.add_argument("--min-reference-flank", dest="min_ref_flank",
                        default=MIN_REF_FLANK, metavar="int", type=int,
                        help="minimum distance between breakpoint and sequence ends [5000]")
    args = parser.parse_args(arguments)

    with pysam.AlignmentFile(args.bam_path, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    all_reads = get_all_reads_parallel(args.bam_path, args.threads, inter_contig=True)
    split_reads = []
    for r in all_reads:
        if len(r) > 1:
            split_reads.append(r)
    logger.info("Parsed %d reads %d split reads", len(all_reads), len(split_reads))

    split_reads = resolve_overlaps(split_reads, args.cluster_size)
    bp_clusters = get_breakpoints(all_reads, split_reads, args.cluster_size, args.bp_min_reads, args.min_ref_flank, ref_lengths)
    all_breaks, balanced_breaks = get_2_breaks(bp_clusters, args.cluster_size, args.double_bp_min_reads)

    out_breaks = os.path.join(args.out_dir, "breakpoints_all.csv")
    out_balanced_breaks = os.path.join(args.out_dir, "breakpoints_balanced.csv")
    out_inversions = os.path.join(args.out_dir, "inversions.bed")

    output_breaks(all_breaks, open(out_breaks, "w"))
    output_breaks(balanced_breaks, open(out_balanced_breaks, "w"))
    output_inversions(balanced_breaks, open(out_inversions, "w"))


def find_breakpoints(input_bam, output_dir, num_threads):
    _run_pipeline(["-b", input_bam, "-o", output_dir, "-t", str(num_threads)])


def main():
    _run_pipeline(sys.argv[1:])


if __name__ == "__main__":
    main()
