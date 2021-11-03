#!/usr/bin/env python3

import os
import sys
import argparse
from distutils import spawn
import subprocess
import threading
import shutil

from hapdup.find_breakpoints import find_breakpoints
from hapdup.bed_liftover import bed_liftover
from hapdup.apply_inversions import apply_inversions
from hapdup.filter_misplaced_alignments import filter_alignments_parallel


MARGIN = "margin"
FLYE = "flye"
SAMTOOLS = "flye-samtools"
MINIMAP = "flye-minimap2"
PEPPER_VARIANT = "pepper_variant"


PEPPER_MODEL = os.environ["PEPPER_MODEL"]
MARGIN_PARAMS = os.environ["MARGIN_MODEL"]

def main():
    parser = argparse.ArgumentParser \
        (description="Reassemble haplotypes from collapsed haploid assmebly")

    parser.add_argument("--assembly", dest="assembly",
                        metavar="path", required=True,
                        help="path to haploid assembly (contigs in fasta format)")
    parser.add_argument("--bam", dest="bam", required=True, metavar="path",
                        default=None, help="path to the alignment of reads on the assembly in bam format")
    parser.add_argument("--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("--overwrite", dest="overwrite",
                        default=False, action="store_true",
                        help="Do not attempt to restart from complete phases, overwrites existing results")

    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=10, metavar="int", help="number of parallel threads [10]")
    args = parser.parse_args()

    for e in [SAMTOOLS, FLYE, MARGIN, MINIMAP, PEPPER_VARIANT]:
        if not spawn.find_executable(e):
            print("Not installed: " + e, file=sys.stderr)
            return 1

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    def file_check(path):
        if not os.path.isfile(path):
            raise Exception("Missing output:", path)

    file_check(args.bam)
    file_check(args.assembly)

    filtered_bam = os.path.join(args.out_dir, "filtered.bam")
    pepper_dir = os.path.join(args.out_dir, "pepper")
    #pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_SNP_OUTPUT.vcf"))
    #pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_OUTPUT_PHASING.vcf"))
    pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_FULL.vcf"))

    margin_dir = os.path.join(args.out_dir, "margin")
    haplotagged_bam = os.path.join(margin_dir, "MARGIN_PHASED") + ".haplotagged.bam"
    phased_blocks_bed = os.path.join(margin_dir, "MARGIN_PHASED") + ".phaseset.bed"

    polished_flye_hap = {1 : os.path.join(args.out_dir, "flye_hap_1", "polished_1.fasta"),
                         2 : os.path.join(args.out_dir, "flye_hap_2", "polished_1.fasta")}

    structural_dir = os.path.join(args.out_dir, "structural")
    inversions_bed = os.path.join(structural_dir, "inversions.bed")

    final_haplotype = {1 : os.path.join(args.out_dir, "haplotype_1.fasta"),
                       2 : os.path.join(args.out_dir, "haplotype_2.fasta")}

    #STAGE 1: filter suspicious alignments, index the resulting bam
    overwrite = args.overwrite
    if (os.path.isfile(filtered_bam) or os.path.isfile(haplotagged_bam)) and not overwrite:
        print("Skipped filtering phase", file=sys.stderr)
    else:
        print("Filtering alignments", file=sys.stderr)
        filter_alignments_parallel(args.bam, filtered_bam, min(args.threads, 30))
        file_check(filtered_bam)

        index_cmd = [SAMTOOLS, "index", "-@4", filtered_bam]
        print("Running:", " ".join(index_cmd), file=sys.stderr)
        subprocess.check_call(" ".join(index_cmd), shell=True)
        overwrite = True

    #STAGE 2: Run PEPPER
    if os.path.isfile(pepper_vcf) and not overwrite:
        print("Skipped pepper phase", file=sys.stderr)
    else:
        pepper_log = os.path.join(pepper_dir, "pepper.log")
        if not os.path.isdir(pepper_dir):
            os.mkdir(pepper_dir)

        model_copy = os.path.join(pepper_dir, "pepper_model.bin")
        shutil.copyfile(PEPPER_MODEL, model_copy)

        pepper_cmd = [PEPPER_VARIANT, "call_variant", "-b", os.path.abspath(filtered_bam), "-f", os.path.abspath(args.assembly),
                      "-o", os.path.abspath(pepper_dir), "-m", model_copy, "-t", str(args.threads), "-s", "Sample", "--ont_r9_guppy5_sup",
                      "--include-supplementary", "2>&1", "|tee", pepper_log]

        print("Running:", " ".join(pepper_cmd), file=sys.stderr)
        subprocess.check_call(" ".join(pepper_cmd), shell=True)
        subprocess.call("rm -r " + os.path.join(pepper_dir, "images*"), shell=True)
        subprocess.call("rm -r " + os.path.join(pepper_dir, "predictions*"), shell=True)
        file_check(pepper_vcf)
        overwrite = True

    #STAGE 3: Phase with Margin
    if os.path.isfile(haplotagged_bam) and not overwrite:
        print("Skipped margin phase", file=sys.stderr)
    else:
        margin_log = os.path.join(margin_dir, "margin.log")
        if not os.path.isdir(margin_dir):
            os.mkdir(margin_dir)

        margin_cmd = [MARGIN, "phase", os.path.abspath(filtered_bam), os.path.abspath(args.assembly), pepper_vcf, MARGIN_PARAMS,
                      "-t", str(args.threads), "-o", os.path.abspath(os.path.join(margin_dir, "MARGIN_PHASED")),
                      "2>&1", "|tee", margin_log]
        print("Running:", " ".join(margin_cmd), file=sys.stderr)
        subprocess.check_call(" ".join(margin_cmd), shell=True)
        file_check(haplotagged_bam)
        #subprocess.call("rm " + os.path.abspath(filtered_bam), shell=True)

        index_cmd = [SAMTOOLS, "index", "-@4", haplotagged_bam]
        print("Running:", " ".join(index_cmd), file=sys.stderr)
        subprocess.check_call(" ".join(index_cmd), shell=True)
        overwrite = True

    #STAGE 4: polish haplotypes with Flye
    if all(map(os.path.isfile, polished_flye_hap.values())) and not overwrite:
        print("Skipped Flye phase", file=sys.stderr)
    else:
        def run_flye_hp(hp):
            threads = max(1, int(args.threads) // 2)
            flye_out = os.path.join(args.out_dir, "flye_hap_{0}".format(hp))
            flye_cmd = [FLYE, "--polish-target", os.path.abspath(args.assembly), "--nano-raw",  haplotagged_bam, "-t", str(threads),
                        "-o", flye_out, "--polish-haplotypes", "0,{}".format(hp), "2>/dev/null"]
            print("Running:", " ".join(flye_cmd), file=sys.stderr)
            subprocess.check_call(" ".join(flye_cmd), shell=True)

        threads = []
        for hp in [1, 2]:
            threads.append(threading.Thread(target=run_flye_hp, args=(hp,)))
            threads[-1].start()
        for t in threads:
            t.join()
        overwrite = True

    #STAGE 5: structural polishing
    if not os.path.isdir(structural_dir):
        os.mkdir(structural_dir)
    print("Finding breakpoints", file=sys.stderr)
    find_breakpoints(haplotagged_bam, structural_dir, args.threads)

    for hp in [1, 2]:
        minimap_out = os.path.join(structural_dir, "liftover_hp{0}.bam".format(hp))
        minimap_cmd = [MINIMAP, "-ax", "asm5", "-t", str(args.threads), "-K", "5G", args.assembly, polished_flye_hap[hp], "2>/dev/null", "|",
                       SAMTOOLS, "sort", "-m", "4G", "-@4", ">", minimap_out]
        print("Running:", " ".join(minimap_cmd), file=sys.stderr)
        subprocess.check_call(" ".join(minimap_cmd), shell=True)
        subprocess.check_call(SAMTOOLS + " index -@ 4 {0}".format(minimap_out), shell=True)

        inversions_hp = os.path.join(structural_dir, "inversions_hp{0}.bed".format(hp))
        bed_liftover(inversions_bed, minimap_out, open(inversions_hp, "w"))

        phased_blocks_hp = os.path.join(args.out_dir, "phased_blocks_hp{0}.bed".format(hp))
        bed_liftover(phased_blocks_bed, minimap_out, open(phased_blocks_hp, "w"))

        apply_inversions(inversions_hp, polished_flye_hap[hp], final_haplotype[hp], hp)


if __name__ == "__main__":
    main()
