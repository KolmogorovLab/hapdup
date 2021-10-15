#!/usr/bin/env python3

import os
import sys
import argparse
from distutils import spawn
import subprocess
import threading

from hap_dup.find_breakpoints import find_breakpoints
from hap_dup.bed_liftover import bed_liftover
from hap_dup.apply_inversions import apply_inversions

pipeline_dir = os.path.dirname(os.path.realpath(__file__))
ALN_FILTER = sys.executable + " " + os.path.join(pipeline_dir, "filter_misplaced_alignments.py")
MARGIN = "/card/tools/margin/build/margin"
MARGIN_PARAMS = "/card/tools/margin/params/misc/allParams.ont_haplotag.sv.json"
FLYE = "/card/projects/Flye/bin/flye"
SAMTOOLS = "samtools"
MINIMAP = "minimap2"

FS_NAME = "/card"
#DOCKER_ID = "kishwars/pepper_deepvariant:test-r0.5-fix"
DOCKER_ID = "kishwars/pepper_deepvariant:test-v0.5"
#MODEL_PATH = "/opt/pepper_models/PEPPER_VARIANT_R941_ONT_V5.pkl"
MODEL_PATH = "/opt/pepper_models/PEPPER_SNP_R941_ONT_V4.pkl"


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

    parser.add_argument("-t", "--threads", dest="threads",
                        default=10, metavar="int", help="number of parallel threads [10]")
    args = parser.parse_args()

    for e in [SAMTOOLS, "docker", FLYE, MARGIN, MINIMAP]:
        if not spawn.find_executable(e):
            print("Not installed: " + e)
            return 1

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    def file_check(path):
        if not os.path.isfile(path):
            raise Exception("Missing output:", path)

    filtered_bam = os.path.join(args.out_dir, "filtered.bam")
    pepper_dir = os.path.join(args.out_dir, "pepper")
    pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_SNP_OUTPUT.vcf"))
    #pepper_vcf = os.path.abspath(os.path.join(pepper_dir, "PEPPER_VARIANT_OUTPUT_PHASING.vcf"))

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
        print("Skipped filtering phase")
    else:
        filter_cmd = [SAMTOOLS, "view", "-h", "-@4", args.bam, "|", ALN_FILTER, "|", SAMTOOLS, "view", "-", "-b", "-1", "-@4",">", filtered_bam]
        print("Running:", " ".join(filter_cmd))
        subprocess.check_call(" ".join(filter_cmd), shell=True)
        file_check(filtered_bam)

        index_cmd = [SAMTOOLS, "index", "-@4", filtered_bam]
        print("Running:", " ".join(index_cmd))
        subprocess.check_call(" ".join(index_cmd), shell=True)
        overwrite = True

    #STAGE 2: Run PEPPER
    if os.path.isfile(pepper_vcf) and not overwrite:
        print("Skipped pepper phase")
    else:
        pepper_log = os.path.join(pepper_dir, "pepper.log")
        if not os.path.isdir(pepper_dir):
            os.mkdir(pepper_dir)
        pepper_cmd = ["docker", "run", "-it", "-v", "{0}:{0}".format(FS_NAME), "-u", "`id -u`:`id -g`", "--ipc", "host",
                      DOCKER_ID, "pepper_variant", "call_variant", "-b", os.path.abspath(filtered_bam), "-f", os.path.abspath(args.assembly),
                      "-o", os.path.abspath(pepper_dir), "-m", MODEL_PATH, "-t", str(args.threads), "-s", "Sample", "-w", "4", "-bs", "64",
                      "2>&1", "|tee", pepper_log]
                      #"-o", os.path.abspath(pepper_dir), "-m", MODEL_PATH, "-t", str(args.threads), "-s", "Sample", "--allow_supplementary"]
        print("Running:", " ".join(pepper_cmd))
        subprocess.check_call(" ".join(pepper_cmd), shell=True)
        subprocess.call("rm -r " + os.path.join(pepper_dir, "images*"), shell=True)
        subprocess.call("rm -r " + os.path.join(pepper_dir, "predictions*"), shell=True)
        file_check(pepper_vcf)
        overwrite = True

    #STAGE 3: Phase with Margin
    if os.path.isfile(haplotagged_bam) and not overwrite:
        print("Skipped margin phase")
    else:
        margin_log = os.path.join(margin_dir, "margin.log")
        if not os.path.isdir(margin_dir):
            os.mkdir(margin_dir)

        margin_cmd = [MARGIN, "phase", os.path.abspath(filtered_bam), os.path.abspath(args.assembly), pepper_vcf, MARGIN_PARAMS,
                      "-t", str(args.threads), "-o", os.path.abspath(os.path.join(margin_dir, "MARGIN_PHASED")),
                      "2>&1", "|tee", margin_log]
        print("Running:", " ".join(margin_cmd))
        subprocess.check_call(" ".join(margin_cmd), shell=True)
        file_check(haplotagged_bam)
        #subprocess.call("rm " + os.path.abspath(filtered_bam), shell=True)

        index_cmd = [SAMTOOLS, "index", "-@4", haplotagged_bam]
        print("Running:", " ".join(index_cmd))
        subprocess.check_call(" ".join(index_cmd), shell=True)
        overwrite = True

    #STAGE 4: polish haplotypes with Flye
    if all(map(os.path.isfile, polished_flye_hap.values())) and not overwrite:
        print("Skipped Flye phase")
    else:
        def run_flye_hp(hp):
            threads = max(1, int(args.threads) // 2)
            flye_out = os.path.join(args.out_dir, "flye_hap_{0}".format(hp))
            flye_cmd = [FLYE, "--polish-target", os.path.abspath(args.assembly), "--nano-raw",  haplotagged_bam, "-t", str(threads),
                        "-o", flye_out, "--polish-haplotypes", "0,{}".format(hp), "2>/dev/null"]
            print("Running:", " ".join(flye_cmd))
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
    print("Finding breakpoints")
    find_breakpoints(haplotagged_bam, structural_dir, args.threads)

    for hp in [1, 2]:
        minimap_out = os.path.join(structural_dir, "liftover_hp{0}.bam".format(hp))
        minimap_cmd = [MINIMAP, "-ax", "asm5", "-t", str(args.threads), "-K", "5G", args.assembly, polished_flye_hap[hp], "2>/dev/null", "|",
                       SAMTOOLS, "sort", "-m", "4G", "-@4", ">", minimap_out]
        print("Running:", " ".join(minimap_cmd))
        subprocess.check_call(" ".join(minimap_cmd), shell=True)
        subprocess.check_call("samtools index -@ 4 {0}".format(minimap_out), shell=True)

        inversions_hp = os.path.join(structural_dir, "inversions_hp{0}.bed".format(hp))
        bed_liftover(inversions_bed, minimap_out, open(inversions_hp, "w"))

        phased_blocks_hp = os.path.join(args.out_dir, "phased_blocks_hp{0}.bed".format(hp))
        bed_liftover(phased_blocks_bed, minimap_out, open(phased_blocks_hp, "w"))

        apply_inversions(inversions_hp, polished_flye_hap[hp], final_haplotype[hp], hp)


if __name__ == "__main__":
    main()
