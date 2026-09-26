"""Write a small synthetic RNA-seq study for exercising the pipeline end to end.

Six paired-end libraries (3 control, 3 heat) are simulated from random transcripts, with the
first genes up- and the next genes down-regulated under heat, plus the tx2gene table, a
RefSeq-style GFF and the two sample sheets the pipeline expects. Output is deterministic.
"""

import argparse
import csv
import gzip
import random
from pathlib import Path

GENE_ID_BASE = 105500000


def revcomp(sequence):
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genes", type=int, default=300)
    parser.add_argument("--mean-pairs", type=int, default=120, help="mean read pairs per gene per library")
    parser.add_argument("--read-length", type=int, default=75)
    args = parser.parse_args()
    rng = random.Random(20260925)
    out = args.output
    (out / "fastq").mkdir(parents=True, exist_ok=True)

    transcripts = {}
    with (out / "transcripts.fa").open("w") as fasta, (out / "tx2gene.tsv").open("w") as tx2gene:
        tx2gene.write("transcript_id\tgene_id\n")
        for index in range(args.genes):
            name = "XM_{:06d}.1".format(index)
            sequence = "".join(rng.choice("ACGT") for _ in range(rng.randint(600, 1400)))
            transcripts[name] = sequence
            fasta.write(">{}\n{}\n".format(name, sequence))
            tx2gene.write("{}\tNCBI:GeneID:{}\n".format(name, GENE_ID_BASE + index))
    with gzip.open(out / "genomic.gff.gz", "wt") as gff:
        gff.write("##gff-version 3\n")
        for index in range(args.genes):
            gene_id = GENE_ID_BASE + index
            gff.write(
                "NC_TEST\tRefSeq\tgene\t{}\t{}\t.\t+\t.\tID=gene-LOC{id};Dbxref=GeneID:{id};Name=LOC{id};"
                "gene_biotype=protein_coding\n".format(index * 2000 + 1, index * 2000 + 1500, id=gene_id)
            )

    up, down = range(0, args.genes // 10), range(args.genes // 10, args.genes // 5)
    baseline = {name: rng.lognormvariate(0, 0.8) for name in transcripts}
    samplesheet, design = [], []
    for condition in ("control", "heat"):
        for replicate in "ABC":
            sample = "{}_{}".format(condition, replicate)
            paths = [out / "fastq" / "{}_{}.fastq.gz".format(sample, mate) for mate in (1, 2)]
            with gzip.open(paths[0], "wt") as mate1, gzip.open(paths[1], "wt") as mate2:
                read_number = 0
                for index, (name, sequence) in enumerate(transcripts.items()):
                    fold = 4.0 if (condition == "heat" and index in up) else 0.25 if (
                        condition == "heat" and index in down) else 1.0
                    expected = args.mean_pairs * baseline[name] * fold * rng.uniform(0.85, 1.15)
                    for _ in range(int(rng.gauss(expected, expected ** 0.5)) if expected > 1 else 0):
                        fragment = rng.randint(200, 350)
                        start = rng.randint(0, len(sequence) - fragment)
                        piece = sequence[start:start + fragment]
                        read_number += 1
                        quality = "I" * args.read_length
                        mate1.write("@r{0}/1\n{1}\n+\n{2}\n".format(read_number, piece[:args.read_length], quality))
                        mate2.write("@r{0}/2\n{1}\n+\n{2}\n".format(
                            read_number, revcomp(piece[-args.read_length:]), quality))
            samplesheet.append({"sample": sample, "fastq_1": str(paths[0].resolve()),
                                "fastq_2": str(paths[1].resolve()), "strandedness": "auto"})
            design.append({"sample": sample, "condition": condition, "replicate": replicate,
                           "run_accession": "SIM{}".format(len(design) + 1)})
    for name, rows in (("samplesheet.csv", samplesheet), ("design.csv", design)):
        with (out / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print("wrote synthetic study to {}".format(out))


if __name__ == "__main__":
    main()
