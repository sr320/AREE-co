"""Test whether RNA-seq "replicates" come from different animals, using expressed-SNP allele frequencies.

Oysters are highly heterozygous, so libraries from different animals disagree at many expressed
SNPs, while libraries made from one RNA sample agree to within sampling noise. The check needs no
aligner beyond Salmon: its selective-alignment SAM (``--writeMappings``) is piled up over a set of
well-expressed single-isoform transcripts.

    python scripts/check_replicate_genotypes.py targets --quant-dir ANALYSIS/salmon \\
        --tx2gene REF/GCF_963853765.1_tx2gene.tsv --output targets.txt
    # per library, on a read subsample (4M pairs is plenty):
    salmon quant -i REF/salmon_index_decoy -l A -p 8 -1 <(gzip -dc R1 | head -16000000) \\
        -2 <(gzip -dc R2 | head -16000000) --writeMappings=/dev/stdout -o /tmp/sq 2>/dev/null \\
      | python scripts/check_replicate_genotypes.py pileup --targets targets.txt \\
        --transcripts REF/GCF_963853765.1_xbMagGiga1.1_rna.fna.gz --output geno/SAMPLE.npz
    python scripts/check_replicate_genotypes.py compare geno/*.npz

``compare`` groups libraries by their name minus the trailing ``_<replicate>`` and reports the mean
allele-frequency difference within and between groups next to a binomial simulation of two
libraries drawn from the same RNA at the observed depths.
"""

import argparse
import csv
import gzip
import itertools
import sys
from collections import Counter
from pathlib import Path

import numpy as np


READ_LENGTH_CIGAR = "150M"
MAX_MISMATCHES = 3


def select_targets(quant_dir, tx2gene_path, output, count=1500, min_length=1000):
    """Pick well-covered transcripts that are the only isoform of their gene, so reads map uniquely."""
    tx2gene = {}
    with open(tx2gene_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tx2gene[row["transcript_id"]] = row["gene_id"]
    isoforms = Counter(tx2gene.values())
    reads, lengths = Counter(), {}
    for quant in sorted(Path(quant_dir).glob("*/quant.sf")):
        with quant.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                reads[row["Name"]] += float(row["NumReads"])
                lengths[row["Name"]] = int(row["Length"])
    candidates = [t for t in reads if t in tx2gene and isoforms[tx2gene[t]] == 1
                  and lengths[t] >= min_length and not tx2gene[t].startswith("MT")]
    candidates.sort(key=lambda t: -reads[t] / lengths[t])
    Path(output).write_text("\n".join(candidates[:count]) + "\n")


def _read_fasta(path, keep):
    sequences, current, chunks = {}, None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    sequences[current] = "".join(chunks)
                name = line[1:].split()[0]
                current, chunks = (name if name in keep else None), []
            elif current:
                chunks.append(line.strip())
    if current:
        sequences[current] = "".join(chunks)
    return sequences


def pileup(sam, targets_path, transcripts_path, output):
    targets = {line.strip() for line in open(targets_path) if line.strip()}
    sequences = _read_fasta(transcripts_path, targets)
    lookup = np.full(256, 4, np.uint8)
    for index, base in enumerate("ACGT"):
        lookup[ord(base)] = index
    reference = {t: lookup[np.frombuffer(s.upper().encode(), np.uint8)] for t, s in sequences.items()}
    counts = {t: np.zeros((len(s), 4), np.uint32) for t, s in sequences.items()}
    used = mismatches = bases = 0
    for line in sam:
        if line.startswith("@"):
            continue
        fields = line.split("\t", 11)
        transcript = fields[2]
        # Salmon reports every alignment as a full-length match, so indels and clipped ends only
        # show up as mismatches; reads with more than a few are dropped rather than trusted.
        if transcript not in counts or fields[5] != READ_LENGTH_CIGAR or int(fields[1]) & 0x904:
            continue
        if "NH:i:1\t" not in fields[11]:
            continue
        start, read = int(fields[3]) - 1, fields[9]
        if start < 0 or start + len(read) > counts[transcript].shape[0]:
            continue
        bases_read = lookup[np.frombuffer(read.encode(), np.uint8)]
        differences = int((bases_read != reference[transcript][start:start + len(read)]).sum())
        if differences > MAX_MISMATCHES:
            continue
        called = bases_read < 4
        counts[transcript][np.arange(start, start + len(read))[called], bases_read[called]] += 1
        used += 1
        mismatches += differences
        bases += len(read)
    np.savez_compressed(output, **counts)
    sys.stderr.write(f"reads used {used}; mismatch rate {mismatches / max(bases, 1):.4f}\n")


def compare(paths, min_depth=30, min_pooled_minor=0.15, seed=1):
    names = [Path(p).name[:-len(".npz")] for p in paths]
    data = [np.load(p) for p in paths]
    frequencies, depths = [], []
    for transcript in data[0].files:
        stacked = np.stack([d[transcript] for d in data]).astype(float)
        depth = stacked.sum(2)
        total = stacked.sum(0)
        minor = np.argsort(-total, 1)[:, 1]
        positions = np.arange(stacked.shape[1])
        pooled_minor = total[positions, minor] / np.maximum(total.sum(1), 1)
        keep = (depth.min(0) >= min_depth) & (pooled_minor >= min_pooled_minor)
        if keep.any():
            frequencies.append(stacked[:, positions, minor][:, keep] / depth[:, keep])
            depths.append(depth[:, keep])
    af, dp = np.concatenate(frequencies, 1), np.concatenate(depths, 1)
    print(f"{af.shape[1]} sites with pooled minor-allele frequency >= {min_pooled_minor} "
          f"and depth >= {min_depth} in all {len(names)} libraries")
    print("per-library allele-frequency histogram (0-1 in tenths); one diploid animal peaks at 0, 0.5 and 1:")
    for name, row in zip(names, af):
        print(f"  {name:20s} {np.histogram(row, bins=10, range=(0, 1))[0]}")
    rng = np.random.default_rng(seed)
    group = lambda name: name.rsplit("_", 1)[0]
    summary = {"within": [], "between": []}
    for i, j in itertools.combinations(range(len(names)), 2):
        shared = (af[i] * dp[i] + af[j] * dp[j]) / (dp[i] + dp[j])
        noise = np.mean(np.abs(rng.binomial(dp[i].astype(int), shared) / dp[i]
                               - rng.binomial(dp[j].astype(int), shared) / dp[j]))
        kind = "within" if group(names[i]) == group(names[j]) else "between"
        summary[kind].append((np.mean(np.abs(af[i] - af[j])), noise, np.corrcoef(af[i], af[j])[0, 1]))
    print("mean |allele-frequency difference| (observed, same-RNA simulation, correlation):")
    for kind, rows in summary.items():
        if rows:
            rows = np.array(rows)
            print(f"  {kind:8s} pairs={len(rows):3d} observed {rows[:, 0].mean():.3f} "
                  f"({rows[:, 0].min():.3f}-{rows[:, 0].max():.3f}) noise {rows[:, 1].mean():.3f} r {rows[:, 2].mean():.2f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    targets = commands.add_parser("targets")
    targets.add_argument("--quant-dir", required=True)
    targets.add_argument("--tx2gene", required=True)
    targets.add_argument("--output", required=True)
    piled = commands.add_parser("pileup")
    piled.add_argument("--targets", required=True)
    piled.add_argument("--transcripts", required=True)
    piled.add_argument("--output", required=True)
    compared = commands.add_parser("compare")
    compared.add_argument("npz", nargs="+")
    args = parser.parse_args()
    if args.command == "targets":
        select_targets(args.quant_dir, args.tx2gene, args.output)
    elif args.command == "pileup":
        pileup(sys.stdin, args.targets, args.transcripts, args.output)
    else:
        compare(args.npz)
