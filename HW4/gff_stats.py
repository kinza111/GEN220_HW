#!/usr/bin/env python3
import argparse
import os

def parse_attributes(attr_str):
    """Parse GFF3 attributes column into a dict."""
    d = {}
    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def mean(xs):
    return sum(xs) / len(xs) if xs else None

def summarize_lengths(xs):
    """Return (min, mean, max) for a list or (None,None,None) if empty."""
    if not xs:
        return (None, None, None)
    return (min(xs), mean(xs), max(xs))

def merge_intervals(intervals):
    """Merge overlapping/adjacent intervals. intervals: list of (start,end) 1-based inclusive."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:  
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged

def total_merged_bp(intervals):
    """Total bp covered by merged intervals (1-based inclusive)."""
    merged = merge_intervals(intervals)
    return sum(e - s + 1 for s, e in merged)

def read_gff3(path):
    """
    Read GFF3 and collect:
    - genes per contig
    - transcripts and their parent gene
    - exon intervals per transcript (or gene if Parent points to gene)
    - CDS intervals per contig (for coding fraction)
    """
    genes_by_seq = {}          # seqid -> list of (start,end,gene_id)
    transcript_to_gene = {}    # transcript_id -> gene_id
    exons_by_transcript = {}   # transcript_id -> list of (start,end)
    exons_by_gene = {}         # gene_id -> list of (start,end)
    cds_by_seq = {}            # seqid -> list of (start,end)

    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = fields
            start = int(start)
            end = int(end)
            a = parse_attributes(attrs)

            if ftype == "gene":
                gid = a.get("ID", None)
                if gid is None:
                    gid = f"gene_{seqid}_{start}_{end}"
                genes_by_seq.setdefault(seqid, []).append((start, end, gid))

            elif ftype in ("mRNA", "transcript"):
                tid = a.get("ID", None)
                parent = a.get("Parent", None)
                if tid and parent:
                    transcript_to_gene[tid] = parent

            elif ftype == "exon":
                parent = a.get("Parent", None)
                if not parent:
                    continue
                if parent in transcript_to_gene:
                    exons_by_transcript.setdefault(parent, []).append((start, end))
                else:
                    exons_by_gene.setdefault(parent, []).append((start, end))

            elif ftype == "CDS":
                cds_by_seq.setdefault(seqid, []).append((start, end))
                parent = a.get("Parent", None)
                if parent:
                    exons_by_transcript.setdefault(parent, [])  
    return genes_by_seq, transcript_to_gene, exons_by_transcript, exons_by_gene, cds_by_seq

def compute_stats(gff_path):
    genes_by_seq, transcript_to_gene, exons_by_transcript, exons_by_gene, cds_by_seq = read_gff3(gff_path)
    gene_lengths = []
    intergenic = []
    annotated_span_bp = 0
    n_genes = 0

    for seqid, genes in genes_by_seq.items():
        genes_sorted = sorted(genes, key=lambda x: (x[0], x[1]))
        n_genes += len(genes_sorted)

        # annotated span for this contig based on genes (twist metric)
        min_s = min(g[0] for g in genes_sorted)
        max_e = max(g[1] for g in genes_sorted)
        annotated_span_bp += (max_e - min_s + 1)

        # gene length
        for s, e, gid in genes_sorted:
            gene_lengths.append(e - s + 1)

        # intergenic distances (gap between neighboring genes; overlaps -> 0)
        for (s1, e1, _), (s2, e2, _) in zip(genes_sorted, genes_sorted[1:]):
            gap = s2 - e1 - 1
            intergenic.append(gap if gap > 0 else 0)

    # Exons per gene and exon lengths / intron lengths
    exon_lengths = []
    intron_lengths = []
    exons_per_gene = {}  # gene_id -> count

    # 1) Use transcript-based exons if present
    for tid, exons in exons_by_transcript.items():
        gid = transcript_to_gene.get(tid, None)
        if not gid:
            continue
        if not exons:
            continue
        exons_sorted = sorted(exons, key=lambda x: x[0])
        exons_per_gene[gid] = exons_per_gene.get(gid, 0) + len(exons_sorted)

        # exon sizes
        for s, e in exons_sorted:
            exon_lengths.append(e - s + 1)

        # introns = gaps between consecutive exons
        for (s1, e1), (s2, e2) in zip(exons_sorted, exons_sorted[1:]):
            gap = s2 - e1 - 1
            if gap > 0:
                intron_lengths.append(gap)
                if not exon_lengths and exons_by_gene:
                    for gid, exons in exons_by_gene.items():
                        exons_sorted = sorted(exons, key=lambda x: x[0])
                        exons_per_gene[gid] = exons_per_gene.get(gid, 0) + len(exons_sorted)
                        for s, e in exons_sorted:
                            exon_lengths.append(e - s + 1)
                            for (s1, e1), (s2, e2) in zip(exons_sorted, exons_sorted[1:]):
                                gap = s2 - e1 - 1
                                if gap > 0:
                                    intron_lengths.append(gap)

    # 3) If STILL no exons, use CDS as a proxy for "exon-like segments"
    # (common in bacteria: exon features often absent; CDS segments represent coding regions)
    if not exon_lengths:
        # we cannot do introns meaningfully here
        # we also can't do exons-per-gene reliably without transcript->gene mapping
        # but we can provide CDS length stats as exon proxy.
        for seqid, cds_intervals in cds_by_seq.items():
            for s, e in cds_intervals:
                exon_lengths.append(e - s + 1)

    # Mean exons per gene (only for genes where we observed exon counts)
    exons_per_gene_values = list(exons_per_gene.values())
    mean_exons_per_gene = mean(exons_per_gene_values)

    # Coding fraction (twist): total merged CDS bp / annotated span
    total_cds_bp = 0
    for seqid, cds_intervals in cds_by_seq.items():
        total_cds_bp += total_merged_bp(cds_intervals)

    coding_fraction = (total_cds_bp / annotated_span_bp) if annotated_span_bp > 0 else None

    genes_per_mb = (n_genes / (annotated_span_bp / 1_000_000)) if annotated_span_bp > 0 else None

    # Summaries
    gene_min, gene_mean, gene_max = summarize_lengths(gene_lengths)
    exon_min, exon_mean, exon_max = summarize_lengths(exon_lengths)
    intr_min, intr_mean, intr_max = summarize_lengths(intron_lengths)
    inter_min, inter_mean, inter_max = summarize_lengths(intergenic)

    return {
        "file": os.path.basename(gff_path),
        "n_genes": n_genes,
        "annotated_span_bp": annotated_span_bp,
        "genes_per_mb": genes_per_mb,
        "gene_len_min": gene_min,
        "gene_len_mean": gene_mean,
        "gene_len_max": gene_max,
        "exon_len_min": exon_min,
        "exon_len_mean": exon_mean,
        "exon_len_max": exon_max,
        "intron_len_min": intr_min,
        "intron_len_mean": intr_mean,
        "intron_len_max": intr_max,
        "mean_exons_per_gene": mean_exons_per_gene,
        "intergenic_min": inter_min,
        "intergenic_mean": inter_mean,
        "intergenic_max": inter_max,
        "coding_fraction": coding_fraction,
    }

def fmt(x):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="Input GFF3 file")
    ap.add_argument("--out", required=True, help="Output TSV file")
    args = ap.parse_args()

    stats = compute_stats(args.gff)

    cols = [
        "file",
        "n_genes",
        "annotated_span_bp",
        "genes_per_mb",
        "gene_len_min", "gene_len_mean", "gene_len_max",
        "exon_len_min", "exon_len_mean", "exon_len_max",
        "intron_len_min", "intron_len_mean", "intron_len_max",
        "mean_exons_per_gene",
        "intergenic_min", "intergenic_mean", "intergenic_max",
        "coding_fraction",
    ]

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(fmt(stats[c]) for c in cols) + "\n")

if __name__ == "__main__":
    main()

