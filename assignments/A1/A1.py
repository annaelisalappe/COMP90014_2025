import argparse
import csv
import resource
import sys
import time
import os


def main(database_path: str, oligos_path: str, outfile: str) -> None:
    """
    Identifies which oligos will bind to gene coding sequences (CDS).

    Inputs
        database_path: path to .tsv file containing gene CDS.
        oligos_path: path to .txt file containing oligos to assess.
        outfile: path to .tsv file which summarises results.

    Output
        Single .tsv (tab seperated) file.
        - Column 1: The oligo being assessed.
        - Column 2: True if the oligo matches the WNT4 gene (up to 1 mismatch) else False.
        - Column 3: Comma-separated string of other genes the oligo matches (up to 1 mismatch).
        Example files are available in the ./test_data/expected directory.
    """
    # Load genes into dictionary
    genes = {}
    with open(database_path, "r") as f:
        r = csv.reader(f, delimiter="\t")
        for row in r:
            if len(row) >= 2:
                gene_name, gene_seq = row[0], row[-1]
                genes[gene_name] = gene_seq.upper()

    # Load oligos into list
    oligos = []
    with open(oligos_path, "r") as f:
        for line in f:
            oligo = line.strip().upper()
            if oligo:
                oligos.append(oligo)

    # Main logic
    oligo_length = len(oligos[0])  # All oligos have the same length
    index = build_kmer_index(genes=genes, k=oligo_length)
    results = lookup_variants_in_index(oligos=oligos, index=index)

    # Write to output file
    if not outfile.lower().endswith(".tsv"):
        outfile = os.path.splitext(outfile)[0] + ".tsv"
        print("INFO: Converting output file to tsv.")

    with open(outfile, "w", newline="") as tsv_file:
        w = csv.writer(tsv_file, delimiter="\t")
        w.writerow(["oligo", "WNT4", "off_target"])
        w.writerows(results)


def build_kmer_index(genes: dict[str, str], k: int) -> dict[str, list[str]]:
    """
    Function to build an index of all k-mers (k should be to the
    length of the oligos) that appear in all genes. Called at the start of the program.
    The k-mer and gene are stored in a dictionary, which by default in python is
    implemented in constant-lookup time through hashing.
    """
    index = {}
    for gene_name, sequence in genes.items():
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i : i + k]
            if kmer not in index:
                index[kmer] = []
            index[kmer].append(gene_name)
    return index


def generate_all_variants_of_seq(seq: str) -> set[str]:
    """
    Return all sequences with zero or one mismatch (Hamming distance smaller or equal to 1).
    """
    bases = ["A", "C", "G", "T"]
    variants = set()
    variants.add(seq)

    for i, _ in enumerate(seq):
        for b in bases:
            variant = seq[:i] + b + seq[i + 1 :]
            variants.add(variant)

    return variants


def lookup_variants_in_index(
    oligos: list[str], index: dict[str, list[str]]
) -> list[list[str, str, str]]:
    """
    Function that loops over all oligos and generate all their variants and looks each variant
    up in index. The complexity of this lookup function is n_genes * (oligo_length * 4) * 2.
    Returns a list containing the oligo, True/ False indicating if oligo binds to WNT4, and
    names of any off-target genes.
    """
    results = []
    for oligo in oligos:
        binds_wnt4 = False
        off_targets = []
        oligo_variants = generate_all_variants_of_seq(seq=oligo)

        # Looping over all variants adds work but only 4^k, where k is length of oligos
        for variant in oligo_variants:
            variant_genes = index.get(variant)

            if variant_genes is not None:
                if "WNT4" in variant_genes:
                    binds_wnt4 = True
                    # Create copy of variant genes to avoid changing the original in index;
                    # remove WNT4 gene to add variant_genes to off_target list
                    variant_genes = list.copy(variant_genes)
                    variant_genes.remove("WNT4")

                off_targets.extend(variant_genes)

        # Add only unique off-target genes, and format result to fit output specifications
        unique_off_targets = list(set(off_targets))
        results.append([oligo, str(binds_wnt4), ",".join(unique_off_targets)])

    return results


### DO NOT ALTER ###
def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identifies the best oligo for a gene which minimises off-target mRNA \
        transcript binding."
    )
    parser.add_argument(
        "--database",
        type=str,
        required=True,
        help="Coding sequence database path (.tsv)",
    )
    parser.add_argument("--oligos", type=str, required=True, help="Oligos path (.txt)")
    parser.add_argument(
        "--outfile", type=str, required=True, help="Results path (.tsv)"
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = load_cmdline_args()
    start_time = time.time()
    main(args.database, args.oligos, args.outfile)
    stop_time = time.time()
    peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":  # macOS -> bytes
        peak_memory_mb = peak_memory / (1024 * 1024)
    else:
        peak_memory_mb = peak_memory / 1024
    print(f"Time={stop_time - start_time:.2f} seconds")
    print(f"Memory={peak_memory_mb:.2f} MB")
### DO NOT ALTER ###
