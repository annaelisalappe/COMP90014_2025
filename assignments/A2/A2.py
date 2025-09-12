import sys
import argparse
import resource
import time
import os
import networkx as nx
import matplotlib.pyplot as plt
from Bio import Seq, Align


def main(indir: str, k_size: int, min_plasmid_len: int, outfile: str) -> None:
    """
    Identifies AMR plasmid sequence given DNA sequencing reads for 5 bacterial organisms ABCDE.
    Organisms A,B,C contain the AMR plasmid, while organisms D,E do not.

    Inputs
        indir:              Path to directory containing sequencing read .fasta files.
        k_size:             Size of K to use (if using kmers).
        min_plasmid_len:    Minimum length of plasmid to consider AMR.
        outfile:            Path to output file where AMR plasmid sequence is written.

    Output
        Single .txt file containing suspected AMR plasmid sequence.
    """
    kmer_counts = build_global_kmer_counts(indir, k_size)
    candidate_kmers, rest_kmers = select_candidate_rest_kmers(kmer_counts, min_diff=100)

    print(f"Candidate k-mers: {len(candidate_kmers)}")
    print(candidate_kmers)
    print(f"Rest k-mers: {len(rest_kmers)}")
    print(rest_kmers)

    candidate_G = build_de_bruijn_graph(candidate_kmers, k_size)
    candidate_contigs = extract_contigs(candidate_G)

    rest_G = build_de_bruijn_graph(rest_kmers, k_size)
    rest_contigs = extract_contigs(rest_G)
    print(f"Rest contigs: {len(rest_contigs)}")
    print(rest_contigs)

    candidate_contigs = [k for k in candidate_contigs if len(k) >= min_plasmid_len]
    print(f"Candidate contigs: {len(candidate_contigs)}")
    print(candidate_contigs)

    result = extract_plasmid_regions(candidate_contigs, rest_contigs, min_plasmid_len)

    for cand, plasmids in result.items():
        print(
            f"Candidate contig: {rotate_sequence(cand)} -> possible plasmids: {rotate_sequence(plasmids)}"
        )


def extract_plasmid_regions(candidate_contigs, rest_contigs, min_plasmid_len):
    """
    Align each candidate contig to each rest contig using global alignment.
    Return plasmid-specific subsequences (uncovered regions) as a dictionary.
    """
    all_possible_plasmids = {}
    aligner = Align.PairwiseAligner(
        match_score=1.0,
        mode="fogsaa",
        open_gap_score=-1.0,
        extend_gap_score=0,
        mismatch_score=-1.0,
    )

    for cand_str in candidate_contigs:
        cand = Seq.Seq(cand_str)
        coverage = [False] * len(cand)  # mark positions covered by any rest contig

        # Align candidate to all rest contigs
        for rest_str in rest_contigs:
            if "CATCC" in rest_str:
                print(rest_str)
            rest = Seq.Seq(rest_str)
            alignments = aligner.align(cand, rest)
            if not alignments:
                print("No alignments found.")
                continue

            best = alignments[0]
            print(best)
            # convert alignment to strings
            cand_algn = str(best[0])
            rest_algn = str(best[1])

            # mark positions covered by this alignment
            idx_cand = 0
            for i in range(len(cand_algn)):
                if cand_algn[i] != "-":  # skip gaps in candidate
                    if rest_algn[i] != "-":
                        coverage[idx_cand] = True
                    idx_cand += 1

        # Extract uncovered regions
        possible_plasmids = []
        plasmid_seq = ""
        for i, is_cov in enumerate(coverage):
            if not is_cov:
                plasmid_seq += cand[i]
            else:
                if len(plasmid_seq) >= min_plasmid_len:
                    possible_plasmids.append(plasmid_seq)
                plasmid_seq = ""
        if len(plasmid_seq) >= min_plasmid_len:
            possible_plasmids.append(plasmid_seq)

        all_possible_plasmids[cand_str] = possible_plasmids

    return all_possible_plasmids


def extract_contigs(G: nx.DiGraph) -> list[str]:
    contigs = []
    visited_edges = set()

    def is_1_in_1_out(node):
        return G.in_degree(node) == 1 and G.out_degree(node) == 1

    # Traverse non 1-in-1-out nodes
    for node in G.nodes:
        if not is_1_in_1_out(node):
            for succ in G.successors(node):
                if (node, succ) not in visited_edges:
                    path = [node, succ]
                    visited_edges.add((node, succ))
                    while is_1_in_1_out(path[-1]):
                        next_node = next(G.successors(path[-1]))
                        visited_edges.add((path[-1], next_node))
                        path.append(next_node)
                    contig = path[0] + "".join(p[-1] for p in path[1:])
                    contigs.append(contig)

    # Handle isolated cycles
    for node in G.nodes:
        if is_1_in_1_out(node):
            succ = next(G.successors(node))
            if (node, succ) not in visited_edges:
                cycle = [node, succ]
                visited_edges.add((node, succ))
                while cycle[-1] != node:
                    next_node = next(G.successors(cycle[-1]))
                    visited_edges.add((cycle[-1], next_node))
                    cycle.append(next_node)
                contig = cycle[0] + "".join(p[-1] for p in cycle[1:])
                contigs.append(contig)

    return contigs


def extract_candidate_kmers_from_mask(kmer_mask):
    """
    Return the set of candidate k-mers present in A,B,C (bits 0,1,2 set)
    and absent in D,E (bits 3,4 unset). Target mask is exactly 0b00111 == 7.
    """
    candidates = [
        k for k, m in kmer_mask.items() if (m & 0b00111) != 0 and (m & 0b11000) == 0
    ]

    return candidates


def get_kmers_from_file(filepath: str, k: int) -> dict[str, int]:
    """
    Reads a single FASTA file and returns a dictionary of k-mer -> count.

    Args:
        filepath: Path to the FASTA file.
        k: Size of k-mers.

    Returns:
        Dictionary mapping each k-mer to the number of times it appears in the file.
    """
    kmers = dict()
    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip()
                for i in range(len(seq) - k + 1):
                    kmer = seq[i : i + k]
                    kmers[kmer] = kmers.get(kmer, 0) + 1
    print(f"{len(kmers)} unique {k}-mers found in {os.path.basename(filepath)}")
    return kmers


def build_global_kmer_counts(indir: str, k: int) -> dict[str, list[int]]:
    """
    Build a dictionary: kmer -> [count_in_A, count_in_B, ..., count_in_E].

    Args:
        indir: Directory containing organism FASTA files.
        k: k-mer size.

    Returns:
        Dictionary mapping k-mer to list of counts for each organism.
    """
    files = [
        os.path.join(indir, "organismA.fasta"),
        os.path.join(indir, "organismB.fasta"),
        os.path.join(indir, "organismC.fasta"),
        os.path.join(indir, "organismD.fasta"),
        os.path.join(indir, "organismE.fasta"),
    ]

    kmer_counts = dict()

    for idx, path in enumerate(files):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing expected file: {path}")

        per_file_kmers = get_kmers_from_file(path, k)

        for km, count in per_file_kmers.items():
            if km not in kmer_counts:
                kmer_counts[km] = [0] * 5
            kmer_counts[km][idx] = count

        per_file_kmers.clear()

    return kmer_counts


def select_candidate_rest_kmers(kmer_counts, min_diff=50):
    """
    Select candidate k-mers that appear much more in A-C than D-E.

    Args:
        kmer_counts: dict of kmer -> [countA, countB, countC, countD, countE]
        min_diff: minimal difference between counts in A-C vs D-E

    Returns:
        candidate_kmers, rest_kmers
    """
    candidate_kmers = []
    rest_kmers = []

    for kmer, counts in kmer_counts.items():
        counts_abc = counts[0:3]
        counts_de = counts[3:5]

        # Candidate: each of A,B,C at least min_diff more than max of D/E
        if all(c >= max(counts_de) + min_diff for c in counts_abc):
            candidate_kmers.append(kmer)
        else:
            rest_kmers.append(kmer)

    print(f"Candidate k-mers: {len(candidate_kmers)}")
    print(f"Rest k-mers: {len(rest_kmers)}")

    return candidate_kmers, rest_kmers


# def build_de_bruijn_graph(kmer_mask: dict, k: int, organism_idx: int) -> dict:
#     bit = 1 << organism_idx
#     G = nx.DiGraph()

#     for km, mask in kmer_mask.items():
#         if mask & bit:  # k-mer belongs to this organism
#             prefix = km[:-1]
#             suffix = km[1:]
#             G.add_edge(prefix, suffix, kmer=km)
#     print(
#         f"Organism {chr(ord('A') + organism_idx)}: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges"
#     )
#     return G


def build_de_bruijn_graph(candidates: list, k: int) -> dict:
    G = nx.DiGraph()

    for candidate in candidates:
        prefix = candidate[:-1]
        suffix = candidate[1:]
        G.add_edge(prefix, suffix, kmer=candidate)
    print(f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def assemble_contigs(G: nx.DiGraph) -> list[str]:
    """
    Extracts all contigs from a de Bruijn graph.

    A contig is a maximal non-branching path.
    """
    contigs = []
    visited = set()

    def is_non_branching(node):
        return G.in_degree(node) == 1 and G.out_degree(node) == 1

    # Pass 1: extract paths from sources or branch points
    for node in G.nodes:
        if G.out_degree(node) > 0 and not is_non_branching(node):
            for succ in G.successors(node):
                contig = [node, succ]
                visited.add(node)
                visited.add(succ)
                while is_non_branching(succ):
                    nxt = next(G.successors(succ))
                    contig.append(nxt)
                    visited.add(nxt)
                    succ = nxt
                # build sequence
                seq = contig[0]
                for n in contig[1:]:
                    seq += n[-1]
                contigs.append(seq)

    # Pass 2: handle isolated cycles of non-branching nodes
    for node in G.nodes:
        if node not in visited and is_non_branching(node):
            contig = [node]
            curr = node
            while True:
                nxt = next(G.successors(curr))
                contig.append(nxt)
                visited.add(curr)
                curr = nxt
                if curr == node:
                    break
            # build sequence
            seq = contig[0]
            for n in contig[1:]:
                seq += n[-1]
            contigs.append(seq)

    # contigs = list(set(contigs))
    contigs.sort(key=len, reverse=True)
    return contigs


def draw_debruijn_graph(G, title="De Bruijn Graph", node_size=800, font_size=8):
    plt.figure(figsize=(12, 8))

    # Use a spring layout for readability (could also try circular/shell)
    pos = nx.spring_layout(G, seed=42, k=0.5)

    # Draw nodes and edges
    nx.draw_networkx_nodes(
        G, pos, node_size=node_size, node_color="lightblue", edgecolors="black"
    )
    nx.draw_networkx_edges(
        G, pos, arrows=True, arrowstyle="->", arrowsize=12, width=1.2, edge_color="gray"
    )

    # Label nodes with (k-1)-mers
    nx.draw_networkx_labels(G, pos, font_size=font_size, font_family="monospace")

    # Optional: label edges with k-mers
    edge_labels = nx.get_edge_attributes(G, "kmer")
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels=edge_labels,
        font_size=font_size - 1,
        font_family="monospace",
    )

    plt.title(title)
    plt.axis("off")
    plt.show()


### DO NOT ALTER ###
def rotate_sequence(seq: str):
    """
    Finds the lexicographically smallest rotation of an input DNA sequence.
    """
    n = len(seq)
    seq_concat = seq + seq
    best_rotation = seq
    for i in range(n):
        rotation = seq_concat[i : i + n]
        if rotation < best_rotation:
            best_rotation = rotation
    return best_rotation


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identifies AMR plasmid sequence given DNA sequencing reads for 5 bacterial organisms ABCDE."
    )
    parser.add_argument(
        "--indir",
        type=str,
        required=True,
        help="Path to directory containing sequencing read .fasta files",
    )
    parser.add_argument(
        "--kmer-size", type=int, required=True, help="Size of K to use (if using kmers)"
    )
    parser.add_argument(
        "--min-plasmid-size",
        type=int,
        required=True,
        help="Minimum length of plasmid to consider AMR",
    )
    parser.add_argument(
        "--outfile",
        type=str,
        required=True,
        help="Path to output file where AMR plasmid sequence is written",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = load_cmdline_args()
    start_time = time.time()
    main(args.indir, args.kmer_size, args.min_plasmid_size, args.outfile)
    stop_time = time.time()
    peak_memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":  # macOS -> bytes
        peak_memory_mb = peak_memory / (1024 * 1024)
    else:
        peak_memory_mb = peak_memory / 1024
    print(f"Time={stop_time - start_time:.2f} seconds")
    print(f"Memory={peak_memory_mb:.2f} MB")
### DO NOT ALTER ###
