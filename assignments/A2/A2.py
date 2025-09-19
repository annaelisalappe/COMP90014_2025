import argparse
import os
import resource
import sys
import time
import networkx as nx


def main(indir: str, k_size: int, min_plasmid_len: int, outfile: str) -> None:
    """
    Identifies AMR plasmid sequence given DNA sequencing reads for 5 bacterial organisms ABCDE.
    Organisms A,B,C contain the AMR plasmid, while organisms D,E do not.

    Args:
        indir:              Path to directory containing sequencing read .fasta files.
        k_size:             Size of K to use (if using kmers).
        min_plasmid_len:    Minimum length of plasmid to consider AMR.
        outfile:            Path to output file where AMR plasmid sequence is written.

    Outputs:
        Single .txt file containing suspected AMR plasmid sequence.
    """
    candidate_plasmids = []
    kmer_counts = True

    while not candidate_plasmids and kmer_counts:
        print(f"Trying to find plasmid using k={k_size}...")

        kmer_counts = build_global_kmer_counts(
            indir, k_size
        )  # will be none if kmer length longer than read length
        candidate_plasmids = try_to_find_plasmid_seq(
            kmer_counts=kmer_counts, k_size=k_size, min_plasmid_len=min_plasmid_len
        )  # will be empty list if no candidate plasmids are found
        if not candidate_plasmids:
            print(
                f"Failed to find circular contigs with sufficient length using k={k_size}."
            )
        k_size += 1

    if not candidate_plasmids:
        return None

    rotated_candidate_contigs = [rotate_sequence(k) for k in candidate_plasmids]

    with open(outfile, "w") as f:
        for seq in rotated_candidate_contigs:
            f.write(seq + "\n")

    print(f"Plasmid discovery ran successfully. Wrote output to {outfile}.")


def try_to_find_plasmid_seq(kmer_counts: dict, k_size: int, min_plasmid_len: int):
    """
    Main pipeline for the plasmid sequence discovery.

    Args:
        kmer_counts (dict): Dictionary mapping k-mer -> list of counts per organism.
        k_size (int): The k-mer size used to build the graph.
        min_plasmid_len (int): Minimum length of candidate plasmid sequences to keep.

    Returns:
        list[str]: List of trimmed candidate plasmid sequences with length ≥ `min_plasmid_len`.
    """
    candidate_kmers = select_candidate_kmers(kmer_counts, diff_factor=1.5)

    candidate_G = build_de_bruijn_graph(candidate_kmers, k_size)
    candidate_G = collapse_linear_paths(candidate_G)
    candidate_G = remove_tips(
        candidate_G, min_length=2 * k_size, min_coverage=2, k=k_size
    )
    candidate_contigs = extract_contigs_from_graph(candidate_G)

    long_candidates_before_trimming = [
        c for c in candidate_contigs if len(c) >= min_plasmid_len
    ]

    trimmed_candidates = validate_and_trim_candidates(
        long_candidates_before_trimming,
        min_overlap=k_size // 2,
        min_plasmid_len=min_plasmid_len,
    )

    long_candidates_after_trimming = [
        c for c in trimmed_candidates if len(c) >= min_plasmid_len
    ]

    return long_candidates_after_trimming


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
            if line.startswith(">"):
                continue
            seq = line.strip()
            for i in range(len(seq) - k + 1):
                kmer = seq[i : i + k]
                kmers[kmer] = kmers.get(kmer, 0) + 1
    return kmers


def build_global_kmer_counts(indir: str, k: int) -> dict[str, list[int]]:
    """
    Builds a dictionary: kmer -> [count_in_A, count_in_B, ..., count_in_E].

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

    for id, path in enumerate(files):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing expected file: {path}")

        per_file_kmers = get_kmers_from_file(path, k)

        for km, count in per_file_kmers.items():
            if km not in kmer_counts:  # first time seeing this k-mer
                kmer_counts[km] = [0] * 5
            kmer_counts[km][id] = count

        per_file_kmers.clear()

    return kmer_counts


def select_candidate_kmers(kmer_counts, diff_factor=2) -> dict:
    """
    Selects candidate k-mers that appear (at least diff_factor) more in A-C than D-E.
    The calculation is conservative: at least one of A,B,C must be >= diff_factor * max(D,E).

    Args:
        kmer_counts: dict of kmer -> [countA, countB, countC, countD, countE]
        diff_factor: minimal difference factor between counts in A-C vs D-E

    Returns:
        candidate_kmers, rest_kmers
    """
    candidate_kmers = {}

    for kmer, counts in kmer_counts.items():
        counts_abc = counts[0:3]
        counts_de = counts[3:5]

        # Assume k-mer must appear in at least one of A,B,C
        if any(c == 0 for c in counts_abc):
            continue

        # Candidate: At least one of A,B,C at least diff_factor times more than max of D/E
        if any(c >= max(counts_de) * diff_factor for c in counts_abc):
            candidate_kmers[kmer] = sum(counts_abc) / 3  # average count in A,B,C

    return candidate_kmers


def build_de_bruijn_graph(kmer_counts: dict, k: int) -> dict:
    """
    Builds a De Bruijn graph from k-mers and their counts.
    Nodes are (k-1)-mers, edges are k-mers with weights as counts.
    Args:
        kmer_counts: dict of kmer -> count
        k: k-mer size
    Returns:
        A NetworkX directed graph representing the De Bruijn graph.
    """
    G = nx.DiGraph()
    for km, counts in kmer_counts.items():
        if counts > 0:  # k-mer appears in this organism
            prefix = km[:-1]
            suffix = km[1:]
            # If edge already exists, increment its weight
            if G.has_edge(prefix, suffix):
                G[prefix][suffix]["cov"] += counts
            else:
                G.add_edge(prefix, suffix, seq=km, cov=counts)
    return G


def get_de_brujn_start_nodes(G: nx.DiGraph):
    """
    Identifies starting nodes in a De Bruijn graph.

    A start node is defined as in tutorial 7:
    - A node with in-degree 0 (source).
    - A node with in-degree > 1 (branch point).
    - Or a child of a node with out-degree > 1.

    Args:
        G (nx.DiGraph): De Bruijn graph.

    Returns:
        list[str]: List of starting nodes.
    """
    start_nodes = []
    for node in G.nodes():
        if G.in_degree(node) == 0 or G.in_degree(node) > 1:
            start_nodes.append(node)
        elif G.out_degree(node) > 1:
            child_nodes = list(G.successors(node))
            start_nodes += child_nodes

    return start_nodes


def traverse_linear_path_in_de_brujn(G: nx.DiGraph, start, visited: set):
    """
    Traverses a non-branching path in a De Bruijn graph starting from a given node.

    Args:
        G (nx.DiGraph): De Bruijn graph
        start (str): Node from which to start the traversal.
        visited (set): Set of nodes that have already been visited; updated in-place.

    Returns:
        tuple:
            contig (str): Assembled contig sequence along the path.
            current_node (str): Last node reached before stopping.
            avg_count (float): Average k-mer coverage along the path.
            visited (set): Updated set of visited nodes including those traversed in this path.
    """
    contig = start
    current_node = start
    counts = 0
    length = 0

    while True:
        out_edges = list(G.out_edges(current_node, data=True))
        if len(out_edges) == 0 or any(
            succ in visited for succ in G.successors(current_node)
        ):
            break
        # assuming only one outgoing edge
        _, next_node, edge_data = out_edges[0]
        contig += edge_data["seq"][-1]
        counts += edge_data["cov"]
        length += 1
        current_node = next_node
        visited.add(current_node)

    avg_count = counts / length if length > 0 else counts

    return contig, current_node, avg_count, visited


def collapse_linear_paths(G: nx.DiGraph) -> nx.DiGraph:
    """
    Collapses unambiguous linear paths (unitigs) into single nodes.
    Each edge must have:
        - 'seq': the k-mer represented by that edge
        - 'cov': coverage (k-mer count)
    Returns a new graph with unitigs collapsed into nodes.

    Args:
        G: Input De Bruijn graph (networkx DiGraph)
    Returns:
        A new DiGraph with linear paths collapsed into single nodes.
    """
    H = nx.DiGraph()

    # Identify start nodes
    start_nodes = get_de_brujn_start_nodes(G)
    visited = set(start_nodes)

    for start in start_nodes:
        contig, current_node, avg_count, visited = traverse_linear_path_in_de_brujn(
            G, start, visited
        )
        prefix = contig[:-1]
        suffix = current_node[1:]

        H.add_edge(prefix, suffix, seq=current_node, cov=avg_count)

    for node in G.nodes:
        if node in visited:
            continue

        contig, current_node, avg_count, visited = traverse_linear_path_in_de_brujn(
            G, node, visited
        )
        prefix = contig[:-1]
        suffix = current_node[1:]

        H.add_edge(prefix, suffix, seq=current_node, cov=avg_count)
    return H


def remove_tips(G: nx.DiGraph, min_length: int, min_coverage: int, k: int):
    """
    Removes tips (dead-end paths) from the De Bruijn graph.
    A tip is defined as a path that starts or ends at a node with in-degree or
    out-degree of 0, and has a total length less than min_length and average
    coverage less than min_coverage.

    Args:
        G: Input De Bruijn graph (networkx DiGraph)
        min_length: Minimum length of tip to retain
        min_coverage: Minimum average coverage of tip to retain
        k: k-mer size (for calculating length)

    Returns:
        The graph with tips removed.
    """
    removed = True
    while removed:
        removed = False
        for node in list(G.nodes):
            if G.in_degree(node) != 0 and G.out_degree(node) != 0:  # Look at tips
                continue

            path_nodes = [node]
            edges = []
            current = node

            # Traverse forward (source tip)
            while G.out_degree(current) == 1 and G.in_degree(current) <= 1:
                _, nxt, data = list(G.out_edges(current, data=True))[0]
                edges.append(data)
                current = nxt
                path_nodes.append(current)

            # Traverse backward (sink tip)
            while G.in_degree(node) == 1 and G.out_degree(node) == 0:
                prev, _, data = list(G.in_edges(node, data=True))[0]
                edges.append(data)
                node = prev
                path_nodes.append(node)

            if not edges:
                continue

            # Total number of characters along tip
            total_contig_len = sum(len(e["seq"]) for e in edges) - (k - 1) * (
                len(edges) - 1
            )
            avg_cov = sum(d["cov"] for d in edges) / len(edges)
            if total_contig_len < min_length and avg_cov < min_coverage:
                G.remove_nodes_from(path_nodes[:-1])  # keep the junction node
                removed = True

    return remove_de_brujn_islands(G)


def remove_de_brujn_islands(G: nx.DiGraph):
    """
    Removes isolated nodes after removing tips.

    Args:
        G: Input De Bruijn graph

    Returns:
        G: Input De Bruijn graph with isolated nodes removed.
    """
    for nodes in list(G.nodes):
        if G.in_degree(nodes) == 0 and G.out_degree(nodes) == 0:
            G.remove_node(nodes)

    return G


def extract_contigs_from_graph(G: nx.DiGraph) -> list[str]:
    """
    Extracts maximal contigs (non-branching paths) from a De Bruijn graph.
    Assumes the graph has already been collapsed into unitigs.

    Args:
        G: Input De Bruijn graph (networkx DiGraph) with unitigs as nodes

    Returns:
        Dictionary mapping contig sequences to their average k-mer coverage.
    """

    # Identify start nodes
    start_nodes = get_de_brujn_start_nodes(G)
    visited = set(start_nodes)

    contigs = {}
    for start in start_nodes:
        contig, _, avg_count, visited = traverse_linear_path_in_de_brujn(
            G, start, visited
        )

        contigs[contig] = avg_count

    # Handle cycles where no start node is present
    for node in G.nodes:
        if node in visited:
            continue

        # Handle self-loop where node sequence is the full sequence
        if G.has_edge(node, node):
            contig = G[node][node]["seq"]
            avg_count = G[node][node]["cov"]
            contigs[contig] = avg_count
            visited.add(node)
            continue

        contig, _, avg_count, visited = traverse_linear_path_in_de_brujn(
            G, node, visited
        )

        contigs[contig] = avg_count

    return contigs


def validate_and_trim_candidates(
    sequences: list[str], min_overlap: int, min_plasmid_len: int
) -> str:
    """
    Checks for a list of sequences if they are circular and trims them if necessary.
    If a sequence is not circular, it is dismissed (plasmid sequence must be circular).

    Args:
        sequences: The DNA sequences to validate and trim.
        min_overlap: The minimum required length for the overlap to be valid.
        min_plasmid_len: The minimum length for a sequence to be considered a valid plasmid.

    Returns:
        A list of valid, trimmed candidate sequences, or None if no valid sequences remain.
    """
    final_seqs = []
    for sequence in sequences:
        n = len(sequence)
        overlap_length = 0
        for i in range(n // 2, min_overlap, -1):  # longest overlap first
            overlap = sequence[-i:]

            if sequence.startswith(overlap):
                overlap_length = i
                break  # found the longest overlap

        if not overlap_length >= min_overlap:
            break  # no overlap found; dismiss this sequence

        print(f"Found circular contig (overlap={i} bp). Trimming...")
        final_seqs.append(sequence[:-overlap_length])

    return final_seqs


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
