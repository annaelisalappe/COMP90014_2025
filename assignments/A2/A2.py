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
    candidate_kmers, rest_kmers = select_candidate_rest_kmers(
        kmer_counts, diff_factor=10000000
    )

    print(f"Candidate k-mers: {len(candidate_kmers)}")
    print(f"Rest k-mers: {len(rest_kmers)}")

    candidate_G = build_de_bruijn_graph(candidate_kmers, k_size)
    candidate_G = collapse_linear_paths(candidate_G)
    candidate_G = remove_tips(
        candidate_G, min_length=2 * k_size, min_coverage=2, k=k_size
    )

    draw_debruijn_graph(
        candidate_G,
        title="Candidate De Bruijn Graph (k-mers in A,B,C but not D,E)",
        node_size=800,
        font_size=8,
    )
    candidate_contigs = extract_contigs_from_graph(candidate_G)
    rotated_candidate_contigs = [
        rotate_sequence(k)
        for k in candidate_contigs.keys()
        if len(k) >= min_plasmid_len
    ]

    print(f"Number of nodes in candidate graph: {candidate_G.number_of_nodes()}")
    print(f"Number of edges in candidate graph: {candidate_G.number_of_edges()}")
    print(f"Edges: {list(candidate_G.edges(data=True))}")
    print(f"Rotated candidate contigs: {rotated_candidate_contigs}")

    # Bar plot of contig coverages
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(candidate_contigs)), candidate_contigs.values())
    plt.xlabel("Contig Index")
    plt.ylabel("Average k-mer Coverage")
    plt.title("Contig Coverage Distribution")
    plt.show()

    # rest_G = build_de_bruijn_graph(rest_kmers, k_size)
    # rest_G = collapse_linear_paths(rest_G)
    # # rest_G = remove_tips(rest_G, min_length=k_size, min_coverage=1, k=k_size)
    # draw_debruijn_graph(
    #     rest_G,
    #     title="Rest De Bruijn Graph (k-mers in D,E or not in A,B,C)",
    #     node_size=800,
    #     font_size=8,
    # )
    # rest_contigs = extract_contigs(rest_G)
    # print(f"Rest contigs: {len(rest_contigs)}")
    # print(rest_contigs)

    # candidate_contigs = [k for k in candidate_contigs if len(k) >= min_plasmid_len]
    # print(f"Candidate contigs: {len(candidate_contigs)}")
    # print(candidate_contigs)

    # result = extract_plasmid_regions(candidate_contigs, rest_contigs, min_plasmid_len)

    # for cand, plasmids in result.items():
    #     print(
    #         f"Candidate contig: {rotate_sequence(cand)} -> possible plasmids: {rotate_sequence(plasmids)}"
    #     )


def extract_contigs_from_graph(G: nx.DiGraph) -> list[str]:
    """
    Extracts maximal contigs (non-branching paths) from a unitig-merged De Bruijn graph.
    Assumes nodes store full unitig sequences as strings.
    """
    visited = set()

    # Identify start nodes
    start_nodes = []
    for node in G.nodes():
        if G.in_degree(node) == 0 or G.in_degree(node) > 1:
            start_nodes.append(node)
        elif G.out_degree(node) > 1:
            child_nodes = list(G.successors(node))
            start_nodes += child_nodes

    contigs = {}
    for start in start_nodes:
        print(f"Starting new contig from node: {start}")
        contig = start
        current_node = start
        counts = 0
        length = 0
        while True:
            out_edges = list(G.out_edges(current_node, data=True))
            if len(out_edges) == 0 or any(
                succ in start_nodes for succ in G.successors(current_node)
            ):
                print(f"Ending contig at node: {current_node}")
                break
            # assuming only one outgoing edge
            _, next_node, edge_data = out_edges[0]
            contig += edge_data["seq"][-1]
            counts += edge_data["cov"]
            length += 1
            current_node = next_node
            visited.add(current_node)

        avg_count = counts / length if length > 0 else counts

        contigs[contig] = avg_count

    for node in G.nodes:
        if node not in visited and node not in start_nodes:
            if G.has_edge(node, node):
                # Handle self-loop
                contig = G[node][node]["seq"]
                avg_count = G[node][node]["cov"]
                contigs[contig] = avg_count
                print(f"Self-loop contig: {contig} with coverage {avg_count}")
                visited.add(node)
                continue
            contig = node
            current_node = node
            counts = 0
            length = 0
            while True:
                out_edges = list(G.out_edges(current_node, data=True))
                if len(out_edges) == 0 or any(
                    succ in visited for succ in G.successors(current_node)
                ):
                    break
                _, next_node, edge_data = out_edges[0]
                contig += edge_data["seq"][-1]
                counts += edge_data["cov"]
                length += 1
                visited.add(current_node)
                current_node = next_node
            avg_count = counts / length if length > 0 else counts

            contigs[contig] = avg_count

    return contigs


# def extract_contigs(G: nx.DiGraph) -> list[str]:
#     contigs = []
#     visited_edges = set()

#     def is_1_in_1_out(node):
#         return G.in_degree(node) == 1 and G.out_degree(node) == 1

#     # Traverse non 1-in-1-out nodes
#     for node in G.nodes:
#         if not is_1_in_1_out(node):
#             for succ in G.successors(node):
#                 if (node, succ) not in visited_edges:
#                     path = [node, succ]
#                     visited_edges.add((node, succ))
#                     while is_1_in_1_out(path[-1]):
#                         next_node = next(G.successors(path[-1]))
#                         visited_edges.add((path[-1], next_node))
#                         path.append(next_node)
#                     contig = path[0] + "".join(p[-1] for p in path[1:])
#                     contigs.append(contig)

#     # Handle isolated cycles
#     for node in G.nodes:
#         if is_1_in_1_out(node):
#             succ = next(G.successors(node))
#             if (node, succ) not in visited_edges:
#                 cycle = [node, succ]
#                 visited_edges.add((node, succ))
#                 while cycle[-1] != node:
#                     next_node = next(G.successors(cycle[-1]))
#                     visited_edges.add((cycle[-1], next_node))
#                     cycle.append(next_node)
#                 contig = cycle[0] + "".join(p[-1] for p in cycle[1:])
#                 contigs.append(contig)

#     return contigs


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


def select_candidate_rest_kmers(kmer_counts, diff_factor=2):
    """
    Select candidate k-mers that appear much more in A-C than D-E.

    Args:
        kmer_counts: dict of kmer -> [countA, countB, countC, countD, countE]
        min_diff: minimal difference between counts in A-C vs D-E

    Returns:
        candidate_kmers, rest_kmers
    """
    candidate_kmers = {}
    rest_kmers = {}

    for kmer, counts in kmer_counts.items():
        counts_abc = counts[0:3]
        counts_de = counts[3:5]

        avg_abc = sum(counts_abc) / 3
        avg_de = sum(counts_de) / 2

        if any(c == 0 for c in counts_abc):
            # print("  -> Skipping k-mer as it is missing in one of A,B,C")
            rest_kmers[kmer] = avg_de
            continue
        # Candidate: each of A,B,C at least min_diff more than min of D/E
        if any(c >= max(counts_de) * diff_factor for c in counts_abc):
            candidate_kmers[kmer] = avg_abc
        else:
            rest_kmers[kmer] = avg_de

    return candidate_kmers, rest_kmers


def build_de_bruijn_graph(kmer_counts: dict, k: int) -> dict:
    G = nx.DiGraph()
    for km, counts in kmer_counts.items():
        if counts > 0:  # k-mer appears in this organism
            prefix = km[:-1]
            suffix = km[1:]
            # If edge already exists, increment its weight
            if G.has_edge(prefix, suffix):
                G[prefix][suffix]["cov"] += counts
            else:
                G.add_edge(prefix, suffix, kmer=km, cov=counts)
    return G


def collapse_linear_paths(G) -> nx.DiGraph:
    """
    Collapse unambiguous linear paths (unitigs) into single nodes.
    Each edge must have:
        - 'seq': the k-mer or sequence represented by that edge
        - 'count': coverage (k-mer count)
    """
    H = nx.DiGraph()
    visited = set()

    start_nodes = []
    for node in G.nodes():
        if G.in_degree(node) == 0 or G.in_degree(node) > 1:
            start_nodes.append(node)
        if G.out_degree(node) > 1:
            child_nodes = list(G.successors(node))
            start_nodes += child_nodes

    for start in start_nodes:
        contig = start
        current_node = start
        counts = 0
        length = 0
        while True:
            out_edges = list(G.out_edges(current_node, data=True))
            if len(out_edges) == 0 or any(
                succ in start_nodes for succ in G.successors(current_node)
            ):
                break
            # assuming only one outgoing edge
            _, next_node, edge_data = out_edges[0]
            contig += edge_data["kmer"][-1]
            counts += edge_data["cov"]
            length += 1
            current_node = next_node
            visited.add(current_node)
        prefix = contig[:-1]
        suffix = current_node[1:]
        # if length == 0:
        #     in_edges = list(G.in_edges(current_node, data=True))

        #     if len(in_edges) == 0:

        #     if len(list(G.in_edges(current_node, data=True))) > 0:
        #         print("Warning: length 0 but has mutliple in-edges")
        #     avg_count = list(G.in_edges(current_node, data=True))[0][2]["cov"]
        # else:
        avg_count = counts / length if length > 0 else counts
        H.add_edge(prefix, suffix, seq=current_node, cov=avg_count)

    for node in G.nodes:
        if node not in visited and node not in start_nodes:
            start = node
            contig = start
            current_node = start
            counts = 0
            length = 0
            while True:
                visited.add(current_node)
                out_edges = list(G.out_edges(current_node, data=True))
                if len(out_edges) == 0 or any(
                    succ in visited for succ in G.successors(current_node)
                ):
                    break
                _, next_node, edge_data = out_edges[0]
                contig += edge_data["kmer"][-1]
                counts += edge_data["cov"]
                length += 1
                current_node = next_node
            prefix = start[:-1]
            suffix = current_node[1:]
            avg_count = counts / (length + 1) if length > 0 else counts
            H.add_edge(prefix, suffix, seq=contig, cov=avg_count)
    return H


def remove_tips(G, min_length, min_coverage, k):
    removed = True
    while removed:
        removed = False
        for node in list(G.nodes):
            # Look at tips (sources or sinks)
            if G.in_degree(node) == 0 or G.out_degree(node) == 0:
                path_nodes = [node]
                edges = []
                current = node

                # Traverse forward if it's a source tip
                while G.out_degree(current) == 1 and G.in_degree(current) <= 1:
                    _, nxt, data = list(G.out_edges(current, data=True))[0]
                    edges.append(data)
                    current = nxt
                    path_nodes.append(current)

                # Traverse backward if it's a sink tip
                while G.in_degree(node) == 1 and G.out_degree(node) == 0:
                    prev, _, data = list(G.in_edges(node, data=True))[0]
                    edges.append(data)
                    node = prev
                    path_nodes.append(node)

                if edges:
                    # total number of characters along tip
                    print(edges)
                    total_contig_len = sum(len(e["seq"]) for e in edges) - (k - 1) * (
                        len(edges) - 1
                    )
                    avg_cov = sum(d["cov"] for d in edges) / len(edges)
                    if total_contig_len < min_length and avg_cov < min_coverage:
                        print(
                            f"Removing tip path: {path_nodes} with length {total_contig_len} and avg coverage {avg_cov:.2f}"
                        )
                        G.remove_nodes_from(path_nodes[:-1])  # keep the junction node
                        removed = True

    for nodes in list(G.nodes):
        if G.in_degree(nodes) == 0 and G.out_degree(nodes) == 0:
            G.remove_node(nodes)

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

    # Pass 1: extract paths from   or branch points
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

    pos = nx.spring_layout(G, seed=42, k=0.5)

    # Draw nodes and edges
    nx.draw_networkx_nodes(
        G, pos, node_size=node_size, node_color="lightblue", edgecolors="black"
    )
    nx.draw_networkx_edges(
        G, pos, arrows=True, arrowstyle="->", arrowsize=12, width=1.2, edge_color="gray"
    )

    # Node labels (k-1-mers)
    nx.draw_networkx_labels(G, pos, font_size=font_size, font_family="monospace")

    # Build edge labels dynamically
    edge_labels = {}
    for u, v, d in G.edges(data=True):
        # Gather all known keys
        parts = []
        if "kmer" in d:
            parts.append(d["kmer"])
        if "seq" in d:
            parts.append(f"seq:{d['seq']}")
        if "cov" in d:
            parts.append(f"cov:{d['cov']}")
        label = " | ".join(parts) if parts else "edge"
        edge_labels[(u, v)] = label

    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels=edge_labels,
        font_size=font_size - 1,
        font_family="monospace",
        label_pos=0.5,
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
