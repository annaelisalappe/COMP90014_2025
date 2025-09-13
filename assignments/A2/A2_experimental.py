import sys
import argparse
import resource
import time
import os
import networkx as nx
import matplotlib.pyplot as plt
from Bio import Seq, Align
from hashlib import blake2b
import Levenshtein
from itertools import combinations


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

    contigs = {}
    for idx in range(5):
        G = build_de_bruijn_graph(kmer_counts, k_size, idx)
        # draw_debruijn_graph(
        #     G, title=f"De Bruijn Graph for organism {chr(ord('A') + idx)}"
        # )
        H = collapse_linear_paths(G)
        # draw_debruijn_graph(
        #     H,
        #     title=f"Collapsed De Bruijn Graph for organism {chr(ord('A') + idx)}",
        #     node_size=300,
        #     font_size=6,
        # )
        S = remove_tips(H, min_length=1.5 * k_size, min_coverage=1, k=k_size)
        draw_debruijn_graph(
            S,
            title=f"Tip-removed De Bruijn Graph for organism {chr(ord('A') + idx)}",
            node_size=300,
            font_size=6,
        )
        contigs[chr(ord("A") + idx)] = extract_contigs(S)
        print(f"Organism {chr(ord('A') + idx)}: {contigs[chr(ord('A') + idx)]}")

    resistant_contigs = list(set(contigs["A"] + contigs["B"] + contigs["C"]))
    sensitive_contigs = list(set(contigs["D"] + contigs["E"]))

    print(f"Resistant contigs: {resistant_contigs}")
    print(f"Sensitive contigs: {sensitive_contigs}")

    unique_contigs = find_unique_sequence(
        resistant_contigs, sensitive_contigs, min_plasmid_len
    )
    print(unique_contigs)
    print("Rotated unique sequences:")
    unique_contigs = [rotate_sequence(uc) for uc in unique_contigs]
    unique_contigs = list(set(unique_contigs))
    unique_contigs.sort(key=len, reverse=True)

    print(unique_contigs)
    with open(outfile, "w") as out:
        for uc in unique_contigs:
            if len(uc) >= min_plasmid_len:
                out.write(uc + "\n")


def get_kmers(seq: str, k: int) -> set[str]:
    return [seq[i : i + k] for i in range(len(seq) - k + 1)]


def find_unique_sequence(
    resistant_contigs: list[str], sensitive_contigs: list[str], k: int
) -> set[str]:
    sensitive_kmers = []
    for seq in sensitive_contigs:
        sensitive_kmers += get_kmers(seq, k)
    sensitive_kmers = set(sensitive_kmers)
    print(sensitive_kmers)

    resistant_kmers = []
    for seq in resistant_contigs:
        resistant_kmers += get_kmers(seq, k)
    resistant_kmers = set(resistant_kmers)
    print(resistant_kmers)

    return resistant_kmers - sensitive_kmers


def hash_seq(seq: str) -> str:
    """Return a stable hash for the sequence."""
    h = blake2b(seq.encode(), digest_size=16)  # 128-bit digest
    return h.hexdigest()


def find_unique_contigs(
    resistant_contigs: list[str], sensitive_contigs: list[str]
) -> list[str]:
    """
    Return contigs that are present in resistant_contigs
    but not in sensitive_contigs (exact string match).
    """
    # Build a set of hashed sensitive contigs for fast lookup
    sensitive_hashes: set[str] = {hash_seq(contig) for contig in sensitive_contigs}

    # Keep only resistant contigs not present in the sensitive set
    return [
        contig
        for contig in resistant_contigs
        if hash_seq(contig) not in sensitive_hashes
    ]

    # get contigs from graph
    # make pools of candidate and rest contigs
    # do sequence alignment to find the uncovered regions
    # maybe reuse some code from a1?

    # candidate_kmers, rest_kmers = select_candidate_rest_kmers(kmer_counts, min_diff=0)
    # print(f"Candidate k-mers: {len(candidate_kmers)}")
    # print(candidate_kmers)
    # print(f"Rest k-mers: {len(rest_kmers)}")
    # print(rest_kmers)

    # candidate_G = build_de_bruijn_graph(candidate_kmers, k_size)
    # candidate_contigs = extract_contigs(candidate_G)

    # rest_G = build_de_bruijn_graph(rest_kmers, k_size)
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


# def extract_plasmid_regions(candidate_contigs, rest_contigs, min_plasmid_len):
#     """
#     Align each candidate contig to each rest contig using global alignment.
#     Return plasmid-specific subsequences (uncovered regions) as a dictionary.
#     """
#     all_possible_plasmids = {}
#     aligner = Align.PairwiseAligner(
#         match_score=1.0,
#         mode="fogsaa",
#         open_gap_score=0,
#         extend_gap_score=0,
#         mismatch_score=-1.0,
#     )

#     for cand_str in candidate_contigs:
#         cand = Seq.Seq(cand_str)
#         coverage = [False] * len(cand)  # mark positions covered by any rest contig

#         # Align candidate to all rest contigs
#         for rest_str in rest_contigs:
#             rest = Seq.Seq(rest_str)
#             alignments = aligner.align(cand, rest)
#             if not alignments:
#                 print("No alignments found.")
#                 continue

#             best = alignments[0]
#             print(best)
#             # convert alignment to strings
#             cand_algn = str(best[0])
#             rest_algn = str(best[1])

#             possible_plasmid = ""
#             for i in range(len(cand_algn)):
#                 if cand_algn[i] != "-":
#             # mark positions covered by this alignment
#             idx_cand = 0
#             for i in range(len(cand_algn)):
#                 if cand_algn[i] != "-":  # skip gaps in candidate
#                     if rest_algn[i] != "-":
#                         coverage[idx_cand] = True
#                     idx_cand += 1

#         # Extract uncovered regions
#         possible_plasmids = []
#         plasmid_seq = ""
#         for i, is_cov in enumerate(coverage):
#             if not is_cov:
#                 plasmid_seq += cand[i]
#             else:
#                 if len(plasmid_seq) >= 3:
#                     possible_plasmids.append(plasmid_seq)
#                 plasmid_seq = ""
#         if len(plasmid_seq) >= min_plasmid_len:
#             possible_plasmids.append(plasmid_seq)

#         all_possible_plasmids[cand_str] = possible_plasmids

#     return all_possible_plasmids


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


# def extract_contigs(G: nx.DiGraph) -> list[str]:
#     """
#     Extracts all contigs from a de Bruijn graph.
#     """
#     start_nodes = []
#     for node in G.nodes():
#         if G.in_degree(node) == 0 or G.in_degree(node) > 1:
#             start_nodes.append(node)
#         elif G.out_degree(node) > 1:
#             child_nodes = list(G.successors(node))
#             start_nodes += child_nodes

#     contigs = []
#     for start in start_nodes:
#         contig = start
#         current_node = start
#         print(f"Start node: {start}")
#         while True:
#             out_edges = list(G.out_edges(current_node, data=True))
#             if len(out_edges) == 0 or any(
#                 succ in start_nodes for succ in G.successors(current_node)
#             ):
#                 print(current_node, out_edges)
#                 break
#             # assuming only one outgoing edge
#             _, next_node, edge_data = out_edges[0]
#             contig += edge_data["kmer"][-1]
#             current_node = next_node
#         contigs.append(contig)
#     return contigs


def extract_contigs(G: nx.DiGraph, min_count=1) -> list[str]:
    """
    Extracts all contigs from a de Bruijn graph.
    """
    start_nodes = []
    for node in G.nodes():
        if G.in_degree(node) == 0 or G.in_degree(node) > 1:
            start_nodes.append(node)
        elif G.out_degree(node) > 1:
            child_nodes = list(G.successors(node))
            start_nodes += child_nodes

    contigs = []
    for start in start_nodes:
        contig = start
        current_node = start
        while True:
            out_edges = list(G.out_edges(current_node, data=True))
            if len(out_edges) == 0 or any(
                succ in start_nodes for succ in G.successors(current_node)
            ):
                break
            # assuming only one outgoing edge
            _, next_node, edge_data = out_edges[0]
            contig += edge_data["label"][-1]
            current_node = next_node
        contigs.append(contig)
    return contigs


# def tour_bus(G: nx.DiGraph) -> nx.DiGraph:
#     visited = set()
#     node = G.anynode

#     for succ in G.successors(node):
#         while (node, succ) not in visited:
#             visited.add((node, succ))
#             continue


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


def build_de_bruijn_graph(kmer_counts: dict, k: int, id: int) -> dict:
    G = nx.DiGraph()
    for km, counts in kmer_counts.items():
        if counts[id] > 0:  # k-mer appears in this organism
            prefix = km[:-1]
            suffix = km[1:]
            # If edge already exists, increment its weight
            if G.has_edge(prefix, suffix):
                G[prefix][suffix]["cov"] += counts[id]
            else:
                G.add_edge(prefix, suffix, kmer=km, cov=counts[id])
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
        prefix = start[:-1]
        suffix = current_node[1:]
        # if length == 0:
        #     in_edges = list(G.in_edges(current_node, data=True))

        #     if len(in_edges) == 0:

        #     if len(list(G.in_edges(current_node, data=True))) > 0:
        #         print("Warning: length 0 but has mutliple in-edges")
        #     avg_count = list(G.in_edges(current_node, data=True))[0][2]["cov"]
        # else:
        avg_count = counts / length if length > 0 else counts
        H.add_edge(prefix, suffix, seq=contig, cov=avg_count)

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
    print(
        f"Before removing: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges. Number of nodes with in-degree 0 and out-degree 0: {sum(1 for n in G.nodes() if G.in_degree(n) == 0 and G.out_degree(n) == 0)}"
    )
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
                        G.remove_nodes_from(path_nodes[:-1])  # keep the junction node
                        removed = True

    print(
        f"After removing: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges. Number of nodes with in-degree 0 and out-degree 0: {sum(1 for n in G.nodes() if G.in_degree(n) == 0 and G.out_degree(n) == 0)}"
    )

    for nodes in list(G.nodes):
        if G.in_degree(nodes) == 0 and G.out_degree(nodes) == 0:
            G.remove_node(nodes)

    return G


def _trace_path(G, start_node, end_node_of_other_path, max_len):
    """Helper function to trace a single path until a merge node or max_len."""
    path_edges = []

    # The first edge comes from the start_node
    try:
        u, v, data = list(G.out_edges(start_node, data=True))[0]
        path_edges.append((u, v, data))
    except IndexError:
        return None  # Dead end

    curr = v
    length = len(data["seq"])

    # Trace until we hit the other path's end, a junction, or max_len
    while (
        curr != end_node_of_other_path
        and G.out_degree(curr) == 1
        and G.in_degree(curr) == 1
        and length < max_len
    ):
        u, v, data = list(G.out_edges(curr, data=True))[0]
        path_edges.append((u, v, data))
        # Add k-1 overlap
        length += len(data["seq"]) - (len(u))
        curr = v

    # If we found a valid end point, return path info
    if G.out_degree(curr) != 1 or length >= max_len:
        return None  # Path ended at a different junction or got too long

    # Successfully found a converging path
    full_seq = path_edges[0][2]["seq"]
    total_cov = path_edges[0][2]["cov"] * len(full_seq)

    for i in range(1, len(path_edges)):
        edge_seq = path_edges[i][2]["seq"]
        k = len(path_edges[i][0])
        full_seq += edge_seq[k:]
        total_cov += path_edges[i][2]["cov"] * (len(edge_seq) - k)

    avg_cov = total_cov / len(full_seq)

    return {
        "end_node": curr,
        "edges": [e[:2] for e in path_edges],
        "seq": full_seq,
        "cov": avg_cov,
    }


def remove_bubbles(
    G: nx.DiGraph, max_bubble_length: int = 1000, similarity_threshold: float = 0.95
) -> nx.DiGraph:
    """
    Finds and removes bubbles from a compacted de Bruijn graph.

    A bubble is a pair of similar paths between two junction nodes. This function
    identifies such structures, compares the paths based on sequence similarity
    and coverage, and removes the lower-coverage path.

    Args:
        G: A compacted de Bruijn graph from collapse_linear_paths.
           Edges must have 'seq' and 'cov' attributes.
        max_bubble_length: The maximum sequence length for a path to be
                           considered part of a bubble.
        similarity_threshold: The sequence similarity required to classify
                              two paths as a bubble (e.g., 0.95 means 95% similar).

    Returns:
        A new NetworkX DiGraph with bubbles removed.
    """
    H = G.copy()
    bubbles_found_and_removed = True

    while bubbles_found_and_removed:
        bubbles_found_and_removed = False
        # Find nodes where bubbles might start (a "fork")
        fork_nodes = [n for n in H.nodes() if H.out_degree(n) > 1]

        for u in fork_nodes:
            # Get all pairs of outgoing paths
            successors = list(H.successors(u))
            for v1, v2 in combinations(successors, 2):
                # Get the sequences and coverages of the two edges leaving the fork
                edge1_data = H.get_edge_data(u, v1)
                edge2_data = H.get_edge_data(u, v2)

                path1_seq = edge1_data["seq"]
                path2_seq = edge2_data["seq"]

                # Check for simple bubbles (single edge pairs)
                # Do they merge at the same node?
                if v1 == v2:
                    continue  # This case is complex, skip for now

                # Do the paths from v1 and v2 merge?
                # We need to find if there is a common successor node
                try:
                    # In a simple bubble, v1 and v2 are the only predecessors of the merge node
                    merge_node = next(
                        w
                        for w in H.successors(v1)
                        if H.has_edge(v2, w) and H.in_degree(w) == 2
                    )
                except StopIteration:
                    continue  # No simple merge node found

                # We found a bubble from u -> v1/v2 -> merge_node
                path1_seq += H.get_edge_data(v1, merge_node)["seq"][len(v1) :]
                path2_seq += H.get_edge_data(v2, merge_node)["seq"][len(v2) :]

                # Check length constraint
                if (
                    len(path1_seq) > max_bubble_length
                    or len(path2_seq) > max_bubble_length
                ):
                    continue

                # Check similarity
                similarity = Levenshtein.ratio(path1_seq, path2_seq)
                if similarity >= similarity_threshold:
                    # Calculate average coverage for each path
                    cov1 = (
                        edge1_data["cov"] + H.get_edge_data(v1, merge_node)["cov"]
                    ) / 2
                    cov2 = (
                        edge2_data["cov"] + H.get_edge_data(v2, merge_node)["cov"]
                    ) / 2

                    # Identify the path to remove
                    if cov1 >= cov2:
                        path_to_remove = [(u, v2), (v2, merge_node)]
                        node_to_remove = v2
                    else:
                        path_to_remove = [(u, v1), (v1, merge_node)]
                        node_to_remove = v1

                    # Remove the lower-coverage path
                    H.remove_edges_from(path_to_remove)
                    # Remove the intermediate node if it's now isolated
                    if (
                        H.in_degree(node_to_remove) == 0
                        and H.out_degree(node_to_remove) == 0
                    ):
                        H.remove_node(node_to_remove)

                    print(
                        f"Removed bubble between {u} and {merge_node}. Kept path with coverage {max(cov1, cov2):.2f}."
                    )
                    bubbles_found_and_removed = True
                    break  # Restart scan since graph has changed
            if bubbles_found_and_removed:
                break

    return H


def assemble_contigs(G: nx.DiGraph) -> list[str]:
    """
    Extracts all contigs from a de Bruijn graph.

    A contig is a maximal non-branching path.
    """
    visited = set()
    H = nx.DiGraph()

    def is_non_branching(node):
        return G.in_degree(node) == 1 and G.out_degree(node) == 1

    # Pass 1: extract paths from sources or branch points
    for node in G.nodes:
        if G.out_degree(node) > 0 and not is_non_branching(node):
            start = node
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
                prefix = start[:-1]
                suffix = succ[1:]
                # avg_count = counts // (length + 1) if length > 0 else counts
                H.add_edge(prefix, suffix, kmer=contig)
    return H

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
