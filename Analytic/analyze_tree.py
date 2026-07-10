"""
Branched Tree Analysis
Polymer topology contraction factor calculation using Theorem 5
"""

import networkx as nx
import numpy as np
from numpy.linalg import eigvals, pinv
from fractions import Fraction

print("="*90)
print("BRANCHED TREE ANALYSIS")
print("1 central node → 3 secondary nodes → 6 leaf nodes")
print("="*90)

print("\nSTEP 1: Build the Branched Tree Graph")
print("-" * 90)
G_tree = nx.Graph()
edges_primary = [(0, 1), (0, 2), (0, 3)]
G_tree.add_edges_from(edges_primary)

edges_secondary = [
    (1, 4), (1, 5),      
    (2, 6), (2, 7),      
    (3, 8), (3, 9)      
]
G_tree.add_edges_from(edges_secondary)
v_tree = G_tree.number_of_nodes()
e_tree = G_tree.number_of_edges()
loops_tree = e_tree - v_tree + 1
print(f"Vertices (v): {v_tree}")
print(f"Edges (e): {e_tree}")
print(f"Loops (cycle rank): {loops_tree} (acyclic tree)")
print(f"\nEdge list:")
print(f"  Primary: {edges_primary}")
print(f"  Secondary: {edges_secondary}")
print(f"\nDegree sequence: {dict(G_tree.degree())}")
print("\nSTEP 2: Compute Adjacency Matrix (A)")
print("-" * 90)

A = nx.to_numpy_array(G_tree)
print("Adjacency Matrix A:")
print(np.array(A, dtype=int))
print("\nSTEP 3: Compute Degree Matrix (D)")
print("-" * 90)

D = np.diag([G_tree.degree(i) for i in G_tree.nodes()])
print("Degree Matrix D:")
print(np.array(D, dtype=int))
print("\nSTEP 4: Compute Laplacian Matrix (L = D - A)")
print("-" * 90)

L = D - A
print("Laplacian Matrix L:")
print(np.array(L, dtype=int))
print("\nSTEP 5: Compute Normalized Laplacian Matrix")
print("-" * 90)

L_norm = nx.normalized_laplacian_matrix(G_tree).toarray()
print("Normalized Laplacian L̄:")
print(np.round(L_norm, 4))

print("\nDefinition: L̄_ij = 1 if i=j")
print("                  = -1/√(d_i × d_j) if i and j are adjacent")
print("                  = 0 otherwise")
print("\nSTEP 6: Compute Eigenvalues of Normalized Laplacian")
print("-" * 90)

eigs_tree = eigvals(L_norm)
eigs_sorted = np.sort(np.real(eigs_tree))

print("Eigenvalues (sorted):")
for i, eig in enumerate(eigs_sorted):
    print(f"  λ_{i} = {eig:.10f}")
print("\nSTEP 7: Compute Moore-Penrose Pseudoinverse and Trace")
print("-" * 90)

L_norm_pinv = pinv(L_norm)
trace_tree = np.trace(L_norm_pinv)

print(f"Trace of L̄† (direct): {trace_tree:.10f}")

nonzero_eigs = [e for e in eigs_sorted if abs(e) > 1e-10]
sum_inv_eigs = sum(1/e for e in nonzero_eigs)

print(f"\nVerification using eigenvalues:")
print(f"Number of nonzero eigenvalues: {len(nonzero_eigs)}")
print(f"Sum of 1/λ_i (nonzero): {sum_inv_eigs:.10f}")
print(f"Match? {np.isclose(trace_tree, sum_inv_eigs)}")

trace_frac = Fraction(trace_tree).limit_denominator(1000)
print(f"\nTrace as fraction: {trace_frac} = {float(trace_frac)}")
print("\nSTEP 8: Apply Theorem 5 - Contraction Factor Formula")
print("-" * 90)

print("\nTheorem 5 (Cantarella et al., 2022):")
print("g(G∞) = 3 × [Tr(L̄†) + (1/3)×Loops(G) - 1/6] / e²")

numerator = 3 * (trace_tree + loops_tree/3 - 1/6)
denominator = e_tree**2

print(f"\nCalculation:")
print(f"  Tr(L̄†) = {trace_tree:.6f}")
print(f"  Loops = {loops_tree}")
print(f"  e = {e_tree}")
print(f"\n  Numerator = 3 × ({trace_tree:.6f} + {loops_tree}/3 - 1/6)")
print(f"            = 3 × ({trace_tree:.6f} + {loops_tree/3:.6f} - 0.166667)")
print(f"            = 3 × {trace_tree + loops_tree/3 - 1/6:.6f}")
print(f"            = {numerator:.6f}")
print(f"\n  Denominator = e² = {e_tree}² = {denominator}")

g_tree = numerator / denominator

print(f"\n  g(∞-tree) = {numerator:.6f} / {denominator}")
print(f"            = {g_tree:.10f}")

g_tree_frac = Fraction(g_tree).limit_denominator(10000)
print(f"            ≈ {g_tree_frac}")
print("\n" + "="*90)
print("SUMMARY: BRANCHED TREE (Group 1)")
print("="*90)

print(f"""
Topology: 1 central → 3 secondary → 6 leaves
Vertices: {v_tree}
Edges: {e_tree}
Loops: {loops_tree}
Trace(L̄†): {trace_tree:.6f}
g(∞-tree): {g_tree:.10f} ≈ {g_tree_frac}
""")

print("="*90)
print("Analysis complete! Results saved for comparison.")
print("="*90)
np.savez('tree_results.npz', 
         g_tree=g_tree, 
         trace_tree=trace_tree,
         v_tree=v_tree,
         e_tree=e_tree,
         loops_tree=loops_tree)

print("\nResults saved to: tree_results.npz")
