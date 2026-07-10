"""
Theta Graph Analysis
Polymer topology contraction factor calculation using Theorem 5
"""

import networkx as nx
import numpy as np
from numpy.linalg import eigvals, pinv
from fractions import Fraction

print("="*90)
print("GROUP 5: THETA GRAPH ANALYSIS")
print("Hexagonal ring with central vertex connected to alternating vertices")
print("="*90)

print("\nSTEP 1: Build the Theta Graph")
print("-" * 90)

G_theta = nx.Graph()

# Hexagon ring: vertices 0-5
edges_ring = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
G_theta.add_edges_from(edges_ring)

# Central vertex (6) connects to vertices 1, 3, 5
edges_branches = [(6,1), (6,3), (6,5)]
G_theta.add_edges_from(edges_branches)

# Graph properties
v_theta = G_theta.number_of_nodes()
e_theta = G_theta.number_of_edges()
loops_theta = e_theta - v_theta + 1

print(f"Vertices (v): {v_theta}")
print(f"Edges (e): {e_theta}")
print(f"Loops (cycle rank): {loops_theta}")
print(f"\nEdge list:")
print(f"  Ring: {edges_ring}")
print(f"  Branches: {edges_branches}")
print(f"\nDegree sequence: {dict(G_theta.degree())}")

print("\nSTEP 2: Compute Adjacency Matrix (A)")
print("-" * 90)

A = nx.to_numpy_array(G_theta)
print("Adjacency Matrix A:")
print(np.array(A, dtype=int))

print("\nSTEP 3: Compute Degree Matrix (D)")
print("-" * 90)

D = np.diag([G_theta.degree(i) for i in G_theta.nodes()])
print("Degree Matrix D:")
print(np.array(D, dtype=int))

print("\nSTEP 4: Compute Laplacian Matrix (L = D - A)")
print("-" * 90)

L = D - A
print("Laplacian Matrix L:")
print(np.array(L, dtype=int))

print("\nSTEP 5: Compute Normalized Laplacian Matrix")
print("-" * 90)

L_norm = nx.normalized_laplacian_matrix(G_theta).toarray()
print("Normalized Laplacian L̄:")
print(np.round(L_norm, 4))

print("\nDefinition: L̄_ij = 1 if i=j")
print("                  = -1/√(d_i × d_j) if i and j are adjacent")
print("                  = 0 otherwise")
print("\nSTEP 6: Compute Eigenvalues of Normalized Laplacian")
print("-" * 90)

eigs_theta = eigvals(L_norm)
eigs_sorted = np.sort(np.real(eigs_theta))

print("Eigenvalues (sorted):")
for i, eig in enumerate(eigs_sorted):
    print(f"  λ_{i} = {eig:.10f}")

print("\nSTEP 7: Compute Moore-Penrose Pseudoinverse and Trace")
print("-" * 90)

L_norm_pinv = pinv(L_norm)
trace_theta = np.trace(L_norm_pinv)

print(f"Trace of L̄† (direct): {trace_theta:.10f}")

nonzero_eigs = [e for e in eigs_sorted if abs(e) > 1e-10]
sum_inv_eigs = sum(1/e for e in nonzero_eigs)

print(f"\nVerification using eigenvalues:")
print(f"Number of nonzero eigenvalues: {len(nonzero_eigs)}")
print(f"Sum of 1/λ_i (nonzero): {sum_inv_eigs:.10f}")
print(f"Match? {np.isclose(trace_theta, sum_inv_eigs)}")

trace_frac = Fraction(trace_theta).limit_denominator(1000)
print(f"\nTrace as fraction: {trace_frac} = {float(trace_frac)}")

print("\nSTEP 8: Apply Theorem 5 - Contraction Factor Formula")
print("-" * 90)

print("\nTheorem 5 (Cantarella et al., 2022):")
print("g(G∞) = 3 × [Tr(L̄†) + (1/3)×Loops(G) - 1/6] / e²")

numerator = 3 * (trace_theta + loops_theta/3 - 1/6)
denominator = e_theta**2

print(f"\nCalculation:")
print(f"  Tr(L̄†) = {trace_theta:.6f}")
print(f"  Loops = {loops_theta}")
print(f"  e = {e_theta}")
print(f"\n  Numerator = 3 × ({trace_theta:.6f} + {loops_theta}/3 - 1/6)")
print(f"            = 3 × ({trace_theta:.6f} + {loops_theta/3:.6f} - 0.166667)")
print(f"            = 3 × {trace_theta + loops_theta/3 - 1/6:.6f}")
print(f"            = {numerator:.6f}")
print(f"\n  Denominator = e² = {e_theta}² = {denominator}")

g_theta = numerator / denominator

print(f"\n  g(∞-theta) = {numerator:.6f} / {denominator}")
print(f"             = {g_theta:.10f}")

g_theta_frac = Fraction(g_theta).limit_denominator(10000)
print(f"             ≈ {g_theta_frac}")

print("\n" + "="*90)
print("SUMMARY: THETA GRAPH")
print("="*90)
print(f"""
Topology: Hexagonal ring + central branch
Vertices: {v_theta}
Edges: {e_theta}
Loops: {loops_theta}
Trace(L̄†): {trace_theta:.6f}
g(∞-theta): {g_theta:.10f} ≈ {g_theta_frac}
""")
print("="*90)
print("Analysis complete! Results saved for comparison with tree graph.")
print("="*90)
np.savez('theta_results.npz', 
         g_theta=g_theta, 
         trace_theta=trace_theta,
         v_theta=v_theta,
         e_theta=e_theta,
         loops_theta=loops_theta)

print("\nResults saved to: theta_results.npz")
