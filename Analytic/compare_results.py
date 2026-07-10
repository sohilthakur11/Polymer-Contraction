"""
Comparison Analysis: Theta vs. Tree
Calculate the relative contraction factor ratio
"""

import numpy as np
from fractions import Fraction

print("="*90)
print("RELATIVE CONTRACTION FACTOR ANALYSIS")
print("Theta Graph (Group 5) vs. Branched Tree (Group 1)")
print("="*90)

theta_data = np.load('theta_results.npz')
tree_data = np.load('tree_results.npz')

g_theta = float(theta_data['g_theta'])
g_tree = float(tree_data['g_tree'])
trace_theta = float(theta_data['trace_theta'])
trace_tree = float(tree_data['trace_tree'])
v_theta = int(theta_data['v_theta'])
e_theta = int(theta_data['e_theta'])
loops_theta = int(theta_data['loops_theta'])
v_tree = int(tree_data['v_tree'])
e_tree = int(tree_data['e_tree'])
loops_tree = int(tree_data['loops_tree'])

print("\nLOADED RESULTS:")
print("-" * 90)
print(f"\nTheta Graph:")
print(f"  g(∞-theta) = {g_theta:.10f}")
print(f"  Trace = {trace_theta:.6f}")
print(f"  v={v_theta}, e={e_theta}, loops={loops_theta}")

print(f"\nBranched Tree:")
print(f"  g(∞-tree) = {g_tree:.10f}")
print(f"  Trace = {trace_tree:.6f}")
print(f"  v={v_tree}, e={e_tree}, loops={loops_tree}")

print("\n" + "="*90)
print("RELATIVE CONTRACTION FACTOR (Definition 21)")
print("="*90)

ratio = g_theta / g_tree
ratio_frac = Fraction(ratio).limit_denominator(1000)

print(f"\nRelative contraction factor:")
print(f"  g(theta)/g(tree) = {g_theta:.10f} / {g_tree:.10f}")
print(f"                   = {ratio:.10f}")
print(f"                   ≈ {ratio_frac}")
print("\n" + "="*90)
print("COMPARISON TABLE")
print("="*90)

print(f"""
┌─────────────────────┬──────────────────┬──────────────────┐
│ Parameter           │ Theta (Group 5)  │ Tree (Group 1)   │
|─────────────────────┼──────────────────┼──────────────────┤
│ Vertices            │ {v_theta:<16} │ {v_tree:<16} │
│ Edges               │ {e_theta:<16} │ {e_tree:<16} │
│ Loops               │ {loops_theta:<16} │ {loops_tree:<16} │
│ Tr(L̄†)              │ {trace_theta:<16.6f} │ {trace_tree:<16.6f} │
│ g(G∞)               │ {g_theta:<16.10f} │ {g_tree:<16.10f} │
└─────────────────────┴──────────────────┴──────────────────┘

FINAL RATIO: {ratio:.6f} ≈ {ratio_frac}""")

print("="*90)
print("Analysis complete!")
print("="*90)

print(f"""
\nCONCLUSION:
The theta polymer is approximately {100*ratio:.2f}% as compact as the tree 
reference polymer in the infinite subdivision limit.

This means the theta topology (with 3 loops) causes more compaction than
the acyclic tree structure, as predicted by topological polymer theory.
""")
