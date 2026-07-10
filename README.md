# Polymer Contraction Factor Analysis

**Analytical Derivation and Molecular Dynamics Validation of Topological Polymer Compaction**

A spectral graph theory approach to computing contraction factors for topological polymers using Python, with validation against molecular dynamics simulations.

---

## Overview

This repository contains the complete analytical and computational framework for deriving **contraction factors** (g-factors) for topological polymers, demonstrating how molecular architecture controls material properties.

### Key Results

- **Theta Polymer:** g(∞) ≈ 0.264 (analytical)
- **Branched Tree Polymer:** g(∞) ≈ 0.605 (analytical)
- **Relative Contraction Factor:** g(θ)/g(tree) = **107/245 ≈ 0.4367**
- **MD Validation:** R² = 0.966 (theory correlates with simulations)

The theta polymer is approximately **44% as compact** as the tree reference polymer. Theoretical predictions show excellent correlation (R² = 0.966) with molecular dynamics simulations, validating that topological structure fundamentally determines polymer behavior.

---

## Theoretical Foundation

### Theorem 5 (Cantarella et al., 2022)

For any connected graph G with v vertices and e edges, the contraction factor in the infinite edge-subdivision limit is:

$$g(G_\infty) = \frac{3\left[\operatorname{Tr}(\bar{\mathbf{L}}^\dagger) + \frac{1}{3}\text{Loops}(G) - \frac{1}{6}\right]}{e^2}$$

where:
- $$\(\bar{\mathbf{L}}\) = normalized Laplacian matrix of the base graph
- $$\(\bar{\mathbf{L}}^\dagger\) = Moore-Penrose pseudoinverse
- $$\(\text{Loops}(G) = e - v + 1\) = cycle rank
- $$\(\operatorname{Tr}(\bar{\mathbf{L}}^\dagger)\) = sum of reciprocal nonzero eigenvalues

### Physical Interpretation

The **contraction factor** quantifies how much more compact a polymer is compared to a linear reference chain of equal contour length. Topological features (loops, branches) create constraints that directly affect polymer size (radius of gyration) and conformational behavior.

---

## Project Structure

```
Polymer-Contraction/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
│
├── ANALYTICAL/
│   ├── analyze_theta.py              # Theta graph contraction factor analysis
│   ├── analyze_tree.py               # Branched tree contraction factor analysis
│   ├── compare_results.py            # Relative contraction factor & validation
│   └── results/
│       ├── theta_results.npz         # Saved theta calculations
│       ├── tree_results.npz          # Saved tree calculations
│       └── comparison_output.txt     # Final ratio & metrics
│
├── DOCUMENTATION/
│   └── polymer_report.pdf               # Compiled report
│
└── lammps/
```
---

## Files Description

### Analytical Scripts

#### `analyze_theta.py`
Computes contraction factor for the **theta polymer** (Group 5).
- **Input:** Theta graph topology (7 vertices, 9 edges, 3 loops)
- **Process:** 
  - Constructs adjacency, degree, and Laplacian matrices
  - Computes normalized Laplacian
  - Calculates eigenvalues and pseudoinverse trace
  - Applies Theorem 5
- **Output:** g(∞-theta) ≈ 0.264, saved to `theta_results.npz`

**Key metrics:**
- Vertices: 7 | Edges: 9 | Loops: 3
- Trace(L†): 6.3
- Contraction factor: 0.264 ≈ 107/405

#### `analyze_tree.py`
Computes contraction factor for the **branched tree** (Group 1).
- **Input:** Branched tree topology (10 vertices, 9 edges, 0 loops)
  - Central node → 3 secondary nodes → 6 leaves
- **Process:** Same analytical workflow as theta
- **Output:** g(∞-tree) ≈ 0.605, saved to `tree_results.npz`

**Key metrics:**
- Vertices: 10 | Edges: 9 | Loops: 0
- Trace(L†): 16.5
- Contraction factor: 0.605 ≈ 49/81

#### `compare_results.py`
Compares analytical results from both polymers.
- **Computation:** Relative contraction factor (theta/tree)
- **Validation:** Comparison with published theoretical value
- **Output:** Final ratio (107/245 ≈ 0.4367), comparison table

---

### Molecular Dynamics

#### `lammps/` folder
Contains all molecular dynamics simulations and analysis:
- **Input scripts:** LAMMPS simulation parameters and execution files
- **Initial configurations:** Polymer starting geometries
- **Simulations:** Trajectory outputs and simulation results
- **Analysis scripts:** Post-processing utilities to extract polymer properties from MD trajectories

---

### Documentation

#### `polymer_report.pdf`
Complete analytical derivation including:
- Mathematical preliminaries (graph theory, spectral analysis)
- Step-by-step derivations for both polymers
- Matrix calculations and eigenvalue analysis
- Theorem 5 application
- Physical interpretation and discussion
- Full references


---

## Installation & Usage

### Requirements

- Python 3.7+
- NumPy
- NetworkX
- SciPy

### Setup

```bash
git clone https://github.com/ShlokP06/Polymer-Contraction.git
cd Polymer-Contraction

pip install -r requirements.txt
```

### Running Analytical Calculations

```bash
# Analyze theta polymer
python ANALYTICAL/analyze_theta.py

# Analyze tree polymer
python ANALYTICAL/analyze_tree.py

# Compare results and compute ratio
python ANALYTICAL/compare_results.py
```

---

## Results Summary

### Analytical Results

| Polymer | Vertices | Edges | Loops | Trace(L†) | g(∞) | Exact Value |
|---------|----------|-------|-------|-----------|------|-------------|
| Theta   | 7        | 9     | 3     | 6.3       | 0.264 | 107/405     |
| Tree    | 10       | 9     | 0     | 16.5      | 0.605 | 49/81       |

### Relative Contraction Factor

$$\boxed{\frac{g(\theta)}{g(\text{tree})} = \frac{107}{245} \approx 0.4367}$$

**Interpretation:** The theta polymer is approximately **44% as compact** as the tree reference polymer, demonstrating that topological structure (cycles vs. branching) fundamentally determines polymer behavior.

---

## Polymer Topologies

### Theta Graph
- **Structure:** Hexagonal ring with central vertex connected to alternating vertices
- **Topological Feature:** 3 independent loops (cycles)
- **Compactness:** g = 0.264 (more extended)

### Branched Tree
- **Structure:** Central hub → 3 secondary branches → 6 leaves
- **Topological Feature:** Acyclic (no loops)
- **Compactness:** g = 0.605 (more compact)

---

## Applications

This framework enables:
- **Rational Polymer Design:** Predict properties before synthesis
- **Drug Delivery Optimization:** Design carriers with target compaction
- **Materials Engineering:** Structure-property relationships for elastomers, coatings
- **Computational Screening:** Compare topologies efficiently using graph theory

---

## References

**Cantarella, J., Deguchi, T., Shonkwiler, C., & Uehara, E. (2022).** 
Radius of gyration, contraction factors, and subdivisions of topological polymers. 
*Journal of Physics A: Mathematical and Theoretical*, 55(47), 475202.
https://doi.org/10.1088/1751-8121/aca300

---

## Project Credits

**Course:** CHE 209 - Soft Matter and Polymers  
**Institution:** Indian Institute of Technology Indore  
**Year:** 2025

---

## Citation

```bibtex
@github{PolymContractionGitHub2025,
  author = {Shlok Parikh, Uday Ranode, Utkarsh Sharma, Sohil Dangi, Shrawani Palange, Shruti Turare},
  title = {Polymer Contraction Factor Analysis: Spectral Graph Theory Approach},
  year = {2025},
  url = {https://github.com/ShlokP06/Polymer-Contraction}
}
```

---

## License

Educational project for CHE 209 coursework at IIT Indore.

---

**Last Updated:** November 2025
