# Edge-Weighted Balanced Connected Partitions

This repository contains the code accompanying the research paper:

**Morteza Davari, Phablo F. S. Moura, and Hande Yaman**  
*Edge-Weighted Balanced Connected Partitions: Hardness and Formulations*, 2025.  
[arXiv:2504.02421](https://arxiv.org/abs/2504.02421)

> 📌 **Note:** This version of the repository corresponds to the codebase used for the initial submission of the paper.

## Overview

This code provides implementations of various exact algorithms and mathematical formulations for the **Edge-Weighted Balanced Connected Partition** problem. The objective is to partition a graph into connected components that are balanced (e.g., in size or weight) while minimizing the total cost of cut edges.

## Directory Structure

- **`TIF/`**  
  Contains the full implementation of all methods described in the paper, including different formulations and algorithms.

- **`HeaviestBalancedTree/`**  
  Includes the code used for generating test instances.

## Compilation

The codebase was developed using **Visual Studio C++ 2025**. For Linux-based systems, you can compile the code using **CMake**. Ensure that the following libraries are installed and properly linked:

- [Boost](https://www.boost.org/)
- [LEMON](https://lemon.cs.elte.hu/trac/lemon)
- [IBM ILOG CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)
- [SFML](https://www.sfml-dev.org/) (Simple and Fast Multimedia Library) — if needed for visualization or I/O

> ⚠️ Make sure CPLEX is licensed and properly set up in your environment.

## Citation

If you use this code in your research, please cite the paper as follows:

```bibtex
@misc{davari2025edgeweightedbalancedconnectedpartitions,
  title        = {Edge-Weighted Balanced Connected Partitions: Hardness and Formulations},
  author       = {Morteza Davari and Phablo F. S. Moura and Hande Yaman},
  year         = {2025},
  eprint       = {2504.02421},
  archivePrefix = {arXiv},
  primaryClass = {cs.D
