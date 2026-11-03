# CS598APE HW3

This repository contains code for homework 3 of CS598APE. It is relatively
simple in comparison to HW1 and HW2 to ensure you have enough time to work on
the course project. 

`sandpile.c` contains the implementation of finding the identity element of a 2D
abelian sandpile.

![Identity sandpile on a 1001x1001 pacman grid](identity-ref.png)

*Identity sandpile on a 1001x1001 pacman grid*


To compile the program, run:

```bash
make all
```

To clean existing build artifacts, run:

```bash
make clean
```

Usage (after compiling):

```bash
./sandpile.exe <N> 
```

where `N` is the grid size

Consider speeding up the generation of this 1001x1001 identity sandpile (which
takes approximately ~3m 40s on my local laptop under the default setup)

```bash
./sandpile.exe 1001
```

For an optional, ungraded challenge, see how large a sandpile you can compute
on the course VMs and report the grid size and runtime. There is no prescribed
upper limit on `N` for this challenge.

## Reproducing the Report Results

The report's numbers were measured on a course-provided virtual machine with 4 vCPUs (Intel Xeon Silver 4216) and 15 GB of RAM. This program assumes `gcc` and `make` are installed:

### Optimizations

Each commit in the git history after the baseline commit represents one optimization from the report.

| # | Report optimization | Commit |
|---|---------------------|--------|
| 0 | Baseline | `15fe837` |
| 1 | Precompute PacMan slices | `a65b7a1` |
| 2 | Topple in batch | `a0c82af` |
| 3 | Improved grid data structure | `e7ea5c9` |
| 4 | Pull grains from neighbors | `bc7e6fb` |
| 5 | Precompute a toppling lower bound | `5636973` |
| 6 | Store tiles as `int8_t`s | `a261a1a` |
| 7 | Iterate the grid in a checkerboard pattern | `46f3593` |
| 8 | Parallelize loops | `89f5a4c` |

### Build

```bash
git checkout <commit hash>
make clean && make
```

### Benchmark

Run the program three times and take the geometric mean of the wall clock times to reproduce the results from the speedup figure in the report. Run from the repository root.

```bash
time ./sandpile.exe 1001
```
