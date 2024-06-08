# Stable randomized least-squares

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=eepperly/Randomized-Least-Squares-Solvers)

This repository contains code for the paper [_Fast randomized least-squares solvers can be just as accurate and stable as classical direct solvers_](https://arxiv.org/abs/2406.03468) by [Ethan N. Epperly](https://www.ethanepperly.com), [Maike Meier](https://www.maths.ox.ac.uk/people/maike.meier), and [Yuji Nakatsukasa](https://people.maths.ox.ac.uk/nakatsukasa/). It also contains code for the earlier paper [_Fast and forward stable randomized algorithms for linear least-squares problems_](https://arxiv.org/abs/2311.04362) by [Ethan N. Epperly](https://www.ethanepperly.com); see also [that paper's repo](https://github.com/eepperly/Iterative-Sketching-Is-Stable).

## Background: Fast, randomized least-squares methods

The overdetermined linear least-squares problem $$x = \text{argmin}_{x \in \mathbb{R}^n} \lVert b - Ax\rVert_2$$ is one of the most fundamental problem in machine learning, statistics, and computational mathematics.
To solve least-squares problems accurately, the standard approach is to use [QR factorization](https://en.wikipedia.org/wiki/QR_decomposition); QR-based least-squares solvers require $O(mn^2)$ operations for an $m\times n$ matrix $A$.

In 2008, Rokhlin and Tygert [introduced](https://www.pnas.org/doi/abs/10.1073/pnas.0804869105) the _sketch-and-precondition method_, a faster least-squares solver that uses _randomized dimensionality reduction_ to _precondition_ an iterative least-squares algorithm.
Sketch-and-precondition solves a least-squares problem in roughly $O(mn + n^3)$ operations, a significant speedup over QR for highly overdetermined least-squares problems $m\gg n\gg 1$.
Another class of randomized least-squares solvers, known as _iterative sketching methods_ (introduced by [Pilanci and Wainwright, 2016](https://www.jmlr.org/papers/v17/14-460.html)), achieve the same $O(mn + n^3)$ operation count as sketch-and-precondition and can have benefits in some situations (e.g., for parallel computing).

Until recently, the _numerical stability_ of randomized least-squares solvers was poorly understood.
How do these algorithms perform when the matrix $A$ is stored on a computer using double-, single-, or even half-precision [floating point numbers](https://en.wikipedia.org/wiki/Floating-point_arithmetic)?
The gold standard stability property for numerical algorithms is [_backward stability_](https://nhigham.com/2020/08/04/what-is-numerical-stability/), a property possessed by the classical QR method for least-squares.
Unfortunately, both sketch-and-precondition and iterative sketching [are not backward stable](https://arxiv.org/abs/2302.07202), though good implementations satisfy a weaker stability property known as forward stability.
Forward stability is enough for many applications, but—to truly use randomized least-squares solvers as a drop-in replacement for the QR method—it is desirable to have a fast randomized least-squares solver that is backward stable.

The paper [_Fast randomized least-squares solvers can be just as accurate and stable as classical direct solvers_](https://arxiv.org/abs/2406.03468) introduces two fast, backward stable randomized least-squares solvers:

- FOSSILS, a backward stable version of iterative sketching
- Sketch-and-precondition with iterative refinement (SPIR), a backward stable version of sketch-and-precondition

These methods are modified versions of iterative sketching and sketch-and-precondition that perform one step of iterative refinement. Fortunately, this one refinement step is enough to obtain backward stable algorithms.

## Code to replicate experiments from [_Fast randomized least-squares solvers can be just as accurate and stable as classical direct solvers_](https://arxiv.org/abs/2406.03468)

Code to reproduce the numerical experiments from this paper are found in `paper2_experiments/`. 

- Figure 1 (`paper2_experiments/test_comparison.m`): Computes forward and backward errors for different randomized least-squares solvers during the course of iteration.
- Table 1: (`paper2_experiments/test_residual_orthogonality.m`): Tests the value of $\lVert A^\top (b-Ax)\rVert` for different forward and backward stable least-squares solvers.
- Figure 2 (`paper2_experiments/test_comparison_2.m`): Computes the final forward and backward errors for different randomized least-squares solvers for problems of increasing difficulty.
- Figure 4 (`paper2_experiments/test_kernel.m` and `paper2_experiments/test_prony.m`): Compares the runtime of FOSSILS, iterative sketching with momentum, and MATLAB's `mldivide` on a kernel regression problem and an application of Prony's method to quantum eigenvalue estimation.
- Table 2: (`paper2_experiments/test_sparse.m`): Compares FOSSILS to `mldivide` for sparse least-squares problems from the [SuiteSparse matrix collection](https://sparse.tamu.edu).
- Figure 6: (`paper2_experiments/test_spir.m`): Numerical evaluation of the sketch-and-precondition with iterative refinement (SPIR) method.

## Code to replicate experiments from [_Fast and forward stable randomized algorithms for linear least-squares problems_](https://arxiv.org/abs/2311.04362)

Code to reproduce the numerical experiments from this paper are found in `paper1_experiments/`. 

- Figure 1 (`paper1_experiments/test_iterative_sketching.m`): Computes the forward, residual, and backward errors for iterative sketching for different condition numbers and residual norms.
- Figure 2 (`paper1_experiments/test_variants.m`): Compares the performance iterative sketching (including versions with damping and momentum) to sketch-and-precondition (with both the zero and sketch-and-solve initializations).
- Figure 3 (`paper1_experiments/test_bad_iterative_sketching.m`): Compare the stable implementation of iterative sketching (in `code/iterative_sketching.m`) to three "bad" implementations.
- Figure 4 (`paper1_experiments/test_timing.m`): Compares the runtime of iterative sketching (including versions with damping and momentum) to MATLAB's `mldivide` on a simplified kernel regression task.
- Figure 5 (`paper1_experiments/test_sparse.m`): Compares the runtime of iterative sketching to MATLAB's `mldivide` for solving random sparse least-squares problems with three nonzeros per row.

## Randomized least-squares solvers

This repository contains code for [iterative sketching](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS4) and [sketch-and-precondition](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS5).
Code for these methods is found in the `code/` directory.

### Sparse sign embeddings

The core ingredient for randomized least-squares solvers is a [_fast random embedding_](https://ar5iv.labs.arxiv.org/html/2002.01387#S9).
Based on the empirical comparison in [this paper (see Fig. 2)](https://arxiv.org/abs/2104.05877) and our [own testing](https://www.ethanepperly.com/index.php/2023/11/27/which-sketch-should-i-use/), our preferred embedding is the [_sparse sign embedding_](https://ar5iv.labs.arxiv.org/html/2002.01387#S9.SS2).

To generate a sparse sign embedding in our code, first build the [mex file](https://www.mathworks.com/help/matlab/ref/mex.html) using the following command:

```
mex sparsesign.c
```

Then, to generate a `d` by `m` sparse sign embedding with `zeta` nonzeros per column, run

```
S = sparsesign(d, m, zeta);
```

### FOSSILS

_FOSSILS_ (**F**ast **O**ptimal **S**table **S**ketchy **I**terative **L**east **S**quares) is a backward stable randomized least-squares solver with a fast roughly $O(mn + n^3)$ runtime.
The first step of the algorithm is to collect a sketch $SA$ of the matrix $A$ and compute its [QR factorization](https://en.m.wikipedia.org/wiki/QR_decomposition) $SA = QR$.
The core engine for FOSSILS is the following _outer solver_: On input $A,b$,

1. Compute $c = R^{-\top} (A^\top b)$.
2. Solve the equation $(R^{-\top} A^\top A R^{-1}) y = c$ using the Polyak heavy ball method.
3. Output $x = R^{-1}y$.

On it's own, the FOSSILS outer solver is _not_ backward stable, but it can be upgraded to backward stability by using _iterative refinement_.

1. Set $x_0 = R^{-1} (Q^\top (Sb))$. This is the "sketch and solve" initialization.
2. Update $x_1 = x_0 + \mathrm{FossilsOuterSolver}(b-Ax_0)$.
3. Output $x_2 = x_1 + \mathrm{FossilsOuterSolver}(b-Ax_1)$.

This simple two-step refinement procedure leads to a backward stable algorithm.

To run FOSSILS using our code, the command is

```
[x, stats] = fossils(A, b, [d, iterations, summary, verbose, reproducible])
```

All but the first two inputs are optional.
The optional inputs are as follows:

- `d`: sketching dimension. (_Default value_: `12*size(A,2)`)
- `iterations`: number of iterations for each refinement step (i.e., `[50,50]` to set fifty steps for each refinement step). Set to `'adaptive'` to automatically set the number of refinement steps. (_Default value_: `'adaptive'`)
- `summary`: a function of `x` to be recorded at every iteration. The results will be outputted in the optional output argument `stats`. (_Default value_: None)
- `verbose`: if true, output at every iteration. (_Default value_: `false`)
- `reproducible`: if true, use a slow, but reproducible implementation of sparse sign embeddings. (_Default value_: `false`)

Inputting a value of `[]` for an optional argument results in the default value.

### Sketch-and-precondition with iterative refinement

Sketch-and-precondition with iterative refinement (SPIR) is a fast, backward stable randomized least-squares solver.
It is similar to FOSSILS, except that it uses a Krylov method in place of the Polyak heavy ball method.

To run SPIR using our code, the command is

```
[x, stats] = spir(A, b, [d, q, summary, verbose, opts, reproducible])
```

All but the first two inputs are optional.
The optional inputs are as follows:

- `d`: sketching dimension. (_Default value_: `2*size(A,2)`)
- `q`: number of iterations for each refinement step. (_Default value_: `[50,50]`)
- `summary`: a function of `x` to be recorded at every iteration. The results will be outputted in the optional output argument `stats`. (_Default value_: None)
- `verbose`: if true, output at every iteration. (_Default value_: `false`)
- `opts`: specifies the initial iterate $x_0$ and iterative method (LSQR, symmetric conjugate gradient, or preconditioned conjugate gradient). If `'cold'` is a substring of `opts`, then the initial iterate is chosen to be $x_0 = 0$. Otherwise, we use a warm start and choose $x_0$ to be the [sketch-and-solve solution](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS3). If `cgne` is a substring of `opts`, then we solve $Ax = b$ using CGNE, and if `sym` is a substring of `opts`, then we perform conjugate gradient on the symmetrically preconditioned system $(R^{-\top} A^\top A R^{-1}) y = c$; otherwise, we use LSQR. We recommend using the LSQR or symmetric conjugate gradient implementations in practice. (_Default value_: `''`)
- `reproducible`: if true, use a slow, but reproducible implementation of sparse sign embeddings. (_Default value_: `false`)

Inputting a value of `[]` for an optional argument results in the default value.

### Iterative sketching

_Iterative sketching_ is a fast, forward stable randomized least-squares solver.
As with FOSSILS, begin by sketching and QR-factorizing $SA = QR$.
We use a sparse sign embedding for the embedding matrix $S$.
After which, iterative sketching produces a sequence of better-and-better least-squares solutions using the iteration

$$
x_{i+1} = x_i + R^{-1} ( R^{-\top} (A^\top(b - Ax_i))).
$$

The main result of [_Fast and forward stable randomized algorithms for linear least-squares problems_](https://arxiv.org/abs/2311.04362) is that iterative sketching is [forward stable](https://nhigham.com/2020/08/04/what-is-numerical-stability/): If you run it for sufficiently many iterations, the forward error $\lVert x - x_i \rVert$ and residual error $\lVert (b-Ax) - (b-Ax_i) \rVert$ are roughly as small as for a standard direct method like (Householder) QR factorization.

Iterative sketching can optionally be accelerated by incorporating _damping_ and _momentum_, resulting in a modified iteration

$$
x_{i+1} = x_i + \alpha R^{-1} ( R^{-\top} (A^\top(b - Ax_i))) + \beta (x_i - x_{i-1}).
$$

We call $\alpha$ and $\beta$ the _damping parameter_ and _momentum parameter_ respectively.
The optimal damping and momentum parameters were computed in [these](https://web.stanford.edu/~pilanci/papers/IHSMomentum18.pdf) [papers](https://arxiv.org/abs/1911.02675).

To run iterative sketching using our code, the command is

```
[x, stats] = iterative_sketching(A, b, [d, q, summary, verbose, damping, momentum, reproducible])
```

All but the first two inputs are optional.
The optional inputs are as follows:

- `d`: sketching dimension. (_Default value_: see paper)
- `q`: number of iterations (if `q>=1`) or tolerance (if `q<1`). If `q>=1`, iterative sketching will be run for `q` iterations. Otherwise, iterative sketching is run for an adaptive number of steps until the norm change in residual is less than `q*(Anorm * norm(x) + 0.01*Acond*norm(r))`. Here, `Anorm` and `Acond` are estimates of the norm and condition number of `A`. (_Default value_: `eps`)
- `summary`: a function of `x` to be recorded at every iteration. The results will be outputted in the optional output argument `stats`. (_Default value_: None)
- `verbose`: if true, output at every iteration. (_Default value_: `false`)
- `damping`: damping coefficient. To use the optimal value, set `damping` to `'optimal'`. (_Default value_: 1, i.e., no damping).
- `momentum`: momentum coefficient. To use the optimal value, set both `damping` and `momentum` to `'optimal'`. (_Default value_: 0, i.e., no momentum).
- `reproducible`: if true, use a slow, but reproducible implementation of sparse sign embeddings. (_Default value_: `false`)

Inputting a value of `[]` for an optional argument results in the default value.

### Sketch-and-precondition

Sketch-and-precondition is also based on a QR factorization $SA = QR$ of a sketch of the matrix $A$.
It then uses $R$ as a preconditioner for solving $Ax = b$ using the [LSQR](https://stanford.edu/group/SOL/software/lsqr/lsqr-toms82a.pdf) or [CGNE](https://en.wikipedia.org/wiki/Conjugate_gradient_method#Conjugate_gradient_on_the_normal_equations) method.
To call sketch-and-precondition, use the following command

```
[x, stats] = sketch_and_precondition(A, b, [d, q, summary, verbose, opts, reproducible])
```

All but the first two inputs are optional.
The optional inputs are as follows:

- `d`: sketching dimension. (_Default value_: `2*size(A,2)`)
- `q`: number of iterations. (_Default value_: `100`)
- `summary`: a function of `x` to be recorded at every iteration. The results will be outputted in the optional output argument `stats`. (_Default value_: None)
- `verbose`: if true, output at every iteration. (_Default value_: `false`)
- `opts`: specifies the initial iterate $x_0$ and iterative method (LSQR, symmetric conjugate gradient, or preconditioned conjugate gradient). If `'cold'` is a substring of `opts`, then the initial iterate is chosen to be $x_0 = 0$. Otherwise, we use a warm start and choose $x_0$ to be the [sketch-and-solve solution](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS3). If `cgne` is a substring of `opts`, then we solve $Ax = b$ using CGNE, and if `sym` is a substring of `opts`, then we perform conjugate gradient on the symmetrically preconditioned system $(R^{-\top} A^\top A R^{-1}) y = c$; otherwise, we use LSQR. We recommend using the LSQR or symmetric conjugate gradient implementations in practice. (_Default value_: `''`)
- `reproducible`: if true, use a slow, but reproducible implementation of sparse sign embeddings. (_Default value_: `false`)

Inputting a value of `[]` for an optional argument results in the default value.

### Computing or estimating the backward error

Our code also provides functions to compute or estimate the backward error, defined to be $$\mathrm{BE}_\theta(\hat{x}) = \min \left\\{ \lVert [\Delta A,\theta\cdot\Delta b]\rVert_F : \hat{x} = \mathrm{argmin}_y \lVert (b+\Delta b) - (A+\Delta A)y \rVert \right\\}.$$The backward error is a quantitative measure of the size of the perturbations to $A$ and $b$ needed to make the numerically computed solution $\hat{x}$ the exact solution to the least-squares problem. The parameter $\theta \in [0,\infty]$ sets the relative importances of perturbations to $A$ and $b$. If $\lVert A \rVert = \lVert b \rVert = 1$, a method is backward stable if and only if $\mathrm{BE}_1(\hat{x})$ is at most a small multiple of the machine precision ($\approx 10^{-16}$ in double precision).
See [this section of our paper](https://arxiv.org/html/2406.03468v1#S3.SS3) for details on the backward error and descriptions of the methods below for estimating it.

The backward error can be computed in $O(m^3)$ operations using the Waldén–Karlson–Sun–Higham formula. 
In our code,

```
backward_error_ls(A,b,xhat,[theta])
```

The default value of `theta` is infinity (`Inf`).

To obtain a cheaper estimate of the backward error, one can use the Karlson–Waldén estimate, which is guaranteed to be within a factor of $\sqrt{2}$ of the true backward error (neglecting rounding errors incurred in evaluating the formulas). 
The Karlson–Waldén estimate can be called in our code using the command

```
kw_estimate(A,b,xhat,[theta,S,V])
```

The default value of `theta` is infinity (`Inf`).
The Karlson–Waldén estimate runs in $O(mn^2)$ operations.
If one has access to an SVD of `A` (i.e., `[U,S,V] = svd(A,"econ")`), then the `S` and `V` matrices can be provided to `kw_estimate` as optional arguments, reducing the runtime to $O(mn)$.

Finally, for an even faster estimate of the backward error, one can use the _sketched_ Karlson–Waldén estimate.
The sketched Karlson–Waldén estimate is also within a small constant factor of the true backward error, and runs in a faster (roughly) $O(mn + n^3)$ operations.
The sketched Karlson–Waldén estimate can be called as

```
kw_estimate(A,b,xhat,theta,"sketched",[d])
```

The parameter `d` sets the embedding dimension for the sketch, defaulting to `2*size(A,2)`.
