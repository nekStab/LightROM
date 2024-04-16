---
project: LightROM
summary: Lightweight implementation of reduced-order modeling of large-scale *input-output linear time invariant dynamical systems*.
author: Jean-Christophe Loiseau and Simon Kern
media_dir: ./imgs
graph: true
license: bsd
project_github: https://github.com/nekStab/LightROM
---

@note
This documentation still is work in progress.
@endnote

This toolbox is aimed at reduced-order modeling and control of large-scale linear time-invariant (LTI) dynamics system for which the system matrix \(\mathbf{A}\) is not explicitly available but is instead defined implicitly via a function performing the matrix-vector product \(\mathbf{A}x\). 

Extending on the core utilities provided in the base package `LightKrylov`, from which `LightROM` derived, this package contains lightweight Fortran implementations based on Krylov methods for reduced-order modeling and control for LTI systems, including:

1. Solver for the differential Lyapunov equation (DALE)
   $$ \dot{\mathbf{X}} = \mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{A}^T + \mathbf{B} \mathbf{B}^T $$

2. Solver for the differential Riccati equation (DARE)
   $$ \dot{\mathbf{X}} = \mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{A}^T + \mathbf{C}^T \mathbf{Q_c} \mathbf{C} - \mathbf{X} \mathbf{B} \mathbf{R}^{-1} \mathbf{B}^T \mathbf{X} $$

For reducted-order modeling and control, we are typically interested in the solutions of the algebraic versions of the Lyapunov and Riccati equations (ALE and ARE, respectively), which correspond to the steady-state solution of the DALE and DARE equations (which exist for stable systems). The solvers implemented in `LightROM` compute low-rank approximations of the solutions which is feasible even in cases when the system is too large for direct solvers.

## Capabilities

`LightROM` is based on `LightKrylov` and thus inherits its ecosystem of abstract types which are extended to the use case in the same way. Additionally, `LightROM` implements the abstract type `abstract_lti_system` corresponding to the standard system representation
   $$ \begin{align} \dot{x} &= \mathbf{A} x + \mathbf{B} u \\ y & = \mathbf{C}^T x + \mathbf{D} u \end{align} $$
where \(\mathbf{A}\) is an `abstract_linop`, \(\mathbf{C}^T, \mathbf{B}\) are instances of the `abstract_vector` type and \(\mathbf{D}\) is a real matrix. `LightROM` then provides the following functionalities:

- DALE-solver based on the dynamic low-rank approximation
- DARE-solver based on the dyanmic low-rank approximation

Both solvers employ a double operator splitting for the integration that takes advantage of the positive-definiteness of the solution and can be run using first and second order precision in time via Lie and Strang splittings, respectively.

### Known limitations

The solvers in `LightROM` are only efficient for problems that are low-rank, i.e. \( \mathbf{C}^T \) and \( \mathbf{B} \) are low-rank.

## Installation

Provided you have `git` installed, getting the code is as simple as:

```
git clone https://github.com/nekStab/LightROM
```

Alternatively, using `gh-cli`, you can type

```
gh repo clone nekStab/LightROM
```

### Dependencies

`LightROM` is build upon `LightKrylov`. 

- [`LightKrylov`](https://github.com/nekStab/LightKrylov)

Apart from this intrinsic dependency, the toolbox depends only on:

- a Fortran compiler,
- [LAPACK]() (or similar),
- [`fpm`](https://github.com/fortran-lang/fpm) or `make` for building the code.

And that's all of it. To date, the tested compilers include:

- `gfortran 12.0.3`

### Building with `fpm`

Provided you have cloned the repo, installing `LightROM` with `fpm` is as simple as

```
fpm build --profile release
```

Please refer to the `fpm` documentation to see how to link `LightROM` with your application.

### Building with `make`

N/A

### Running the tests

To see if the library has been compiled correctly, a set of unit tests are provided in [test](). If you use `fpm`, running these tests is as simple as

```
fpm test
```

If everything went fine, you should see

```
All tests successfully passed!
```

If not, please feel free to open an Issue.

### Running the examples

To run the examples:

```
fpm run --example
```

This command will run all of the examples sequentially. You can alternatively run a specific example using e.g.

```
fpm run --example DLRA_laplacian2D_lti_lyapunov
```

For more details, please refer to each of the examples.

### Documentation

Online documentation is available [here](https://nekstab.github.io/LightROM/).
If you want to generate the documentation locally, you can do so by using [`ford`](https://github.com/Fortran-FOSS-Programmers/ford), an automatic documentation generator for modern Fortran programs.
Provided you have `ford` installed, you can build the documentation locally by running

```{bash}
ford project-file.md
```

Using Github Actions, the online documentation is automatically updated on every pull request to the `main` branch.
Documentation related to features included in the `dev` branch is only available in the source code.

## Contributing

### Current developers

`LightROM` is currently developed and maintained by a team of three:

- [Jean-Christophe Loiseau](https://loiseaujc.github.io/) : Assistant Professor of Applied maths and Fluid dynamics at [DynFluid](https://dynfluid.ensam.eu/), Arts et Métiers Institute of Technology, Paris, France.
- [Simon Kern](https://github.com/Simkern/) : PhD in Fluid dynamics (KTH, Sweden, 2023) and currently postdoctoral researcher at DynFluid.

Anyone else interested in contributing is obviously most welcome!

## Acknowledgment

The development of `LightROM` is part of an on-going research project funded by [Agence Nationale pour la Recherche](https://anr.fr/en/) (ANR) under the grant agreement ANR-22-CE46-0008.

### Related projects

`LightROM` is based on `LightKrylov`, the base package of our ecosystem. If you like it, you may also be interested in :

- [`neklab`]() : a bifurcation and stability analysis toolbox based on `LightKrylov` for the massively parallel spectral element solver [`Nek5000`]().
