# Quasi-normal modes of non-commutative geometry-inspired dirty black holes

[![Proc. R. Soc. A](https://img.shields.io/badge/Proc.%20R.%20Soc.%20A-481%20%282025%29%2020250021-005a9c)](https://royalsocietypublishing.org/doi/10.1098/rspa.2025.0021)
[![DOI](https://img.shields.io/badge/DOI-10.1098%2Frspa.2025.0021-1a7f37)](https://doi.org/10.1098/rspa.2025.0021)
[![arXiv](https://img.shields.io/badge/arXiv-2507.19107-b31b1b)](https://arxiv.org/abs/2507.19107)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-4c1)](LICENSE)
[![Maple](https://img.shields.io/badge/Maple-2022.1%2B-ff6f00)](https://www.maplesoft.com/products/Maple/)
[![MATLAB](https://img.shields.io/badge/MATLAB-Advanpix%20MCT-0076a8)](https://www.advanpix.com/)
[![QNMs Hall of Fame](https://img.shields.io/badge/QNMs-Hall%20of%20Fame-8a2be2)](https://qnms.denys-dutykh.com/)
[![GitHub](https://img.shields.io/badge/GitHub-dutykh%2Fdirtybh-181717?logo=github)](https://github.com/dutykh/dirtybh)
[![Last commit](https://img.shields.io/github/last-commit/dutykh/dirtybh)](https://github.com/dutykh/dirtybh/commits/main)

**Computational supplemental material**

**Authors**: Dr. Davide Batic, Dr. Denys Dutykh and Zeinabou Ahmed Babou
(Mathematics Department, Khalifa University of Science and Technology,
Abu Dhabi, UAE)

**Repository**: <https://github.com/dutykh/dirtybh>

**Numerical method**: collocation-type spectral method based on Chebyshev
polynomials, coupled with a polynomial (quadratic) eigenvalue solver

**Programming languages**: Maple&trade; and MATLAB&trade;

This repository contains the computational supplemental material for the
following peer-reviewed article:

> D. Batic, D. Dutykh and Z. A. Babou, *Quasi-normal modes of non-commutative
> geometry-inspired dirty black holes*, **Proc. R. Soc. A** 481 (2318):
> 20250021 (2025). Open Access.
> DOI: [10.1098/rspa.2025.0021](https://doi.org/10.1098/rspa.2025.0021).
> Preprint: [arXiv:2507.19107](https://arxiv.org/abs/2507.19107) [gr-qc].

---

## 1. Overview

The codes collected here compute the quasi-normal mode (QNM) spectrum of a
non-commutative geometry-inspired *dirty* black hole, for scalar,
electromagnetic and gravitational (axial) perturbations, in both the
non-extremal and the extremal configuration.

The computation is split into a symbolic stage and a numerical stage. Maple
carries out the derivation of the master equation, the change of variables that
compactifies the black hole exterior onto a finite interval, the ansatz that
encodes the quasi-normal boundary conditions, the classification of the
singularities of the resulting ordinary differential equation, and finally the
assembly of the spectral matrices. MATLAB, driven by the Advanpix
Multiprecision Computing Toolbox, then solves the resulting quadratic
eigenvalue problem in arbitrary precision and produces the convergence diagram.

The abstract of the article reads:

> We investigate the quasinormal modes (QNMs) of noncommutative
> geometry-inspired dirty black holes, focusing on both non-extremal and
> extremal configurations. These gravitational objects, characterized by
> smeared energy distributions within a modified de Sitter-like equation of
> state, modify the classical Schwarzschild metric and regularize central
> singularities. We employ a spectral method based on Chebyshev polynomials to
> solve the eigenvalue problem for scalar, electromagnetic, and gravitational
> perturbations. Our results reveal new overdamped modes indicative of rapid
> decay without oscillation, particularly prominent in near-extremal and
> extremal regimes. Additionally, we establish that the QNMs converge to
> classical Schwarzschild values for large mass parameters, validating our
> method's robustness. Our findings highlight the impact of dirtiness and
> noncommutative effects on black hole QNM spectra, offering potential
> observational signatures for distinguishing these objects in
> gravitational-wave detections.

---

## 2. The physical model

### 2.1 Metric

The background is a static, spherically symmetric *dirty* black hole, that is,
a black hole surrounded by a non-trivial matter distribution. Its line element
carries two independent metric functions,

$$
\mathrm{d}s^{2} = - f(r) e^{-2g(r)} \mathrm{d}t^{2} + \frac{\mathrm{d}r^{2}}{f(r)} + r^{2} \bigl(\mathrm{d}\vartheta^{2} + \sin^{2}\vartheta \mathrm{d}\varphi^{2}\bigr).
$$

The function $f$ is fixed by the enclosed mass, $f(r) = 1 - 2m(r)/r$, while the
factor $e^{-2g(r)}$ measures the *dirtiness*: it is identically one for a clean
vacuum black hole and departs from one wherever matter is present. In the
Schwarzschild case $g \equiv 0$ and $m \equiv M$.

### 2.2 Non-commutative smeared source

Non-commutativity of spacetime coordinates is implemented, following the
standard prescription, by replacing the point-like source with a Gaussian
smeared energy density of minimal width $\sqrt{\theta}$,

$$
\rho(r) = \frac{M}{(4\pi\theta)^{3/2}} \exp\left(-\frac{r^{2}}{4\theta}\right),
$$

where $\theta$ is the non-commutativity parameter (an area). Integrating the
density gives the smeared mass function, which regularizes the central
singularity.

### 2.3 Dimensionless variables used throughout the codes

Every routine in this repository works with the two dimensionless quantities

$$
x = \frac{r}{2M}, \qquad \mu = \frac{M}{\sqrt{\theta}},
$$

so that $x = 1$ is the Schwarzschild radius and $\mu$ is the dimensionless
black hole mass. In these variables the two metric functions read

$$
f(x) = 1 - \frac{\mathrm{erf}(\mu x)}{x} + \frac{2\mu}{\sqrt{\pi}} e^{-\mu^{2}x^{2}},
$$

$$
g(x) = \frac{\mu}{2} \left[ 1 - \mathrm{erf}(\mu x) + \frac{2\mu x}{\sqrt{\pi}} e^{-\mu^{2}x^{2}}\right] = \frac{\mu}{2}\left(1 - \frac{m(r)}{M}\right).
$$

The second form makes the physical content transparent: the dirtiness is
sourced by the fraction of the total mass that has *not* yet been enclosed at
radius $r$. Consequently $g \to 0$ as $x \to \infty$, so the spacetime is
asymptotically flat, while $g \to \mu/2$ at the centre.

These are exactly the expressions coded in `matrixassembler_nonext.mpl` and
`matrixassembler_ext.mpl`.

### 2.4 Horizon structure: non-extremal versus extremal

The horizon is the root of $f$, and it does not involve $g$, so the horizon
structure coincides with that of the non-commutative Schwarzschild geometry.
The number of roots depends on $\mu$:

* **Non-extremal case.** For $\mu$ above the critical value there are two
  horizons; the codes track the outer one, $x_h$, obtained from
  `fsolve(f(x) = 0, x = 0.98)`. As $\mu$ grows, $x_h \to 1$ and the classical
  Schwarzschild results are recovered.
* **Extremal case.** At the critical mass the two horizons merge into a
  degenerate one. It is found by solving the double-root conditions
  $f = 0$ and $\partial f / \partial x = 0$ simultaneously, which yields

  $\mu_{e} \simeq 1.90$ and $x_{e} \simeq 0.79$, that is
  $M \simeq 1.90\sqrt{\theta}$ and $r_{e} \simeq 3.0\sqrt{\theta}$, the
  well-known critical values of the non-commutative Schwarzschild geometry.
  Below $\mu_{e}$ no horizon forms.

A useful identity, used implicitly by the codes, is that $f(r_h) = 0$ implies
$m(r_h) = M x_h$, whence $g(r_h) = \tfrac{\mu}{2} (1 - x_h)$.

The distinction matters numerically, not just physically: at a non-degenerate
horizon the tortoise coordinate diverges logarithmically, whereas at a
degenerate horizon it diverges as a simple pole. The two cases therefore
require different ansatzes, and this is why the repository carries two parallel
sets of files throughout.

---

## 3. Perturbation equations

Linear perturbations of spin $s$ separate, and the radial part obeys a
Regge-Wheeler-type master equation

$$
\frac{\mathrm{d}^{2}\psi}{\mathrm{d}r_{\ast}^{2}} + \bigl[\omega^{2} - V_{s}(r)\bigr] \psi = 0, \qquad \frac{\mathrm{d}r_{\ast}}{\mathrm{d}r} = \frac{e^{g(r)}}{f(r)},
$$

with the effective potential for a dirty background

$$
V_{s}(r) = f(r) e^{-2g(r)} \left[\frac{\ell(\ell+1)}{r^{2}} + (1-s^{2}) \frac{f'(r)}{r} - (1-s) \frac{f(r) g'(r)}{r}\right].
$$

The three physically relevant spins are covered:

| $s$ | Perturbation | Potential |
|:---:|:---|:---|
| 0 | massless scalar | $f e^{-2g}\bigl[\ell(\ell+1)/r^{2} + f'/r - f g'/r\bigr]$ |
| 1 | electromagnetic | $f e^{-2g} \ell(\ell+1)/r^{2}$ (purely centrifugal) |
| 2 | gravitational, axial | $f e^{-2g}\bigl[\ell(\ell+1)/r^{2} - 3f'/r + f g'/r\bigr]$ |

Setting $g \equiv 0$ and $f = 1 - 2M/r$ recovers the classical Schwarzschild
potentials, namely $+2M/r^{3}$ for the scalar field and $-6M/r^{3}$ for the
Regge-Wheeler equation. This limit is a useful sanity check when modifying the
codes. The angular momentum must satisfy $\ell \geqslant s$, and the assembler
routines enforce this.

### Quasi-normal boundary conditions

Quasi-normal modes are the solutions that are purely ingoing at the horizon and
purely outgoing at spatial infinity,

$$
\psi \sim e^{-i\omega r_{\ast}} \quad (r \to r_{h}), \qquad \psi \sim e^{+i\omega r_{\ast}} \quad (r \to \infty).
$$

The admissible frequencies $\omega$ form a discrete set; they are complex, and
$\mathrm{Im}\omega \lt 0$ signals damping. Purely imaginary frequencies
correspond to overdamped, non-oscillatory modes, which the article finds to be
prominent in the near-extremal and extremal regimes.

---

## 4. Numerical method

### 4.1 Compactification

The black hole exterior $x \in [x_h, \infty)$ is mapped onto the Chebyshev
interval $y \in [-1, 1]$ by

$$
x = \frac{2 x_{h}}{1-y} \qquad\Longleftrightarrow\qquad y = 1 - \frac{2 x_{h}}{x},
$$

so that $y = -1$ is the event horizon and $y = +1$ is spatial infinity. The
same map with $x_h$ replaced by $x_e$ is used in the extremal case.

### 4.2 Encoding the boundary conditions

Rather than imposing the boundary conditions as extra algebraic constraints,
the divergent oscillatory behaviour at both ends is factored out analytically.
The remaining amplitude is regular on the closed interval, so it can be expanded
in Chebyshev polynomials with no boundary rows at all. Only solutions with the
correct asymptotics survive, and every computed eigenvalue is a genuine QNM.

The factorization requires the near-horizon scale of the tortoise coordinate,
and this is where the two cases part company.

**Non-extremal case.** The relevant constant is

$$
\alpha = \frac{e^{-\mu(x_h-1)/2}} {1 - 4\mu^{3}x_{h}^{2}e^{-\mu^{2}x_{h}^{2}}/\sqrt{\pi}} = \frac{1}{2 x_{h} \kappa}, \qquad \kappa = \tfrac{1}{2} f'(r_h) e^{-g(r_h)},
$$

with $\kappa$ the surface gravity. The logarithmic divergence of $r_{\ast}$ is
governed by $1/(2\kappa)$.

**Extremal case.** Here $f'(r_e) = 0$ and $\kappa = 0$, so $\alpha$ is
meaningless. The code instead introduces an auxiliary function $\eta(y)$ that
carries an explicit $1/(1+y)$ pole, together with its first two derivatives and
a constant $a$ built from the double-root data. This reproduces the power-law
divergence of the tortoise coordinate at a degenerate horizon.

### 4.3 Spectral discretisation

The regular amplitude is expanded in a truncated Chebyshev series,

$$
P(y) = \sum_{j=0}^{n-1} a_{j} T_{j}(y),
$$

and the transformed ordinary differential equation is collocated at the $n$
Chebyshev roots

$$
y_{i} = \cos\left(\frac{(2i-1)\pi}{2n}\right), \qquad i = 1,\dots,n .
$$

The transformed equation is quadratic in the frequency and second order in $y$,

$$
\sum_{k=0}^{2} \sum_{m=0}^{2} L_{km}(y) \Omega^{k} \frac{\mathrm{d}^{m}P}{\mathrm{d}y^{m}} = 0,
$$

with $L_{12} = L_{21} = L_{22} = 0$, so that the second derivative appears only
at order $\Omega^{0}$ and the term of order $\Omega^{2}$ multiplies $P$ alone.
Collecting the coefficients of each $a_j$ at each collocation point produces
three real $n \times n$ matrices $M_0$, $M_1$, $M_2$.

### 4.4 The quadratic eigenvalue problem

The discrete problem is the quadratic (polynomial) eigenvalue problem

$$
\bigl[ M_{0} + i \omega M_{1} + \omega^{2} M_{2} \bigr] \mathbf{a} = \mathbf{0},
$$

solved in MATLAB with `polyeig`. The three matrices are real; the imaginary
unit is attached to the linear term at solve time, which is why the Maple side
exports $M_1$ without it.

### 4.5 Why multiprecision is unavoidable

The Chebyshev coefficients of a function with endpoint singularities of the
type occurring here span many orders of magnitude, and the quadratic eigenvalue
problem is severely ill-conditioned: in double precision the computed spectrum
is dominated by round-off well before the physical modes are resolved. The
codes therefore run at a working precision matched to the truncation order, by
default `Digits` and `mp.Digits` equal to $n$ itself. Convergence is assessed
by superposing the spectra obtained for several values of $n$: the physical
QNMs are the points that do not move.

---

## 5. Repository structure

```
dirtybh/
├── 01-NONEXTREME CASE_NCGDBHs-SUPPLEMENTAL MATERIAL.mw   Maple worksheet: full symbolic derivation, non-extremal case
├── 02-EXTREME CASE_NCDBHs-SUPPLEMENTAL MATERIAL.mw       Maple worksheet: full symbolic derivation, extremal case
├── Effective_Potential_in_r_non_extreme.mw               Maple worksheet: effective potential in r, non-extremal case
├── Effective_Potential_in_r_extreme.mw                   Maple worksheet: effective potential in r, extremal case
├── matrixassembler_nonext.mpl                            Maple procedure: assembles M0, M1, M2, non-extremal case
├── matrixassembler_ext.mpl                               Maple procedure: assembles M0, M1, M2, extremal case
├── chaseQNMs.m                                           MATLAB driver: solves the eigenvalue problem and plots
├── LICENSE                                               GNU Lesser General Public License, version 2.1
└── README.md                                             This file
```

The matrices exchanged between Maple and MATLAB live in a `data/`
subdirectory, which is not tracked by git and must be created before the first
run (see [Section 7](#7-how-to-run)).

---

## 6. The codes, file by file

### 6.1 `01-NONEXTREME CASE_NCGDBHs-SUPPLEMENTAL MATERIAL.mw`

Maple 2022.1 worksheet. Complete symbolic derivation for the non-extremal
configuration, annotated step by step with the corresponding equation numbers
of the manuscript. It proceeds through:

1. definition of the master perturbation equation;
2. introduction of the quasi-normal ansatz;
3. derivation of the transformed equation, and an independent verification that
   the transformed and original forms are equivalent;
4. the change of independent variable, performed with `PDEtools[dchange]`,
   onto the compactified interval;
5. classification of the points $y = \pm 1$ for the coefficient functions
   $F(y)$, $V(y)$, $S_2(y)$, $S_1(y)$, $S_0(y)$: the worksheet establishes, for
   instance, that $S_1$ and $S_0$ each have a pole of order 2 at $y = 1$ and a
   zero of order 3 at $y = -1$. These results constitute the singularity table
   of the article;
6. definition of $M_0(y)$, $N_1(y)$, $N_0(y)$ and $C_2(y)$, $C_1(y)$, $C_0(y)$,
   and evaluation of the limits that define the coefficients $\Lambda_1$,
   $\Lambda_0$ and $B_2$, $B_1$, $B_0$, each checked against the closed forms
   quoted in the article;
7. after an explicit `restart`, definition of the entries
   `L00hat` through `L22hat` that the numerical assembler
   implements.

**Note.** The worksheet contains a deliberate `restart` partway through. The
second half must be executed from that fresh start, as the annotations in the
worksheet indicate.

### 6.2 `02-EXTREME CASE_NCDBHs-SUPPLEMENTAL MATERIAL.mw`

Maple 2022.1 worksheet. The counterpart of the above for the extremal
configuration, following the same seven steps with the extremal ansatz. The
singularity structure differs, as the worksheet documents explicitly: $S_{2e}$
is finite at $y = 1$ and has a zero of order 6 at $y = -1$, while $S_{1e}$ and
$S_{0e}$ each have a pole of order 2 at $y = 1$ and a zero of order 4 at
$y = -1$. It also rewrites $M_{2e}$, $M_{1e}$, $M_{0e}$ as polynomials in
$\Omega$, which is the form required by the spectral method, and computes the
limits collected in the corresponding table of the article. It likewise
contains an explicit `restart` before the definition of
`L00ehat` through `L22ehat`.

### 6.3 `Effective_Potential_in_r_non_extreme.mw` and `Effective_Potential_in_r_extreme.mw`

Maple 18 worksheets. They construct the effective potential directly in the
radial variable $r$, compare it against the potential of the non-commutative
Schwarzschild geometry (denoted `UNCSCHW` in the source), and perform the change
of variables to the compactified coordinate. They are self-contained and serve
to validate the potential used by the assembler routines; they are not part of
the eigenvalue pipeline.

These two worksheets were written under an older Maple release than the two
above, and open in any later version.

### 6.4 `matrixassembler_nonext.mpl`

Defines the Maple procedure

```maple
MatrixAssembler(d, n, s, L, mu, p)
```

| Argument | Type | Meaning |
|:---|:---|:---|
| `d` | integer | number of digits used in the computation |
| `n` | integer | number of Chebyshev modes, and the size of the matrices |
| `s` | numeric | spin: 0 (scalar), 1 (electromagnetic), 2 (gravitational) |
| `L` | numeric | angular momentum $\ell$, must satisfy $\ell \geqslant s$ |
| `mu` | numeric | dimensionless black hole mass $\mu = M/\sqrt{\theta}$ |
| `p` | string | path under which the assembled matrices are written |

The procedure locates the event horizon, prints the value of $x_h$ it found,
builds $\alpha$, $f$, $f'$, $g$, $g'$, $F = f e^{-g}$ and the effective
potential $V$, assembles the nine coefficient functions $L_{km}$, evaluates
them at the Chebyshev roots, extracts the coefficient of each $a_j$, and writes
the three matrices to

```
<p>/data/M0_<n>.mat
<p>/data/M1_<n>.mat
<p>/data/M2_<n>.mat
```

in MATLAB ASCII format via `ExportMatrix`. An error is raised if $\ell \lt s$.

### 6.5 `matrixassembler_ext.mpl`

The extremal counterpart, with the signature

```maple
MatrixAssembler(d, n, s, L, p)
```

The argument `mu` is absent: in the extremal case the mass is not free. The
procedure solves the double-root system for $x_e$ and $\mu$ simultaneously and
prints both values it found. It then builds the constant $a$ and the auxiliary
functions $\eta$, $\eta'$, $\eta''$ described in Section 4.2, assembles the
same nine coefficient functions in their extremal form, and exports the three
matrices under exactly the same file names as the non-extremal routine.

Because the file names coincide, the extremal and non-extremal matrices must be
kept in separate `data/` directories, or generated and consumed one case at a
time.

### 6.6 `chaseQNMs.m`

The MATLAB driver. It loops over a list of truncation orders, by default

```matlab
list = [200, 150, 100];
```

and for each order sets the working precision to match, reads the three
matrices produced by Maple, and solves the quadratic eigenvalue problem:

```matlab
mp.Digits(n);
M0 = mp.read(strcat('data/M0_', nstr, '.mat'));
M1 = mp.read(strcat('data/M1_', nstr, '.mat'));
M2 = mp.read(strcat('data/M2_', nstr, '.mat'));
e  = polyeig(M0, mp('1i')*M1, M2);
```

It finally overlays the spectra for all truncation orders in the complex
frequency plane, $\mathrm{Re}\omega$ against $\mathrm{Im}\omega$,
with one marker style per order. Eigenvalues that coincide across the three
resolutions are the converged quasi-normal modes; the rest are numerical
artefacts of the truncation. This plot is the convergence diagram of the
article.

One line must be edited before the first run:

```matlab
addpath('/path/to/advanpix/');
```

It has to point at the local installation of the Advanpix toolbox.

---

## 7. How to run

The pipeline runs Maple first, then MATLAB.

**Step 1. Create the exchange directory.** From the repository root:

```bash
mkdir -p data
```

This is required: `ExportMatrix` writes into `<p>/data/` but does not create
the directory.

**Step 2. Assemble the matrices in Maple.** For the non-extremal case, load the
assembler and call it once per truncation order, keeping `s`, `L` and `mu`
fixed:

```maple
read "matrixassembler_nonext.mpl":
MatrixAssembler(200, 100, 0, 2, 5.0, "/absolute/path/to/dirtybh"):
MatrixAssembler(200, 150, 0, 2, 5.0, "/absolute/path/to/dirtybh"):
MatrixAssembler(200, 200, 0, 2, 5.0, "/absolute/path/to/dirtybh"):
```

The example computes the scalar ($s = 0$), $\ell = 2$ spectrum at $\mu = 5$.
For the extremal case the call drops the mass argument, which the routine
determines itself:

```maple
read "matrixassembler_ext.mpl":
MatrixAssembler(200, 100, 0, 2, "/absolute/path/to/dirtybh"):
MatrixAssembler(200, 150, 0, 2, "/absolute/path/to/dirtybh"):
MatrixAssembler(200, 200, 0, 2, "/absolute/path/to/dirtybh"):
```

Choose `d` at least as large as the highest `mp.Digits` that MATLAB will
request, since `d` fixes the number of digits actually written to the ASCII
files. With the default driver this means `d` of at least 200.

**Step 3. Solve and plot in MATLAB.** Edit the `addpath` line in
`chaseQNMs.m`, make sure `list` matches the truncation orders assembled in
Step 2, then run

```matlab
>> chaseQNMs
```

The script reports its progress order by order and produces the complex-plane
convergence diagram.

**Changing the physical parameters.** The spin and angular momentum are
arguments of `MatrixAssembler`, so a different perturbation type or multipole
requires re-running Step 2. The same applies to the mass $\mu$ in the
non-extremal case. Nothing in `chaseQNMs.m` needs to change other than `list`.

---

## 8. Requirements

| Component | Requirement |
|:---|:---|
| Maple | version 2022.1 or later for `matrixassembler_*.mpl` and the two supplemental worksheets; the two effective-potential worksheets were written under Maple 18 and open in any later release |
| MATLAB | any reasonably recent release providing `polyeig` |
| MATLAB toolbox | [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/), mandatory: the `mp` class, `mp.Digits`, `mp.read` and multiprecision `polyeig` are all used |

No other third-party package is needed, and no compilation step is involved.

**On cost.** Both stages are expensive, and the cost grows sharply with the
truncation order `n`. The Maple stage evaluates the coefficient functions
symbolically at every one of the `n` collocation points and extracts `n`
coefficients at each, so its work grows quadratically in `n` on top of a high
constant from multiprecision special-function evaluation. The MATLAB stage
solves a dense quadratic eigenvalue problem of size `n` at several hundred
decimal digits, which linearises to a dense $2n \times 2n$ problem. At
$n = 200$ and 200 digits both stages take a substantial time, and it is worth
generating the smaller orders first to confirm that the parameters are sensible
before committing to the largest one.

---

## 9. Reproducing the results of the article

* **The tables of quasi-normal frequencies.** Each table corresponds to one
  choice of the perturbation spin $s$, the multipole number $\ell$ and, in the
  non-extremal case, the dimensionless mass $\mu$. Set those three as arguments
  of `MatrixAssembler` in Step 2, then read off the converged eigenvalues in
  Step 3. The overtone index $n$ orders the converged modes by increasing
  $\lvert \mathrm{Im}\omega \rvert$; it is not an input of the computation, since the
  spectral method returns the whole spectrum at once.
* **The large-mass limit.** Increasing $\mu$ drives the non-extremal spectrum
  towards the classical Schwarzschild values, which is the validation test
  reported in the article.
* **The extremal and near-extremal spectra.** Use `matrixassembler_ext.mpl` for
  the exactly extremal configuration, and `matrixassembler_nonext.mpl` with
  $\mu$ slightly above $\mu_e \simeq 1.90$ for the near-extremal one. The
  purely imaginary, overdamped modes reported in the article appear in this
  regime.
* **The figure.** The convergence diagram is the direct output of
  `chaseQNMs.m` with the default `list = [200, 150, 100]`.

---

## 10. Quasi-Normal Modes Hall of Fame

This article is part of the **Quasi-Normal Modes Hall of Fame**, subtitled
*Celestial Music*:

### <https://qnms.denys-dutykh.com/>

The site is maintained by Dr. Davide Batic and Dr. Denys Dutykh at the
Mathematics Department of Khalifa University of Science and Technology, with
contributions from international collaborators. It gathers the group's
peer-reviewed manuscripts on quasi-normal modes, each entry carrying its full
citation, the DOI of the published version, a link to the arXiv preprint, and a
link to the companion repository holding the computational material.

The present repository is the companion code for the dirty black hole entry.
Readers looking for the group's other quasi-normal mode studies, and for the
codes accompanying them, will find them indexed there.

---

## 11. How to cite

If you use these codes, in whole or in part, please acknowledge this work by
citing the article:

> D. Batic, D. Dutykh and Z. A. Babou, *Quasi-normal modes of non-commutative
> geometry-inspired dirty black holes*, Proc. R. Soc. A **481** (2318):
> 20250021 (2025). DOI: 10.1098/rspa.2025.0021.

In BibTeX:

```bibtex
@article{Batic2025,
  author  = {Batic, Davide and Dutykh, Denys and Babou, Zeinabou Ahmed},
  title   = {Quasi-normal modes of non-commutative geometry-inspired
             dirty black holes},
  journal = {Proceedings of the Royal Society A},
  volume  = {481},
  number  = {2318},
  pages   = {20250021},
  year    = {2025},
  doi     = {10.1098/rspa.2025.0021},
  url     = {https://doi.org/10.1098/rspa.2025.0021}
}
```

The preprint, should you need it:

```bibtex
@misc{Batic2025arXiv,
  author       = {Batic, Davide and Dutykh, Denys and Babou, Zeinabou Ahmed},
  title        = {Quasinormal modes of noncommutative geometry-inspired
                  dirty black holes},
  year         = {2025},
  eprint       = {2507.19107},
  archivePrefix= {arXiv},
  primaryClass = {gr-qc},
  url          = {https://arxiv.org/abs/2507.19107}
}
```

---

## 12. Licence

The codes in this repository are distributed under the **GNU Lesser General
Public License, version 2.1**. The full text is in [`LICENSE`](LICENSE).

They are provided in the hope that they will be useful, but without any
warranty, without even the implied warranty of merchantability or fitness for a
particular purpose. See the licence for details.

---

## 13. Authors and contact

* **Dr. Davide Batic**, Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE (corresponding author of the article)
* **Dr. Denys Dutykh**, Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE
* **Zeinabou Ahmed Babou**, Mathematics Department, Khalifa University of
  Science and Technology, Abu Dhabi, UAE

Questions about the codes, bug reports and suggestions are welcome through the
[issue tracker](https://github.com/dutykh/dirtybh/issues) of this repository.
