---
title: Multiplierless FIR Filter Design with Parallel-Cut Ellipsoid Method
bibliography:
  [
    "fir-ref.bib",
    "ellipsoid.bib",
  ]
csl: "applied-mathematics-letters.csl"
abstract: |
  A general-purpose multiplier is the most expensive arithmetic element in a
  digital finite impulse response (FIR) filter, so a filter whose coefficients
  are represented in canonical signed digit (CSD) form can be realized with a
  small number of adders and shifts instead of multipliers. The design problem,
  however, is awkward: the magnitude constraints are not convex in the impulse
  response, the coefficient quantizer is not convex at all, and the constraint
  set is effectively infinite because the frequency response must be bounded at
  every frequency. This article shows that all three difficulties are handled
  naturally by the ellipsoid method. Spectral factorization turns the magnitude
  specifications into convex constraints on the autocorrelation, the two-sided
  magnitude bounds are returned by the oracle as parallel cuts that shrink the
  search ellipsoid faster than a single cut, and the CSD quantizer is embedded
  in the oracle so that the optimizer reasons about a realizable coefficient
  vector rather than a relaxed one. The complete pipeline, from a
  frequency-domain specification to a synthesizable shift-and-add description,
  is described and implemented in three languages. Experiments quantify the
  iteration economy of the parallel-cut formulation --- a reduction of $77$ to
  $97$ percent --- and compare the running time across the three languages over
  four filter orders, and the generated RTL is synthesized to report
  technology-independent cell, lookup-table, and flip-flop counts.
---

## Introduction

A finite impulse response (FIR) filter computes a weighted sum of the recent
input samples,
$$y[t] = \sum_{k=0}^{n-1} h[k]\, u[t-k],$$
and is specified by its frequency response $H(\omega) = \sum_{m=0}^{n-1}
h(m) e^{-j m \omega}$. Its design consists of choosing the coefficients
$h[0], \ldots, h[n-1]$ so that $|H(\omega)|$ stays inside a prescribed mask.
The classical tools for this task --- windowing and frequency sampling
[@oppenheim1989discrete], the Parks--McClellan exchange algorithm
[@park1972chebyshev], METEOR [@steiglitz1992meteor], and
peak-constrained least-squares design [@selesnick1996constrained;
@adams1998peak] --- all assume that a full multiplier is available for each
coefficient.

In application-specific integrated circuits and field-programmable gate arrays
the multiplier is precisely the element one wants to avoid: it dominates area
and power, and it becomes disproportionately expensive as the word length
grows. A *multiplierless* filter replaces each product by a signed sum of
shifted copies of the input. Because a shift by a power of two is free in
dedicated hardware, a coefficient that is a sum of a few signed powers of two
can be realized with a few adders and no multiplier [@george1960csd;
@samueli1989improved; @markovic2012dsp]. The design problem is therefore not
merely to find good coefficients, but to find coefficients that are *cheap to
realize*.

The design of filters whose coefficients are confined to a discrete set has a
long history, and the multiplierless case is its most hardware-driven instance.
The finite-wordlength problem --- choosing coefficients from a fixed grid rather
than from the continuum --- was among the earliest digital-filter design
questions, and the first systematic attacks cast it as an integer program over a
discrete coefficient space [@kodek1980design; @lim1982finite]. Because
enumeration grows exponentially with the number of taps, the literature turned
to branch-and-bound [@lawler1966branch; @kodek1980algorithm] and to local search
seeded from a continuous relaxation [@kodek1981comparison]. Restricting the
coefficients to powers of two, or to a signed sum of a few powers of two, proved
especially attractive in hardware, and specialized search procedures were
developed for it [@lim1983fir; @samueli1989improved; @lim1999signed]. Later work
exploited structure in the coefficient space: a trellis formulation solved the
signed-power-of-two problem by dynamic programming [@chen1999trellis],
metaheuristics such as genetic algorithms were applied where the search space
resisted exact methods [@xu1995design; @glover2003handbook], and lattice
reduction was used to obtain optimal finite-wordlength designs [@kodek2012lll].
A parallel line of work characterized the performance achievable under a
word-length constraint [@kodek2005performance; @kodek1997limits]. These methods
are exact, or nearly so, only under assumptions that are difficult to satisfy in
practice: integer programming becomes prohibitive as the order grows, and the
heuristics return a feasible design only when the search happens to fall into a
feasible basin.

A second line of work keeps the coefficients continuous and reshapes the problem
so that convex optimization applies. Spectral factorization turns the magnitude
specification into a linear constraint on the autocorrelation [@wu1999fir;
@goodman1997spectral], and the resulting convex program is solved by
second-order cone or semidefinite methods, with the spectral mask encoded as a
linear matrix inequality [@davidson2002linear; @tuan2007efficient] and
implemented in disciplined-convex-programming tools such as CVX
[@grant2008cvx]. The survey of [@davidson2010enriching] collects many such
reformulations. Because these methods return a continuous optimum, however, they
answer a different question from the one the hardware asks. The optimum
typically lies on the boundary of the feasible set, where the active magnitude
constraints are tight, so rounding it to the nearest CSD vector can push part of
the response below its lower bound and destroy feasibility
[@kodek1981comparison; @davidson2010enriching]. A relaxation can supply a lower
bound or a warm start, but it cannot by itself certify a realizable design.

The ellipsoid method has meanwhile been applied well beyond linear programming,
and it is the natural algorithm whenever the constraint set is convex but its
description is implicit or infinite. Prior applications of the same principle
include robust analog circuit sizing with affine arithmetic [@liu2007robust],
parametric network optimization through optimal matrix scaling
[@orlin1985computing] and minimum mean-cycle computation [@dasdan1998faster],
and multi-parameter clock-skew scheduling [@zhou2015multi]. What these
applications share, and what multiplierless FIR design inherits, is a moderate
number of design variables combined with a very large number of constraints ---
the regime in which the ellipsoid method is competitive with interior-point
methods [@boyd2009convex; @boyd2008ellipsoid]. The distinguishing feature of the
present article is that it keeps the quantizer inside the oracle, so that the
optimizer reasons about a realizable coefficient vector rather than a relaxed
one.

Two obstacles stand in the way. First, the magnitude constraint
$$L(\omega) \le |H(\omega)| \le U(\omega), \qquad \forall \omega \in [0,\pi],
$$ {#eq:mag_cons}
is not convex in the coefficient vector $\mathbf{h}$; and it must hold at every
frequency, so it is a semi-infinite constraint. Second, the set of coefficients
representable with a bounded number of signed powers of two is discrete, and
the natural procedure of optimizing first and rounding afterwards can destroy
feasibility. This article addresses both by treating the whole design as a
feasibility problem for the ellipsoid method, whose only requirement is a
separation oracle.

The article is organized as follows. @sec:convex reformulates the magnitude
specifications as convex constraints on the autocorrelation via spectral
factorization. @sec:method recalls the ellipsoid method and develops the
parallel-cut update in the form used throughout. @sec:csd introduces the
canonical signed digit representation and the associated adder cost.
@sec:discrete embeds the quantizer in the oracle. @sec:spectral compares two
spectral-factorization algorithms. @sec:experiments describes the
implementation and reports cross-language results, and @sec:conclusion
concludes.

## Deriving a Convex Design Problem {#sec:convex}

Work with the squared magnitude rather than the magnitude. The *autocorrelation*
of the impulse response is
$$r[k] = \sum_{i} h[i]\, h[i+k], \qquad k = 0, 1, \ldots, n-1,
$$ {#eq:autocorr}
and the squared magnitude is its cosine series
$$R(\omega) = |H(\omega)|^{2} = r[0] + 2\sum_{k=1}^{n-1} r[k] \cos(k\omega)
  = \mathbf{a}(\omega)^\mathsf{T} \mathbf{r},
$$ {#eq:power_spectrum}
where $\mathbf{a}(\omega) = (1, 2\cos\omega, \ldots, 2\cos((n-1)\omega))^\mathsf{T}$.
The map $\mathbf{r} \mapsto R(\omega)$ is affine, so the *squared* magnitude
bounds
$$L^{2}(\omega) \le R(\omega) \le U^{2}(\omega), \qquad
  R(\omega) \ge 0, \qquad \forall \omega \in [0,\pi],
$$ {#eq:cvx_design}
define a convex set in $\mathbf{r}$. This is the classical spectral-decomposition
device of [@wu1999fir], building on the spectral factorization of
[@goodman1997spectral]; the same convexification is surveyed in
[@davidson2010enriching].

The problem solved in this article is
$$\begin{array}{ll}
    \text{minimize}   & \gamma, \\
    \text{subject to} & L^{2}(\omega) \le R(\omega) \le U^{2}(\omega),
                        \quad \omega \in \Omega_{p}, \\
                      & R(\omega) \le \gamma, \quad \omega \in \Omega_{s}, \\
                      & R(\omega) \ge 0, \quad \omega \in [0,\pi],
  \end{array}
$$ {#eq:design}
where $\Omega_{p}$ and $\Omega_{s}$ are the passband and stopband, and
$\gamma$ bounds the stopband energy. The optimization variable is the
autocorrelation $\mathbf{r}$; the frequency response is recovered afterwards by
spectral factorization (@sec:spectral). Because the constraints must hold on a
continuum, and because the number of design variables is moderate while the
number of constraints is large, this is exactly the regime in which the
ellipsoid method is competitive with interior-point methods
[@boyd2009convex; @boyd2008ellipsoid].

![A typical structure of an FIR filter.](ellipsoid.files/fir_strctr.pdf){#fig:fir-strctr}

## The Parallel-Cut Ellipsoid Method {#sec:method}

### Separation Oracles and Cutting Planes

Let $\mathcal{K} \subset \mathbb{R}^{n}$ be compact and convex. A *separation
oracle* $\Omega$, queried at a point $x_{0}$, either certifies that
$x_{0} \in \mathcal{K}$ or returns a hyperplane that separates $x_{0}$ from
$\mathcal{K}$:
$$g^\mathsf{T}(x - x_{0}) + \beta \le 0, \qquad \beta \ge 0, \quad g \neq 0,
  \qquad \forall x \in \mathcal{K}.
$$ {#eq:cut}
The pair $(g,\beta)$ is a *cutting plane*: it discards the half-space where
$g^\mathsf{T}(x-x_{0}) + \beta > 0$. If $\beta = 0$ the cut is *central*, if
$\beta > 0$ it is *deep*, and if $\beta < 0$ it is a *shadow cut*. When
$\mathcal{K}$ is defined by $f_{j}(x) \le 0$, a cut is obtained for free from a
subgradient, $(g,\beta) = (\partial f_{j}(x_{0}), f_{j}(x_{0}))$; for a
differentiable $f_{j}$ the subgradient is the gradient. The method was
introduced independently by Shor [@shor1977cutoff] and by Yudin and Nemirovskii
[@yudin1976informational], and used by Khachiyan to prove that linear programming
is polynomial-time solvable [@khachiyan1979polynomial]; see the survey of
[@BGT81; @bland1981ellipsoid].

### The Ellipsoid Update

The search space is the ellipsoid
$$\mathcal{E} = \{\, x \mid (x - x_{c})^\mathsf{T} Q^{-1} (x - x_{c})
  \le \kappa \,\},
$$ {#eq:ellipsoid}
with center $x_{c}$, shape matrix $Q \succ 0$, and scale $\kappa > 0$. Writing
$$\tilde g = Q g, \qquad \omega = g^\mathsf{T} \tilde g, \qquad
  \tau = \sqrt{\kappa\,\omega},$$
every cut type applies the same update
$$x_{c}^{+} = x_{c} - \frac{\rho}{\omega}\tilde g, \qquad
  Q^{+} = Q - \frac{\sigma}{\omega}\tilde g\,\tilde g^\mathsf{T}, \qquad
  \kappa^{+} = \delta\,\kappa,
$$ {#eq:update}
and differs only in the scalar triple $(\rho,\sigma,\delta)$. For a central cut
$(\beta = 0)$,
$$\rho = \frac{\tau}{n+1}, \qquad \sigma = \frac{2}{n+1}, \qquad
  \delta = \frac{n^{2}}{n^{2}-1},
$$ {#eq:central}
and for a deep cut with offset $\beta > 0$,
$$\begin{aligned}
\rho &= \frac{\tau + n\beta}{n+1}, \\
\sigma &= \frac{2(\tau + n\beta)}{(n+1)(\tau+\beta)}, \\
\delta &= \frac{n^{2}}{n^{2}-1}\cdot\frac{\tau^{2}-\beta^{2}}{\tau^{2}}.
\end{aligned}$$
{#eq:deep}
The application of the update to $Q$ alone, rather than to a general positive
definite matrix, saves $n^{2}$ floating-point operations per iteration; the
split into $\kappa$ and $Q$ is what makes the scalar $\tau$ cheap to compute
[@bland1981ellipsoid].

### Parallel Cuts

The magnitude specifications (@eq:mag_cons) are two-sided, and this is where the
method pays off. Whenever a constraint has the form
$$l \le a^\mathsf{T} x + b \le u,$$
the oracle can return a *pair* of parallel cuts that share a normal $g$:
$$\begin{aligned}
g^\mathsf{T}(x - x_{c}) + \beta_{0} &\le 0, \\
g^\mathsf{T}(x - x_{c}) + \beta_{1} &\ge 0, \qquad \forall x \in \mathcal{K}.
\end{aligned}$$
{#eq:parallel}
The pair removes a slab from the ellipsoid rather than a half-space, which
shrinks the volume more than either cut alone and therefore reduces the number
of iterations [@frenk1994deep]. A single deep cut is the special case in which
one of the two planes is tangent.

The geometry is worth stating explicitly. Parameterize the current ellipsoid
along the cut normal by $s = g^\mathsf{T}(x - x_{c})/\tau$, so that $s$ ranges
over $[-1,1]$ on $\mathcal{E}$; the two cuts retain only the band
$-\beta_{1}/\tau \le s \le -\beta_{0}/\tau$. The complement of the band consists
of two caps, one where $s$ is too large and one where it is too small. A single
cutting plane removes only one cap, whereas the pair removes both in the same
iteration, so the least-volume ellipsoid that must contain the survivors is
tighter. For a magnitude specification the two caps have a physical meaning:
they are the frequency regions where the response exceeds the upper bound and
falls below the lower bound, and the parallel cut discards both violations at
once. Two degenerate cases anchor the formulas below. When the band is centered,
$\beta_{0} + \beta_{1} = 0$, the new center coincides with the old one and the
update must remain finite. When one plane is tangent, $\beta_{1} = \tau$, the
band collapses to a half-space and the update must reduce to the deep cut
@eq:deep.

The auxiliary quantities have direct interpretations. Each
$\zeta_{i} = \tau^{2} - \beta_{i}^{2}$ is the room that remains along $g$ beyond
plane $i$: it is positive exactly when that plane cuts the ellipsoid
($|\beta_{i}| \le \tau$) and vanishes when the plane is tangent. The
combination $\eta = \tau^{2} + n\beta_{0}\beta_{1}$ is the admissibility test,
because $\eta > 0$ is equivalent to $\beta_{0}\beta_{1} > -\tau^{2}/n$, the
condition under which a strictly smaller enclosing ellipsoid exists; the case
$\eta \le 0$ is precisely the "no smaller ellipsoid" branch quoted below. The
term $\xi$ is the positive root of the quadratic that remains after the two
tangency conditions are combined; it is the quantity that couples the two
planes, and it collapses to the single-cut value when they coalesce. Finally,
the displacement $\rho = \sigma(\beta_{0}+\beta_{1})/2$ moves the center toward
the middle of the band, which is why a centered band leaves the center fixed.

The update itself is obtained by choosing the least-volume ellipsoid that
contains the intersection of $\mathcal{E}$ with the band. Equivalently, it
minimizes $\ln\det Q^{+}$ over the ellipsoids that contain the band, and the
optimum is characterized by both planes being supporting hyperplanes of the new
ellipsoid. The stationarity conditions force the rank-one form @eq:update with
two scalar multipliers; imposing tangency at both planes eliminates one
multiplier and leaves a quadratic whose positive root is $\xi$, which yields
@eq:parallel_sigma and @eq:parallel_delta. The central and deep cuts are the
limiting tangency cases, and the complete algebra is worked out in
[@frenk1994deep].

The complete update uses the auxiliary quantities
$$\begin{aligned}
\zeta_{0} &= \tau^{2} - \beta_{0}^{2}, \\
\zeta_{1} &= \tau^{2} - \beta_{1}^{2}, \\
\xi &= \sqrt{\zeta_{0}\zeta_{1}
        + \left(\tfrac{n}{2}(\beta_{1}^{2}-\beta_{0}^{2})\right)^{2}},
\end{aligned}$$
{#eq:zeta}
and, with $\eta = \tau^{2} + n\beta_{0}\beta_{1}$,
$$\sigma = \frac{2\eta}{\tau^{2} + \beta_{0}\beta_{1}
            + \tfrac{n}{2}(\beta_{0}+\beta_{1})^{2} + \xi}, \qquad
  \rho = \sigma\cdot\frac{\beta_{0}+\beta_{1}}{2},
$$ {#eq:parallel_sigma}
$$\delta = \frac{n^{2}}{(n^{2}-1)\,\tau^{2}}\left(
          \frac{\zeta_{0}+\zeta_{1}}{2} + \frac{\xi}{n}\right).
$$ {#eq:parallel_delta}

This formulation is finite even when the two offsets are symmetric,
$\beta_{0} + \beta_{1} = 0$, which is precisely the case that defeats a
mean-offset expression dividing by $(\beta_{0}+\beta_{1})^{2}$. It also
degenerates correctly: setting $\beta_{1} = \tau$ makes the second plane
tangent and reduces the parameters to the deep-cut values @eq:deep with
$\beta = \beta_{0}$. The intersection is empty when $\beta_{0} > \beta_{1}$; no
smaller ellipsoid exists when $\beta_{0}\beta_{1} < -\tau^{2}/n$; and the update
reduces to a deep cut when $\beta_{1}^{2} > \tau^{2}$.

```{=latex}
\begin{algorithm}[t]
\caption{Cutting-plane optimization with parallel cuts}
\begin{algorithmic}[1]
\Require oracle $\Omega$, ellipsoid $\mathcal{E} \supseteq \mathcal{K}$, tolerance $\epsilon$
\Ensure value $\gamma^\star$, coefficients $\mathbf{h}_{\mathrm{csd}}$
\State $\gamma \gets +\infty$
\Repeat
    \State $x_c \gets \mathrm{center}(\mathcal{E})$
    \State $(g, \beta_0, \beta_1, t) \gets \Omega(x_c, \gamma)$
    \If{$t < \gamma$}
        \State $\gamma \gets t$;\quad $\mathcal{E} \gets \textsc{CentralCut}(\mathcal{E}, g, 0)$
    \ElsIf{$\beta_1^{2} \le \tau^{2}$}
        \State $\mathcal{E} \gets \textsc{ParallelCut}(\mathcal{E}, g, \beta_0, \beta_1)$
    \Else
        \State $\mathcal{E} \gets \textsc{DeepCut}(\mathcal{E}, g, \beta_0)$
    \EndIf
\Until{$\operatorname{vol}(\mathcal{E}) < \epsilon$}
\State \Return $\gamma^\star, \mathbf{h}_{\mathrm{csd}}$
\end{algorithmic}
\end{algorithm}
```

![Parallel cuts.](ellipsoid.files/parallel_cut.pdf){#fig:parallel_cut}

### Volume Reduction

Applying the matrix-determinant lemma to @eq:update,
$$\begin{aligned}
\det Q^{+} &= (1-\sigma)\det Q, \\
\frac{\operatorname{vol}(\mathcal{E}^{+})}
     {\operatorname{vol}(\mathcal{E})} &= \delta^{n/2}(1-\sigma)^{1/2}
  \le e^{-1/(2n)}.
\end{aligned}$$
{#eq:volume}
After $k$ iterations the volume is at most $e^{-k/(2n)}$ of the initial volume,
so reaching a fraction $\epsilon$ requires $k \approx 2n\ln(1/\epsilon)$
iterations. The determinant of $Q$ therefore decreases monotonically, and the
factor $\delta$ stays in $(0, n^{2}/(n^{2}-1))$.

A robust implementation stores a factorization $Q = \kappa L D L^\mathsf{T}$
with $L$ unit lower triangular and $D$ diagonal, and applies the rank-one
update to $(L,D)$ directly. Positive definiteness is then preserved by
construction, which prevents the silent loss of definiteness that repeated
rank-one downdates can cause in the unfactored form.

## Canonical Signed Digit and Adder Cost {#sec:csd}

A coefficient is realized without a multiplier when it is a signed sum of
powers of two,
$$c = \sum_{j} s_{j}\, 2^{e_{j}}, \qquad s_{j} \in \{-1, 0, +1\},
$$ {#eq:csd_coeff}
because each term is a shift and each pair of terms is an add or subtract. A
canonical signed digit (CSD) representation is one in which no two adjacent
digits are non-zero; it is unique and contains the fewest non-zero digits of
any signed-digit representation of the same number [@george1960csd]. For
example,
$$28.5 = \texttt{+00-00.+0}_{\mathrm{csd}} = 32 - 4 + 0.5,$$
which uses three non-zero digits rather than the four set bits of the binary
form. A power of two, such as $h_{i} = 2^{-3}$, needs no adder at all.

The number of adders per coefficient is one less than the number of non-zero
digits, before any sharing. Because several coefficients are realized
simultaneously, the relevant problem is *multiple constant multiplication*
(MCM): the goal is to build a network of adders and shifts that realizes all
coefficients at minimum cost. *Common subexpression elimination* (CSE) searches
for bit patterns shared across coefficients and reuses their computation,
reducing the adder count and, importantly, the *adder depth* --- the number of
adder stages, which controls the critical path
[@samueli1989improved; @markovic2012dsp]. The number of non-zero digits is
therefore a design parameter, not an accident of the coefficient values.

To expose the trade-off, a coefficient is quantized to the nearest CSD number
with at most a prescribed number of non-zero digits, written
$$h_{\mathrm{csd}}[k] = \operatorname{csd}(h[k], \mathrm{nnz}).
$$ {#eq:csd_quant}
Enlarging the budget $\mathrm{nnz}$ improves the approximation and enlarges the
set of feasible filters, at the price of more adders; shrinking it to one or
two digits makes many otherwise realizable filters infeasible. The budget is a
constraint of the design problem, and it is what makes the problem discrete.

```{=latex}
\begin{algorithm}[t]
\caption{CSD quantization with a non-zero-digit budget}
\begin{algorithmic}[1]
\Require value $x$, budget $\mathrm{nnz}$
\Ensure $c$ with at most $\mathrm{nnz}$ non-zero digits
\State $c \gets 0$;\quad $b \gets 2^{\lceil \log_2(1.5\,|x|)\rceil - 1}$
\While{$\mathrm{nnz} > 0$ \textbf{and} $|x| > 0$}
    \If{$|1.5\,x| > b$}
        \State $c \gets c + \operatorname{sign}(x)\,b$
        \State $x \gets x - \operatorname{sign}(x)\,b$;\quad $\mathrm{nnz} \gets \mathrm{nnz}-1$
    \EndIf
    \State $b \gets b/2$
\EndWhile
\State \Return $c$
\end{algorithmic}
\end{algorithm}
```

### Shift-Add Synthesis and Common Subexpressions

A CSD coefficient is realized by a shift for every non-zero digit and an add or
subtract for every subsequent digit, so a coefficient with $d$ non-zero digits
costs $d-1$ adders before any sharing. Two forms of sharing reduce this cost.
Within a coefficient, a signed pattern that occurs at positions $p$ and $q > p$
satisfies
$$\operatorname{pat}_{q}(x) = \operatorname{pat}_{p}(x) \gg (q-p),$$
so one shift-and-add network serves both occurrences, and the longest repeated
substring of the digit string identifies the pattern worth sharing. Across
coefficients, the same idea selects the pattern that maximizes
$$\text{score} = (\mathrm{nnz}-1)(\text{occurrences}-1),$$
which is exactly the number of adders saved by sharing that pattern. The
percentage saved is a property of the pattern, not of the filter: how much of
the saving a whole design realizes depends on how many patterns are common
across coefficients. The reference generator of @sec:experiments shares one
dominant cross-coefficient pattern, so the whole-design reduction it achieves is
modest, and a search over all patterns would remove more.

### Filter Architecture

The same coefficient set admits a direct or a transposed realization. The direct
form forms each product $h[k]\,x[n-k]$ and accumulates them; the transposed form
pipelines the accumulator and applies the coefficients in reverse order,
producing one output with a single-cycle latency. Both forms use the same CSD
strings and therefore the same per-tap hardware, but they differ in register
count, adder depth, and external wiring. The transposed form is the natural
default for a self-contained module, whereas the direct form exposes the
individual products for a custom accumulator tree or a polyphase decomposition.

## Quantization-Aware Design {#sec:discrete}

Designing a multiplierless filter is not a matter of optimizing first and
rounding afterwards. At the relaxed optimum the active magnitude constraints
are saturated, so the optimum lies on the boundary of the feasible set; the
composite map $S^{-1} \circ \operatorname{csd} \circ S$, which rounds a
coefficient vector and returns its autocorrelation, is neither a projection nor
exact, and even a small perturbation can push part of the passband below its
lower bound. The rounded design can therefore be substantially infeasible.

The remedy is to move the quantizer inside the oracle. The pipeline is
$$\mathbf{r} \;\longrightarrow\; \mathbf{h} = S(\mathbf{r})
  \;\longrightarrow\; \mathbf{h}_{\mathrm{csd}}
  \;\longrightarrow\; \mathbf{r}_{\mathrm{csd}} = S^{-1}(\mathbf{h}_{\mathrm{csd}}),
$$ {#eq:pipeline}
where $S$ is spectral factorization and $S^{-1}$ recomputes the autocorrelation
of the quantized response. The oracle is queried at a continuous center
$\mathbf{r}$ and proceeds as follows.

- **Quantize.** Compute the continuous factor $\mathbf{h} = S(\mathbf{r})$ and
  replace every coefficient by its nearest CSD value under the budget:
  $h_{\mathrm{csd}}[k] = \operatorname{csd}(h[k], \mathrm{nnz})$.
- **Re-evaluate exactly.** Recompute $\mathbf{r}_{\mathrm{csd}} =
  S^{-1}(\mathbf{h}_{\mathrm{csd}})$ and test the sampled constraints against
  $\mathbf{r}_{\mathrm{csd}}$, not against $\mathbf{r}$. A feasible point then
  corresponds to a genuinely realizable coefficient vector. Because the cut is
  generated at the quantized point rather than at the queried center, the oracle
  re-anchors it with the first-order correction $\beta \leftarrow \beta +
  g^\mathsf{T}(\mathbf{r}_{\mathrm{csd}} - \mathbf{r})$, so that the cut is valid
  at $\mathbf{r}$ and can be applied directly.
- **Classify.** The candidate is feasible (improve the best-so-far value and
  take a central cut), violated (take a deep cut from the gradient at
  $\mathbf{r}_{\mathrm{csd}}$), or ineffective (the cut does not shrink the
  ellipsoid).
- **Retry.** The ineffective case is peculiar to the discrete setting: because
  the admissible set is finite, several consecutive centers can map to the same
  CSD pattern, producing a cut that removes no new volume. The oracle replies
  with a retry flag and a budget equal to the number of sampled frequency rows,
  $m = c_{\mathrm{disc}} n$; each retry advances a round-robin cursor over the
  constraints, so that a different violated row is examined instead of the one
  that produced the ineffective cut. When the budget is exhausted the oracle
  returns the last cut even if it is shallow, which lets the ellipsoid update
  proceed; termination is then governed by the usual volume test. Bounding the
  retries by the grid size guarantees that each oracle call makes progress or
  reports its inability to do so, so the outer loop cannot stall on a single
  discrete pattern.

```{=latex}
\begin{algorithm}[t]
\caption{Quantization-aware oracle $\Omega_Q$}
\begin{algorithmic}[1]
\Require center $\mathbf{r}$, budget $\mathrm{nnz}$, retry flag, count $t$, grid size $m$
\Ensure a cut $(g,\beta)$, possibly an improved $\gamma$, and whether a retry is allowed
\If{not retry \textbf{and} $\mathbf{r} \notin \mathcal{K}$}
    \State \Return cut at $\mathbf{r}$
\EndIf
\State $\mathbf{h} \gets \textsc{SpectralFact}(\mathbf{r})$
\State $\mathbf{h}_{\mathrm{csd}} \gets \textsc{CsdQuantize}(\mathbf{h}, \mathrm{nnz})$
\State $\mathbf{r}_{\mathrm{csd}} \gets \textsc{InverseSpectralFact}(\mathbf{h}_{\mathrm{csd}})$
\If{$\mathbf{r}_{\mathrm{csd}} \in \mathcal{K}$}
    \State $\gamma \gets f_0(\mathbf{r}_{\mathrm{csd}})$
    \State \Return $(\partial f_0(\mathbf{r}_{\mathrm{csd}}), 0)$ \Comment{central cut}
\Else
    \State $g \gets \partial f_j(\mathbf{r}_{\mathrm{csd}})$;\quad $\beta \gets f_j(\mathbf{r}_{\mathrm{csd}}) + g^\mathsf{T}(\mathbf{r}_{\mathrm{csd}} - \mathbf{r})$
    \State \Return $(g,\beta)$, retry allowed iff $t < m$ \Comment{deep cut, or retry}
\EndIf
\end{algorithmic}
\end{algorithm}
```

One pitfall deserves emphasis. The oracle certifies feasibility only on the
finite grid of $m = c_{\mathrm{disc}} n$ sampled frequencies. Between samples
the squared magnitude can dip below its lower bound even though the design
passes at every sample; this is a *discretization artifact*, the familiar
consequence of replacing a semi-infinite constraint by finitely many samples. It
is mitigated by enlarging the sampling factor $c_{\mathrm{disc}}$ when the
constraints are tight, and by re-checking the final design on a finer grid. The
fine-grid check is advisory, because the optimizer never certified that
resolution.

## Spectral Factorization {#sec:spectral}

Both the convex formulation and the oracle require a spectral factorization:
given the autocorrelation @eq:autocorr, find the minimum-phase impulse response
$h$ that produces it. The minimum-phase factor is unique when $r$ corresponds
to a positive spectrum, and it is the factor a causal implementation stores.

*Transform-based method (Kolmogorov).* Evaluate the power spectrum
$R(\omega)$ on an oversampled grid, form the log-magnitude
$\alpha(\omega) = \tfrac{1}{2}\ln R(\omega)$, obtain the minimum-phase phase as
its Hilbert transform $\phi = \mathcal{H}[\alpha]$, and return
$$h = \mathcal{F}^{-1}\!\left[\exp\!\left(\alpha(\omega) + j\,\phi(\omega)
  \right)\right].
$$ {#eq:kolmogorov}
The method is non-iterative, deterministic, and stable, and it has no tuning
parameter; its costs are an FFT dependency and a memory footprint proportional
to the oversampling factor, which is customarily $100$ times the filter order.

*Root-based method (Aberth--Ehrlich).* The autocorrelation defines the palindromic
polynomial
$$P(z) = z^{N-1}\left(r[0] + \sum_{k=1}^{N-1} r[k]\,(z^{k} + z^{-k})\right)$$
of degree $2N-2$, whose roots occur in reciprocal pairs: if $\zeta$ is a root then
$1/\bar{\zeta}$ is also a root. All roots are found simultaneously by the
Aberth--Ehrlich iteration
$$\zeta_{i}^{(k+1)} = \zeta_{i}^{(k)}
  - \frac{P(\zeta_{i}^{(k)})}{P'(\zeta_{i}^{(k)})}
  \Big/ \left(1 - \sum_{j \neq i}
  \frac{P(\zeta_{j}^{(k)})}{(\zeta_{i}^{(k)} - \zeta_{j}^{(k)})\,P'(\zeta_{j}^{(k)})}\right),$$
the roots inside the unit circle are retained, and the factor is reconstructed
from them. The method needs no FFT
library, uses memory linear in the order, and exposes a convergence tolerance;
convergence is not guaranteed for pathological inputs, and the smallest
coefficients are recovered with somewhat lower relative accuracy.

Reconstruction from roots is sensitive to the order in which the linear factors
are multiplied. Multiplying well-separated roots first --- a Leja ordering ---
keeps the intermediate coefficients from growing, and reduces the reconstruction
error from roughly $10^{-2}$ to machine precision. A second numerical point is
the convergence test: an absolute residual $|P(z)|$ is not scale-free, and when
the coefficients are large the threshold it implies can lie below the
attainable floating-point floor, so the iteration never reports convergence.
Scaling the polynomial by $s = \max(1, \max_{k} |a_{k}|)$ leaves the iteration
invariant --- both the correction and its denominator scale the same way ---
while making the tolerance relative, which restores convergence without
changing the roots.

The natural measure of a factorization is its round-trip relative error, the
agreement between the prescribed $r$ and the autocorrelation of the returned
factor. The transform method attains close to machine precision, on the order of
$10^{-5}$ in double precision for representative designs, whereas the root-based
method is around $10^{-3}$, concentrated in the smallest coefficients. Because
CSD quantization discards precisely those coefficients once the budget is
applied, the discrepancy is usually absorbed by the quantizer, but it must be
accounted for in verification.

Which method is faster is a property of the implementation rather than of the
algorithm. When the root finder is compiled and the FFT is not, the root-based
method can be the faster choice; when both are compiled, the two methods are
comparable and the crossover depends on the order. This observation motivates
the cross-language experiments of @sec:experiments. The transform method is the
default for production designs, for high-order filters, and for reproducible
benchmarks, because it needs no tuning and is the most robust; the root-based
method is preferable for exploratory work and for moderate orders, where its
speed and tunability dominate.

## Implementation and Experiments {#sec:experiments}

### The Software Pipeline

The pipeline of @eq:pipeline has been implemented three times, in Python, Rust,
and C++, from a single algorithmic specification. The heavy components are
shared: an ellipsoid engine supplies the search space and the generic
cutting-plane loop; a spectral-factorization routine supplies $S$; and a CSD
quantizer supplies $\operatorname{csd}(\cdot,\mathrm{nnz})$. The oracle checks
the passband, stopband, and non-negativity constraints in round-robin order and
returns the first violation it finds, which preserves the round-robin cursor
that prevents the same cut from being reintroduced. The implementations are
checked by property-based testing: a fast implementation is compared against a
reference on randomly generated inputs, which exposes silent conversion and
factorization defects that example-based unit tests miss.

```{=latex}
\begin{algorithm}[t]
\caption{Multiplierless FIR design pipeline}
\begin{algorithmic}[1]
\Require specification $(n, \omega_p, \omega_s, L, U)$, budget $\mathrm{nnz}$
\Ensure a realizable coefficient vector $\mathbf{h}_{\mathrm{csd}}$
\State initialize $\mathcal{E}$ and $\gamma \gets +\infty$
\Repeat
    \State $\mathbf{r} \gets \mathrm{center}(\mathcal{E})$
    \State $(\mathbf{h}_{\mathrm{csd}}, \text{cut}) \gets \Omega_Q(\mathbf{r}, \mathrm{nnz})$
    \State $\mathcal{E} \gets \textsc{Update}(\mathcal{E}, \text{cut})$
\Until{$\operatorname{vol}(\mathcal{E}) < \epsilon$}
\State \Return $\mathbf{h}_{\mathrm{csd}}$
\end{algorithmic}
\end{algorithm}
```

### Experimental Setup

The reference design is a lowpass filter with passband edge $0.12\pi$, stopband
edge $0.20\pi$, a CSD budget of $\mathrm{nnz} = 7$ non-zero digits, a
discretization factor $c_{\mathrm{disc}} = 15$, a tolerance of $10^{-14}$, and an
initial ellipsoid radius of $40$. Four orders are used, $n \in \{16, 32, 64,
128\}$, which range from a small design to one whose dimension makes the
$O(n^{2})$ iteration budget visible. All three implementations read the same
specification and share the same defaults, so any difference in the results is
attributable to the language and its standard libraries rather than to the
algorithm. Measurements were taken on a Windows x64 machine with release builds
of the three implementations; each figure is the mean of five measured runs
taken after two warm-up runs.

### Cross-Language Results

```{=latex}
\begin{table*}[t]
\centering
\caption{Running time and iteration count by filter order, averaged over five
measured runs. The three implementations share one specification; the figures
are indicative wall-clock times from the authors' implementations.}
\label{tbl:results}
\begin{tabular}{lrrrrrr}
\hline
Order $n$ & \multicolumn{3}{c}{Mean time (ms)} & \multicolumn{3}{c}{Iterations} \\
\cline{2-4}\cline{5-7}
          & Python & C++ & Rust & Python & C++ & Rust \\
\hline
$16$  & $1977$ & $61$   & $50$   & $620$  & $625$  & $625$ \\
$32$  & $2786$ & $268$  & $252$  & $1931$ & $1850$ & $2530$ \\
$64$  & $4224$ & $997$  & $1041$ & $1555$ & $1585$ & $1558$ \\
$128$ & $6713$ & $2322$ & $2303$ & $1784$ & $1869$ & $1842$ \\
\hline
\end{tabular}
\end{table*}
```

Two observations stand out. First, the compiled implementations are an order of
magnitude or more faster than the interpreted one at small orders --- about
thirty-two times at $n = 16$ and ten times at $n = 32$ --- but the gap narrows
as the order grows, to about $2.9$ times at $n = 128$. The narrowing is
informative: the interpreted cost is the dispatch of the Python-level loop over
the constraints, which is amortized against an ever-larger amount of FFT work as
the order rises, whereas the compiled cost is dominated by the arithmetic
itself. The two compiled implementations are within ten percent of each other at
every order, and neither dominates consistently: C++ is faster at $n = 32$ and
$64$, Rust at $n = 16$ and $128$.

Second, the iteration counts differ across languages even though the inputs are
identical, and they are not monotone in the order. The ellipsoid method is a
first-order method with no mechanism for converging to a unique optimum, so
floating-point differences in the FFT libraries and in the accumulation order
lead to different, equally valid trajectories; the solutions satisfy the same
specifications and have comparable energy $\sum h[k]^{2}$. The counts also show
that the number of iterations is not a proxy for the order: $n = 32$ requires
more iterations than $n = 64$ in every language, because the parallel-cut
geometry, not the dimension, determines how much volume each cut removes.

### Iteration Economy of Parallel Cuts

To isolate the effect of the parallel-cut update, the same designs were run with
the parallel-cut path disabled, so that the oracle could remove only one
half-space per iteration.

```{=latex}
\begin{table*}[t]
\centering
\caption{Iteration count with a single cut and with parallel cuts. The
single-cut variant did not reach the tolerance within the $50{,}000$-iteration
budget at the two largest orders.}
\label{tbl:parallel_econ}
\begin{tabular}{lrrr}
\hline
Order $n$ & Single cut & Parallel cut & Reduction \\
\hline
$16$  & $2666$            & $625$  & $76.6\%$ \\
$32$  & $19{,}038$        & $1850$ & $90.3\%$ \\
$64$  & $\ge 50{,}000$    & $1585$ & $\ge 96.8\%$ \\
$128$ & $\ge 50{,}000$    & $1869$ & $\ge 96.3\%$ \\
\hline
\end{tabular}
\end{table*}
```

The parallel cut reduces the iteration count by between $77$ and $97$ percent,
and the gap widens with the order: at $n = 64$ and $n = 128$ the single-cut
variant exhausts the $50{,}000$-iteration budget without reaching the tolerance,
whereas the parallel-cut variant converges in fewer than $1{,}900$ iterations.
This is the quantitative form of the observation of @sec:method: the two-sided
magnitude constraints are what the parallel update is for, and a single-cut
method pays for the lower bound and the upper bound in separate iterations.

### Post-Synthesis Hardware Cost

The emitted shift-and-add RTL was synthesized to confirm that the design is
realizable in hardware and to give a technology-independent cost. Two flows were
run in Yosys: a generic gate-level flow, which maps the datapath to elementary
logic and flip-flops and reports a target-independent cell count, and the Xilinx
flow, which maps the same RTL to six-input lookup tables and flip-flops.

```{=latex}
\begin{table*}[t]
\centering
\caption{Post-synthesis cost of the emitted RTL for the transposed
cross-CSE design. The lookup-table and flip-flop counts come from a
six-input-LUT mapping; the generic cell count is target-independent. The design
uses no multiplier or DSP cells.}
\label{tbl:hardware}
\begin{tabular}{lrrr}
\hline
Order $n$ & Lookup tables & Flip-flops & Generic cells \\
\hline
$16$ & $4850$     & $704$  & $27{,}115$ \\
$32$ & $10{,}661$ & $1472$ & $57{,}074$ \\
$64$ & $22{,}946$ & $3264$ & $130{,}225$ \\
\hline
\end{tabular}
\end{table*}
```

The design contains no multiplier or DSP cells, and the cost is dominated by the
carry chains of the adders and by the pipeline registers: the flip-flop count is
essentially $n$ times the output word length. Both the lookup-table and the
generic cell counts grow a little faster than linearly in the order, because the
word length itself grows with $n$. The $n = 128$ design was not carried through
the synthesis flow, which did not complete within a practical time budget on the
available machine; the trend across the three smaller orders is consistent.

These figures should be read as a synthesis-level estimate, not as a power or
area measurement. No place-and-route was performed, no dynamic power was
measured, and the counts are specific to the mapping used; the software running
times of the three implementations are a proxy for algorithmic cost and carry no
direct hardware meaning. A power comparison against a multiplier-based
realization would require a full synthesis, place-and-route, and simulation
flow, which is outside the scope of this article.

### Profiling

Profiling the interpreted implementation produced a counter-intuitive result:
the spectral factorization accounted for only about six percent of the running
time, while the constraint scan accounted for about eighty-three percent. The
scan issued one short dot product per constraint --- on the order of half a
million tiny calls --- and the interpreter dispatch, not the arithmetic,
dominated. Replacing the per-constraint dot products with a single
matrix--vector product followed by boolean masks, while preserving the
round-robin cursor, made the scan primitive about $94$ times faster and the
end-to-end run about $1.6$ times faster. After the change the run was dominated
by library import time.

The same change did not transfer to C++. Porting the vectorization literally
made the compiled implementation about $1.9$ times *slower*, because the
original compile-time per-row loop stops at the first violation and therefore
performs about half the work of a full matrix--vector product. The lesson is
that a vectorization win observed in an interpreted language must be
re-measured after compilation; the measurement, not the conclusion, is what
ports. Exposing raw pointers to the compiler and using a short chain of
independent partial sums (rather than a single scalar reduction, which the
compiler may not reassociate) recovered a further small margin in the compiled
implementations.

### Optimization and Sharing Results

The same kernel changes paid off very differently across the three languages:

```{=latex}
\begin{table*}[t]
\centering
\caption{Effect of the kernel changes in three languages.}
\begin{tabular}{llrrr}
\hline
Change & Language & Before & After & Speedup \\
\hline
Vectorized constraint scan & Python & 3.05\,s & 0.082\,s & $37\times$ \\
End-to-end after the scan & Python & 3.01\,s & 1.87\,s & $1.6\times$ \\
Raw-pointer dot product & C++ & 0.172\,s & 0.160\,s & $1.07\times$ \\
Four-accumulator dot & Rust & 23.3\,ns & 10.1\,ns & $2.3\times$ \\
Oracle scan & Rust & $13.9\,\mu\mathrm{s}$ & $10.7\,\mu\mathrm{s}$ & $1.3\times$ \\
\hline
\end{tabular}
\end{table*}
```

The CSD coefficient set is compressed further by common-subexpression sharing:

```{=latex}
\begin{table*}[t]
\centering
\caption{Adder count with and without common-subexpression sharing.}
\begin{tabular}{lrrr}
\hline
CSD pattern & Flat adders & Shared adders & Saving \\
\hline
\texttt{+0-0+0-0} & 4 & 2 & 50\% \\
\texttt{+0-0+0-0+0-0} & 6 & 2 & 67\% \\
\texttt{+00-00+00-00} & 4 & 2 & 50\% \\
\hline
\end{tabular}
\end{table*}
```

The per-pattern figures above are the savings from a single shared pattern, not
from a whole design. For $\mathrm{nnz} = 7$ the flat count is six adders per tap
--- $96$ adders at $n = 16$ and $768$ at $n = 128$ --- and after the cross-pattern
sharing performed by the reference generator the counts are $91$ and $745$,
respectively. Closing the remaining gap requires a
multiple-constant-multiplication search over all patterns, which is left to
future work (@sec:conclusion).

### Multiplierless Result

![Magnitude response of a lowpass design.](ellipsoid.files/lowpass.pdf){#fig:lowpass}

![Magnitude response after CSD quantization.](ellipsoid.files/csdlowpass.pdf){#fig:csdlowpass}

The output of the pipeline is the quantized coefficient vector
$\mathbf{h}_{\mathrm{csd}}$ together with its adder network. Because the
quantizer was consulted during the optimization rather than after it, the
returned design satisfies the sampled specifications by construction, and the
only remaining verification is to confirm that the design also holds between
the samples of the grid. The shift-and-add description is then emitted directly
for synthesis, and the post-synthesis cell and register counts of the preceding
subsection confirm that the result is realizable in hardware without a
multiplier.

## Concluding Remarks {#sec:conclusion}

The ellipsoid method is an unusually good fit for multiplierless FIR design.
Its separation oracle accommodates the semi-infinite magnitude constraints
without ever enumerating them; the two-sided nature of those constraints is
returned as parallel cuts, which reduce the iteration count; and the discrete
CSD quantizer can be embedded in the oracle so that the optimizer searches over
realizable coefficient vectors rather than over a relaxation. The same
formulation applies unchanged to the robust and parametric variants of the
design problem, and the numerical experiments confirm that a compiled
implementation of the whole pipeline runs in under $2.4$ seconds for orders up to
$128$, that the parallel-cut update keeps the iteration count below $1{,}900$
where a single-cut method does not converge, and that the emitted RTL synthesizes
to logic and registers without any multiplier cells.

Three directions remain open. The first is a sharper analysis of the parallel
cut under a discretization grid, so that the fine-grid advisory check can be
replaced by a guarantee. The second is to fold the adder-sharing stage --- the
multiple-constant-multiplication network that realizes the CSD coefficients ---
into the oracle, so that adder cost, and not only coefficient error, is a
first-class design objective. The third is a systematic study of the filter
architecture: the direct and transposed forms realize the same transfer
function but differ substantially in adder depth and register count, and the
choice interacts with the coefficient ordering that the CSD quantizer produces;
the synthesis counts reported in @sec:experiments are a first datapoint, not a
substitute for that study.

## Acknowledgments

The author is grateful to the contributors of the open-source packages on which
this work depends.

## References {-}

\
