---
title: Ellipsoid Method and the Amazing Oracles
subtitle: Separation oracles, parallel cuts, and the practice of convex optimization
author: Wai-Shing Luk
institute: Fudan University
date: \today
---

# Agenda

## 📋 Agenda

:::: {.columns}
::: {.column width="47%"}

**🧱 Part 1 — Foundations**
why the ellipsoid method? · separation oracles · cut types · termination

**🧠 Part 2 — Amazing Oracles**
robust optimization · parametric networks · matrix inequalities

:::
::: {.column width="47%"}

**🥚 Part 3 — Ellipsoid Revisited**
search space · deep & central cuts · parallel cuts · stable factorization

**🔬 Part 4 — Applications & Practice**
FIR design · discrete optimization · implementation · benchmarks

:::
::::

# Foundations

## 🤔 Why the Ellipsoid Method?

- The interior-point method evaluates **all** constraints; the ellipsoid method needs only a **separation oracle** 🔎
- Naturally suited to a **moderate number of variables** with a **huge or infinite** number of constraints ♾️
- It cannot exploit **sparsity** — but the oracle can exploit **structure** 🧩
- One algorithm serves feasibility, optimization, discrete design, and bisection 🎯

> The oracle, not the ellipsoid, does the intellectual work. 🧠

## 🔎 Separation Oracle: the Real Star

When queried at $x_0 \in \mathbb{R}^n$, the oracle $\Omega$ answers one of two ways:

1. “$x_0 \in \mathcal{K}$” ✅ — stop, the point is feasible, or
2. a hyperplane separating $x_0$ from $\mathcal{K}$:

$$g^\mathsf{T}(x - x_0) + \beta \le 0, \quad \beta \ge 0, \; g \neq 0, \; \forall x \in \mathcal{K}.$$

```{=latex}
\begin{center}
\begin{tikzpicture}[scale=1.0]
  \draw[fill=green!12, draw=green!45!black, line width=1pt]
    plot[smooth cycle, tension=0.8]
    coordinates {(-2,0.6) (-1.4,1.6) (-0.2,1.8) (0.6,1.0) (0.4,-0.2) (-0.8,-0.6) (-1.9,-0.2)};
  \node at (-0.7,0.8) {$\mathcal{K}$};
  \draw[nordblue, line width=1pt] (-0.5,2.6) -- (2.9,0.4);
  \fill[nordred] (2.2,2.4) circle (1.8pt);
  \node[nordred, above right] at (2.2,2.4) {$x_0$};
  \draw[ar] (1.15,1.55) -- (1.85,2.1) node[above, font=\tiny] {$g$};
\end{tikzpicture}
\end{center}
```

## ✂️ Three Flavours of Cutting Planes

| cut | condition | geometry |
|:--|:--|:--|
| **central** ⭕ | $\beta = 0$ | the plane passes through $x_0$ |
| **deep** 🔻 | $\beta > 0$ | $x_0$ is strictly cut off |
| **shadow** 🌓 | $\beta < 0$ | $x_0$ already lies inside |

- A convex $f$ gives a cut for free: $(g, \beta) = (\partial f(x_0),\, f(x_0))$ 🪄
- Differentiable $f$ → the subgradient is the gradient.

## 🎯 From Feasibility to Optimization

- Add a level constraint $f_0(x) \le \gamma$ and shrink the **sublevel set** as soon as a better point is found:

$$
\begin{array}{ll}
    \text{minimize} & \gamma, \\
    \text{subject to} & \Phi(x, \gamma) \le 0, \; x \in \mathcal{K}.
  \end{array}
$$

- A **binary search** on $\gamma$ is just a nested feasibility problem 🔁
- `no solution` ⇒ infeasible; `no effect` ⇒ the cut was too shallow.

## 🔁 The Cutting-plane Loop

```{=latex}
\begin{center}
\begin{tikzpicture}[node distance=8mm and 11mm]
  \node[nblue] (init) {initial ellipsoid\\$\mathcal{E}_0 \supseteq \mathcal{K}$};
  \node[nyellow, right=of init] (q) {query oracle\\at $x_c$};
  \node[nred, right=of q] (chk) {feasible?};
  \node[ngreen, right=of chk] (done) {return $x_c$ ✅};
  \node[nblue, below=11mm of q] (upd) {update ellipsoid\\$\mathcal{E} \leftarrow$ smaller};
  \node[nred, left=of upd] (small) {volume small\\or empty?};
  \draw[ar] (init) -- (q);
  \draw[ar] (q) -- (chk);
  \draw[ar] (chk) -- node[above, font=\tiny] {yes} (done);
  \draw[ar] (chk) -- node[right, font=\tiny] {no, cut $(g,\beta)$} (upd);
  \draw[ar] (upd) -- (small);
  \draw[ar] (small) -- node[left, font=\tiny] {no} (q);
\end{tikzpicture}
\end{center}
```

## ⏱️ Termination and the Stall Guard

- Comparing a **scale-dependent width** against an **absolute** tolerance can never succeed 🛑
- Once the bracket is at machine resolution, the midpoint is the only representable point
- The counter pins and inner subproblems are solved for nothing ⏳
- **Remedy:** terminate when the midpoint is no longer strictly inside the bracket — a *stall guard*, not a tighter tolerance 🛡️

> An absolute stopping threshold should always be paired with a stall test. 📌

# Amazing Oracles

## 🧠 Three Oracles, One Method

```{=latex}
\begin{center}
\begin{tikzpicture}[node distance=8mm and 14mm]
  \node[ngreen] (m) {cutting-plane\\method};
  \node[nblue, above right=10mm and 20mm of m] (r) {robust\\optimization};
  \node[nblue, right=20mm of m] (n) {parametric\\network};
  \node[nblue, below right=10mm and 20mm of m] (l) {matrix\\inequalities};
  \node[nyellow, right=of r, font=\tiny] (ra) {affine\\arithmetic};
  \node[nyellow, right=of n, font=\tiny] (nn) {negative\\cycle};
  \node[nyellow, right=of l, font=\tiny] (lc) {Cholesky};
  \draw[ar] (m) -- (r);
  \draw[ar] (m) -- (n);
  \draw[ar] (m) -- (l);
  \draw[ar] (r) -- (ra);
  \draw[ar] (n) -- (nn);
  \draw[ar] (l) -- (lc);
\end{tikzpicture}
\end{center}
```

## 🛡️ Robust Convex Optimization

- The worst case over an uncertainty set $\mathcal{Q}$ — a **robust counterpart** stays convex, but its constraints become infinite ♾️

$$
\begin{array}{ll}
    \text{minimize} & \sup_{q \in \mathcal{Q}} f_0(x, q), \\
    \text{subject to} & f_j(x, q) \le 0, \; \forall q \in \mathcal{Q}, \; j = 1, \ldots, m.
  \end{array}
$$

- Oracle: maximize over $q$ and return the cut at the maximizer 🎯
- **Affine arithmetic** propagates uncertainty symbolically, avoiding the $2^K$ vertices of a box 🧮
- Certified feasible iff the affine upper bound $\bar{f}_j(x) \le 0$ for every $j$.

## 🕸️ Parametric Network Problems

- A directed graph $G = (V, E)$ with $u_i - u_j \le h_{ij}(x, \gamma)$ 🕸️
- Feasible **iff** $G$ has **no negative cycle** — the oracle is negative-cycle detection
- Equivalently: $W_k(x, \gamma) = \sum_{(i,j) \in C_k} h_{ij}(x, \gamma) \ge 0$ for every cycle $C_k$
- Fast detectors: Bellman–Ford, Tarjan (1976), Howard (policy iteration) ⚡
- Application: optimal **matrix scaling** under the min–max-ratio criterion ⚖️

## 🔷 Matrix Inequalities and SDP

- Seek $x$ with $F(x) \succeq 0$, where $F(x) = F_0 + x_1 F_1 + \cdots + x_n F_n$
- $A \succeq 0 \iff v^\mathsf{T} A v \ge 0$ for all $v$
- **Cholesky / $LDL^\mathsf{T}$ witness:** if the factorization fails at row $p$, then

$$v = R_{:p,:p}^{-1} e_p, \qquad v^\mathsf{T} F_{:p,:p}(x_0)\, v < 0,$$

$$\text{cut} \;=\; \bigl(-v^\mathsf{T} \partial F_{:p,:p}\, v, \; -v^\mathsf{T} F_{:p,:p}\, v\bigr).$$

- Examples: **minimum eigenvalue** (EVP) and matrix-norm minimization 📐

## 🔬 Estimating a Correlation Function

- The Gaussian criterion combines a **convex** trace term with a **concave** log-determinant — a difference of convex functions 🎢
- It is convex exactly on the trust region $0 \prec \Omega \preceq 2Y$
- Beyond it, apply the **convex–concave procedure (CCP)**: majorize the concave part by its tangent, minimize the convex surrogate, repeat 🔁
- Do **not** impose $0 \preceq \Omega \preceq 2Y$ as a hard constraint — it biases the estimate ⚠️

# Ellipsoid Revisited

## 🥚 The Ellipsoid as a Search Space

- The search space must be cheap to store, cheap to update, and guaranteed to contain $\mathcal{K}$ — the ellipsoid is all three 🥚

$$\mathcal{E} = \{\, x \mid (x - x_c)^\mathsf{T} Q^{-1} (x - x_c) \le \kappa \,\}$$

- $O(n^2)$ parameters · closed-form update · containment by induction 🎁
- Eigenvectors give the axes, eigenvalues the lengths; the volume is

$$\operatorname{vol}(\mathcal{E}) = \kappa^{n/2} \sqrt{\det Q}\; \operatorname{vol}(\mathcal{B}^n).$$

## 🔄 Updating the Ellipsoid (Deep Cut)

- Let $\tilde g = Q g$, $\omega = g^\mathsf{T} \tilde g$, and $\tau = \sqrt{\kappa\, \omega}$ 🔧

$$x_c^+ = x_c - \frac{\rho}{\omega}\tilde g, \qquad
  Q^+ = Q - \frac{\sigma}{\omega}\tilde g\tilde g^\mathsf{T}, \qquad
  \kappa^+ = \delta\,\kappa.$$

- **Deep cut** ($\beta > 0$): $\rho = \dfrac{\tau + n\beta}{n+1}$, $\sigma = \dfrac{2\rho}{\tau+\beta}$, $\delta = \dfrac{n^2}{n^2-1}\cdot\dfrac{\tau^2 - \beta^2}{\tau^2}$
- Empty if $\beta > \tau$; too shallow if $n\beta < -\tau$ ⚠️

## ⭕ Central Cut

- The special case $\beta = 0$ is much simpler — worth its own implementation ⭕

$$\rho = \frac{\tau}{n+1}, \qquad
  \sigma = \frac{2}{n+1}, \qquad
  \delta = \frac{n^2}{n^2-1}.$$

- Same update map, cheaper coefficients 🎯

## 🪜 Parallel Cuts

- The oracle returns a **pair** of cuts sharing a normal $g$:

$$g^\mathsf{T}(x - x_c) + \beta_1 \le 0, \qquad g^\mathsf{T}(x - x_c) + \beta_2 \ge 0, \qquad \forall x \in \mathcal{K}.$$

- Produced by any two-sided linear constraint $l \le a^\mathsf{T} x + b \le u$ 📏
- The update removes a **slab**, not a half-space → faster convergence 🚀

![Parallel cuts](ellipsoid.files/parallel_cut.pdf){height="3.0cm"}

## 🪜 Parallel-cut Parameters

With $\tilde g = Q g$, $\tau^2 = \kappa\omega$, and the auxiliary quantities

$$\zeta_0 = \tau^2 - \beta_0^2, \quad
  \zeta_1 = \tau^2 - \beta_1^2, \quad
  \xi = \sqrt{\zeta_0\zeta_1 + \Bigl(\tfrac{n}{2}(\beta_1^2 - \beta_0^2)\Bigr)^2},$$

the scalar triple is

$$\sigma = \frac{2\eta}{\tau^2 + \beta_0\beta_1 + \tfrac{n}{2}(\beta_0+\beta_1)^2 + \xi}, \qquad
  \rho = \sigma\cdot\frac{\beta_0+\beta_1}{2},$$

$$\delta = \frac{n^2}{(n^2-1)\,\tau^2}\left(\frac{\zeta_0+\zeta_1}{2} + \frac{\xi}{n}\right), \qquad \eta = \tau^2 + n\beta_0\beta_1.$$

> This form stays finite when $\beta_0 + \beta_1 = 0$ — the symmetric case a mean-offset formula cannot handle. 🧩

## 📉 How Fast Does the Volume Shrink?

- Apply the matrix-determinant lemma to the rank-one downdate:

$$\det Q^+ = (1 - \sigma)\det Q
  \quad\Longrightarrow\quad
  \frac{\operatorname{vol}(\mathcal{E}^+)}{\operatorname{vol}(\mathcal{E})} \le e^{-1/(2n)}.$$

- After $k$ iterations the volume is at most $e^{-k/(2n)}$ of the initial volume 📉
- To reach a fraction $\epsilon$ of it: $k \approx 2n \ln(1/\epsilon)$ iterations 🎯

## 🛡️ Numerical Drift and the Stable $LDL^\mathsf{T}$ Variant

- Repeated rank-one downdates can lose positive definiteness — **silently** 😱
- Keep a factored shape instead:

$$Q = \kappa\, L\, D\, L^\mathsf{T}, \qquad
  w = L^{-1} g, \quad z = D^{-1} w, \quad \omega = w^\mathsf{T} z, \quad q = L^{-\mathsf{T}} z.$$

- Three triangular sweeps + a factor update; positive definiteness **by construction** 🛡️
- Compiled cost is comparable to the direct form — near-free insurance ✅

# Applications

## 🎛️ FIR Filter Design

- Magnitude constraints $L(\omega) \le |H(\omega)| \le U(\omega)$ are **not** convex in $\mathbf{h}$ 🚫
- **Spectral factorization** rewrites them convexly in the autocorrelation $\mathbf{r}$ 🪄

$$L^2(\omega) \le R(\omega) \le U^2(\omega), \qquad R(\omega) = |H(\omega)|^2.$$

- Two-sided bounds ⇒ **parallel cuts** ⇒ markedly fewer iterations 🚀

![Lowpass result](ellipsoid.files/lowpass.pdf){height="2.7cm"}

## 🏭 The Multiplierless FIR Pipeline

```{=latex}
\begin{center}
\begin{tikzpicture}[node distance=7mm and 5mm]
  \node[nblue, font=\tiny] (s) {filter\\spec};
  \node[ngreen, font=\tiny, right=of s] (c) {convex\\magnitude\\design};
  \node[nyellow, font=\tiny, right=of c] (e) {parallel-cut\\ellipsoid};
  \node[npurple, font=\tiny, right=of e] (f) {spectral\\factorization};
  \node[nred, font=\tiny, right=of f] (q) {CSD\\quantization};
  \node[ngreen, font=\tiny, right=of q] (v) {synthesizable\\Verilog};
  \draw[ar] (s) -- (c);
  \draw[ar] (c) -- (e);
  \draw[ar] (e) -- (f);
  \draw[ar] (f) -- (q);
  \draw[ar] (q) -- (v);
\end{tikzpicture}
\end{center}
```

- A CSD coefficient is a signed sum of powers of two → shifts and adds, **no multipliers** 🔢
- The non-convex quantization constraint is folded **into the oracle** 🎯

## 🌊 Spectral Factorization: FFT vs Roots

| | FFT (Kolmogorov) | root-finding (Aberth–Ehrlich) |
|:--|:--|:--|
| transform needed | FFT library | none |
| memory | $O(\text{oversample}\cdot N)$ | $O(N)$ |
| tuning | none | tolerance |
| round-trip error | $\sim 10^{-5}$ | $\sim 10^{-3}$ |
| use for | production, high order | exploration, moderate order |

- Both return a valid minimum-phase factor — they differ in the numerical profile, not the mathematics 🔬
- The optimizer may run fewer iterations with the better-conditioned root-based factor 🔁

## 🔢 Discrete Optimization

- Some variables live in a discrete set $\mathbb{D}$ — mixed-integer convex programming 🎲
- Relaxation + branch-and-bound **ignores convexity** ❌
- The oracle looks for a nearby discrete point $x_d$ and cuts there:

$$g^\mathsf{T}(x - x_d) + \beta \le 0, \quad \beta \ge 0, \; g \neq 0.$$

- The cut can be **shallow**; on “no effect”, perturb and **retry** 🔁

| | quantize-after | discrete in oracle |
|:--|:--|:--|
| constraint | after optimization | inside the oracle |
| domain | $\mathbb{R}^N$ | CSD-realizable |
| guarantee | none | at sampled frequencies |
| cost | single pass | multiple passes with retries |

## 🔬 Maximum Likelihood and CCP

- The criterion $\log\det V + \operatorname{Tr}(V^{-1}Y)$ is convex on $0 \preceq V \preceq 2Y$ ✅
- Beyond the region, **CCP** turns one non-convex problem into a sequence of convex ones 🔁
- A triangular factor must be inverted by back-substitution — a symmetric routine silently corrupts the gradient ⚠️

# Implementation & Practice

## 🏗️ Architecture: Algorithm · SearchSpace · Oracle

```{=latex}
\begin{center}
\begin{tikzpicture}[node distance=9mm and 16mm]
  \node[nyellow] (alg) {generic\\algorithm};
  \node[nblue, above right=9mm and 18mm of alg] (ss) {search space\\(ellipsoid)};
  \node[ngreen, below right=9mm and 18mm of alg] (or) {oracle\\(the problem)};
  \draw[ar] (alg) -- node[above, font=\tiny] {center $x_c$} (ss);
  \draw[ar] (ss) -- node[right, font=\tiny] {cut $(g,\beta)$} (or);
  \draw[ar] (or) -- node[below, font=\tiny] {feedback} (alg);
\end{tikzpicture}
\end{center}
```

- One problem-agnostic loop serves **feasibility**, **optimization**, **quantized**, and **bisection** oracles 🧩
- Swap the search space (dense · factored · interval) without touching the algorithm 🔧

## ⚡ Performance Lessons

- **Profile first:** in an interpreted oracle, the *number of calls* dominated, not the arithmetic 🔍
- Replace per-constraint callbacks with **one matrix-vector product + masks**, preserving the round-robin cursor 🎯
- **Per-iteration allocation is the tax** — pre-allocate and reuse scratch buffers 🗑️
- Early exit can beat vectorization in compiled code — **port the measurement, not the conclusion** 📌
- An unchanged iteration count on fixed input is strong evidence that a rewrite preserved behaviour ✅

## ⚖️ Ellipsoid vs Interior-Point

| problem | interior-point | ellipsoid (interp.) | ellipsoid (compiled) |
|:--|--:|--:|--:|
| Chebyshev center, $n=10$ | $1.0\times$ | $2.9\times$ | $0.03\times$ |
| min-eigenvalue LMI, $m=8$ | $1.0\times$ | $14.6\times$ | $0.08\times$ |
| FIR lowpass, $n=48$ | $1.0\times$ | $7.7\times$ | $0.26\times$ |

- Runtimes relative to the interior-point solver; **below $1.0$ is faster** ⚖️
- The per-iteration constant (~language overhead) can flip the ranking 🔁
- Use the ellipsoid method when constraints are reachable **only through an oracle** 🔎

# Closing

## 🎯 Key Takeaways

:::: {.columns}
::: {.column width="47%"}

**The idea** 💡

- the **oracle** is the intellectual centre — the ellipsoid is bookkeeping 🧠
- **parallel cuts** pay off whenever bounds are two-sided 🪜
- discrete and robust problems fit the same loop 🎯

:::
::: {.column width="47%"}

**The practice** 🔧

- factor the ellipsoid ($LDL^\mathsf{T}$) for **free insurance** 🛡️
- profile, pre-allocate, and re-measure 📉
- pick the method by **constraint access**, not by reputation 🧭

:::
::::

## 📚 References

- **Boyd & Vandenberghe**, *Convex Optimization* — interior-point methods 📕
- **Bland, Goldfarb & Todd (1981)** — the ellipsoid method: a survey 📐
- **Wu et al. (1999)** — FIR design via spectral decomposition 🎛️
- **Goodman (1997)** — spectral factorization for FIR design 🌊
- **Liu et al. (2007)** — affine arithmetic for robust geometric programming 🛡️

> Full citations and the complete argument live in the paper **“Ellipsoid Method and the Amazing Oracles.”** 📄

## 🙋 Q&A

**Ellipsoid Method and the Amazing Oracles**

Questions? Discussion? 💬

## 👏 Thank You

The oracle does the thinking; the ellipsoid just shrinks. 🥚🔻

Slides built with Beamer · TikZ 🧩 · LuaLaTeX 📐 · Nord 🌙
