---
title: Ellipsoid Method and the Amazing Oracles
bibliography: ['ellipsoid.bib', 'fir-ref.bib', 'Geostatistics.bib', 'mpcss1.bib', 'mpcss2.bib']
abstract: |
  'The ellipsoid method is revisited. In addition to this, three separation oracles are investigated. They are robust optimization, semidefinite programming and network optimization. Stability issues are discussed. Finally, parallel-cut is introduced.'
---

# Introduction

The ellipsoid method has a bad reputation. The method is often considered slow for large-scale convex problems compared to the interior-point method [@unknown]. And this is not fair. First, unlike the interior-point method, the ellipsoid method does not require an explicit evaluation of all constraint functions at each iteration. All it needs is a *separation oracle* that provides a *cutting-plane* (@sec:cutting-plane). This can make the method attractive for problems with a large number or even an infinite number of constraints. Second, while the ellipsoid method itself cannot exploit the sparsity of the problem, the separation oracle can exploit certain types of structure of the problems...

In @sec:robust, robust optimization...

In @sec:network, we show that for the network parametric problems, the cutting-plane can be obtained by finding a negative cycle of a directed graph. Efficient algorithms exist in which the locality of the network as well as other properties can be exploited.

In @sec:lmi, problems involving matrix inequalities are discussed. Recall that the positive definiteness of a symmetric matrix can be checked efficiently using Cholesky factorization, or more precisely $LDL^\mathsf{T}$ factorization. Let a symmetric matrix $A \in \mathbb{R}^{m \times m}$. Recall that if the factorization process stops at row $p$ due to encountering a non-positive diagonal entry of $A$, then $A$ is not positive-definite. Furthermore, by using lazy evaluation techniques, it is possible to construct a cutting plane in $O(p^3)$ instead of $O(m^3)$. Thus, it can be used for efficient oracle implementations.

The implementation of the ellipsoid method is discussed in @sec:ellipsoid. This method is a cutting-plane method in which the search space is an ellipsoid, usually expressed as:
$$\{ x \mid (x-x_c)P^{-1}(x-x_c) \le 1 \},$$
where $x_c \in \mathbb{R}^n$ is the center of the ellipsoid. The matrix $P \in \mathbb{R}^{n \times n}$ is symmetric positive-definite. At each iteration, $x_c$ and $P$ are updated according to the oracle. While the updating of ellipsoids is simple, and methods have been found to implement it for decades, we show that by splitting $P$ into multiplying $\alpha$ by $Q$, i.e.,
$$\{ x \mid (x-x_c)Q^{-1}(x-x_c) \le \alpha \},$$
the updating cost can be further reduced by $n^2$ flops (floating-point operations).

Besides, the use of parallel cuts is discussed in Section 4.2. When a pair of parallel inequalities, one of which is violated, two constraints can be used simultaneously to update the ellipsoid. Some papers reported that this technique did not provide significant improvements. However, we show that for cases where the upper and lower bounds of some constraints are tight, such as in the case of some filter designs, the use of parallel cuts can dramatically improve the runtime (of what?) . Furthermore, we show that if the method is implemented correctly, each update, whether using a deep cut or a parallel cut, at most one square root operation will be required.

In many practical engineering problems, some design variables may be restricted to discrete forms. Since the cutting-plane method requires only a separation oracle, it can also be used for discrete problems...

Cutting-plane Method Revisited {#sec:cutting-plane}
==============================

Convex Feasibility Problem
--------------------------

Let $\mathcal{K} \subseteq \mathbb{R}^n$ be a convex set. Consider the feasibility problem:

1.  find a point $x^* \in \mathbb{R}^n$ in $\mathcal{K}$, or

2.  determine that $\mathcal{K}$ is empty (i.e., has no feasible solution).

When a *separation oracle* $\Omega$ is *queried* at $x_0$, it either

1.  asserts that $x_0 \in \mathcal{K}$, or

2.  returns a separating hyperplane between $x_0$ and $\mathcal{K}$:

$$g^\mathsf{T} (x - x_0) + \beta \le 0, \beta \ge 0, g \neq 0, \; \forall x \in \mathcal{K}.$$
{#eq:cut}

The pair of $(g, \beta)$ is called a *cutting-plane*, because it eliminates the half-space $\{x \mid g^\mathsf{T} (x - x_0) + \beta > 0\}$ from our search. We have the following observations:

-   If $\beta=0$ ($x_0$ is on the boundary of the half-space), the cutting-plane is called *neutral-cut*.

-   If $\beta>0$ ($x_0$ lies in the interior of the half-space), the cutting-plane is called *deep-cut*.

-   If $\beta<0$ ($x_0$ lies in the exterior of the half-space), the cutting-plane is called *shadow-cut*.

$\mathcal{K}$ is usually given by a set of inequalities $f_j(x) \le 0$ or $f_j(x) < 0$ for $j = 1 \cdots m$, where $f_j(x)$ is a convex function.
A vector $g \equiv \partial f(x_0)$ is called the *sub-gradient* of a convex function $f$ at $x_0$ if $f(z) \ge f(x_0) + g^\mathsf{T} (z - x_0)$.
Thus, the cut $(g, \beta)$ is given by $(\partial f(x_0), f(x_0))$.

Note that if $f(x)$ is differentiable, we can simply take $\partial f(x_0) = \nabla f(x_0)$.
The cutting-plane method consists of two key components: separation oracle $\Omega$ and a search space $\mathcal{S}$ initially sufficiently large to cover $\mathcal{K}$.
For example,

-   Polyhedron $\mathcal{P}$ = $\{z \mid C z \preceq d \}$.

-   Interval $\mathcal{I}$ = $[l, u]$ (for one-dimensional problem).

-   Ellipsoid $\mathcal{E}$ = $\{z \mid (z-x_c)P^{-1}(z-x_c) \le 1 \}$.

Generic cutting-plane method:

-   **Given** initial $\mathcal{S}$ known to contain $\mathcal{K}$.

-   **Repeat**

    1.  Select a point $x_0$ in $\mathcal{S}$.

    2.  Query the separation oracle at $x_0$.

    3.  **If** $x_0 \in \mathcal{K}$, exit.

    4.  **Else**, update $\mathcal{S}$ to a smaller set which covers:
        $$\mathcal{S}^+ = \mathcal{S} \cap \{z \mid g^\mathsf{T} (z - x_0) + \beta \le 0\}.$$

    5.  **If** $\mathcal{S}^+ = \emptyset$ or it is small enough, exit.

Todo: What if the search space is not large enough?

From Feasibility to Optimization
--------------------------------

Consider:
$$\begin{array}{ll}
    \text{minimize}     & f_0(x), \\
    \text{subject to}   & x \in \mathcal{K}.
  \end{array}$$
{#eq:convex-optimization}
We treat the optimization problem as a feasibility problem with an additional constraint $f_0(x) \le t$.
Here, $f_0(x)$ can be a convex function or a quasi-convex function.
$t$ is the best-so-far value of $f_0(x)$.
We can reformulate the problem as:
$$\begin{array}{ll}
    \text{minimize}   & t, \\
    \text{subject to} & \Phi(x, t) \le 0, \\
                      & x \in \mathcal{K},
  \end{array} $$
{#eq:cvx-in-feasibility-form}
where $\Phi(x, t) \le 0$ is the $t$-sublevel set of $f_0(x)$.

For every $x$, $\Phi(x, t)$ is a non-increasing function of $t$, *i.e.*, $\Phi(x, t’) \le \Phi(x, t)$ whenever $t’ \ge t$. Note that $\mathcal{K}_t \subseteq \mathcal{K}_u$ if and only if $t \le u$ (monotonicity). An easy way to solve the optimization problem is to apply a binary search on $t$.

Another possibility is to update the best-so-far $t$ whenever a feasible solution $x_0$ is found such that $\Phi(x_0, t) = 0$. We assume that the oracle takes responsibility for that.

Generic Cutting-plane method (Optim)

-   **Given** initial $\mathcal{S}$ known to contain $\mathcal{K}_t$.

-   **Repeat**

    1.  Choose a point $x_0$ in $\mathcal{S}$

    2.  Query the separation oracle at $x_0$

    3.  **If** $x_0 \in \mathcal{K}_t$, update $t$ such that $\Phi(x_0, t) = 0$.

    4.  Update $\mathcal{S}$ to a smaller set that covers:
        $$\mathcal{S}^+ = \mathcal{S} \cap \{z \mid g^\mathsf{T} (z - x_0) + \beta \le 0\} $$

    5.  **If** $\mathcal{S}^+ = \emptyset$ or it is small enough, exit.

Example: Profit Maximization {#sec:profit}
------------------------------------------

This example is taken from [@Aliabadi2013Robust]. Consider the following *short-run* profit maximization problem:
$$\begin{array}{ll}
   \text{maximize}   & p(A x_1^\alpha x_2^\beta) - v_1 x_1 - v_2 x_2, \\
   \text{subject to} & x_1 \le k, \\
                     & x_1 > 0, x_2 > 0,
  \end{array}$$
{#eq:profit-max-in-original-form}
where $p$ is the market price per unit, $A$ is the scale of production, $\alpha$ and $\beta$ are output elasticities, $x_i$ and $v_i$ are the i-th input quantity and output price, $A x_1^\alpha x_2^\beta$ is the Cobb-Douglas production function, and $k$ is a constant that limits the quantity of $x_1$. The above formulation is not in convex form. First, we rewrite the problem:
$$\begin{array}{ll}
    \text{maximize}   & t, \\
    \text{subject to} & t  + v_1 x_1  + v_2 x_2 \le p A x_1^{\alpha} x_2^{\beta}, \\
                      & x_1 \le k, \\
                      & x_1 > 0, x_2 > 0.
  \end{array} $$
By the change of variables, we can obtain the following convex form of\ @eq-profit-max-in-orginal-form:
$$\begin{array}{ll}
    \text{maximize}   & t, \\
    \text{subject to} & \log(t + v_1 e^{y_1} + v_2 e^{y_2}) -
                    (\alpha y_1 + \beta y_2) \le \log(p\,A), \\
                      & y_1 \le \log k,
  \end{array}$$
{#eq:profit-in-cvx-form}
where $y_1 = \log x_1$ and $y_2 = \log x_2$. Some readers may recognize that we can also write the problem in a geometric program by introducing one additional variable [@Aliabadi2013Robust].

Amazing Oracles {#sec:oracles}
================

-   Robust convex optimization
    -   oracle technique: affine arithmetic

-   Parametric network potential problem
    -   oracle technique: negative cycle detection

-   Semidefinite programming
    -   oracle technique: Cholesky factorization

Robust Convex Optimization {#sec:robust}
--------------------------

Consider:
$$\begin{array}{ll}
    \text{minimize}   & \sup_{q \in \mathcal Q} f_0(x, q), \\
    \text{subject to} & f_j(x, q) \le 0, \;
            \forall q \in \mathcal{Q}, \; j = 1,2,\cdots, m,
  \end{array}
$$ {#eq:robust-optim}
where $q$ represents a set of varying parameters.
We can reformulate the problem as:
$$\begin{array}{ll}
    \text{minimize}   & t, \\
    \text{subject to} & f_0(x, q) \le t,  \\
                      & f_j(x, q) \le 0, \;
            \forall q \in \mathcal{Q}, \; j = 1,2,\cdots,m.
  \end{array} $$

### Algorithm

The oracle only needs to determine:

-   If $f_j(x_0, q) > 0$ for some $j$ and $q = q_0$, then

-   the cut $(g, \beta)$ = $(\partial f_j(x_0, q_0), f_j(x_0, q_0))$

-   If $f_0(x_0, q) \ge t$ for some $q = q_0$, then

-   the cut $(g, \beta)$ = $(\partial f_0(x_0, q_0), f_0(x_0, q_0) - t)$

-   Otherwise, $x_0$ is feasible, then

-   Let $q_{\max} = \text{argmax}_{q \in \mathcal Q} f_0(x_0, q)$.

-   $t := f_0(x_0, q_{\max})$.

-   The cut $(g, \beta)$ = $(\partial f_0(x_0, q_{\max}), 0)$


### Example: Robust Profit Maximization {#sec:profit-rb}

Consider again the profit maximization problem in @sec:profit. Now suppose that the parameters $\{\alpha, \beta, p, v_1, v_2, k\}$ are subject to interval uncertainties:
$$\begin{array}{rcl}
\alpha - \varepsilon_1 \le & \hat{\alpha} & \le \alpha + \varepsilon_1 \\
\beta  - \varepsilon_2 \le & \hat{\beta}  & \le \beta  + \varepsilon_2 \\
p   - \varepsilon_3 \le  & \hat{p}    & \le p   + \varepsilon_3 \\
v_1 - \varepsilon_4 \le  & \hat{v}_1  & \le v_1 + \varepsilon_4 \\
v_2 - \varepsilon_5 \le  & \hat{v}_2  & \le v_2 + \varepsilon_5 \\
k   - \varepsilon_6 \le  & \hat{k}    & \le k   + \varepsilon_6
\end{array}$$
The problem formulation of the robust counterpart considering the worst-case scenario is:
$$\begin{array}{ll}
    \text{max}  & t \\
    \text{s.t.} & \log(t + \hat{v}_1 e^{y_1} + \hat{v}_2 e^{y_2}) -
                        (\hat{\alpha} y_1 + \hat{\beta} y_2) \le \log(\hat{p}\,A)  \\
                & y_1 \le \log \hat{k}.
  \end{array}$$
In [@Aliabadi2013Robust], the authors present a *piecewise linear approximation* approach. It involves a lot of programming work, but the results are inaccurate. However, this can easily be solved using the cutting-plane method. Note that in this simple example, the worst-case scenario occurs when:

-   $\hat{p} = p - e_3$, $k = \bar{k} - e_3$

-   $v_1 = \bar{v}_1 + e_3$, $v_2 = \bar{v}_2 + e_3$,

-   if $y_1 > 0$, $\alpha = \bar{\alpha} - e_1$, else
    $\alpha = \bar{\alpha} + e_1$

-   if $y_2 > 0$, $\beta = \bar{\beta} - e_2$, else
    $\beta = \bar{\beta} + e_2$

We can even reuse the original oracle to compose the robust counterpart.

```python
class profit_rb_oracle:
    def __init__(self, params, a, v, vparams):
        p, A, k = params
        e1, e2, e3, e4, e5 = vparams
        params_rb = p - e3, A, k - e4
        self.a = a
        self.e = [e1, e2]
        self.P = profit_oracle(params_rb, a, v + e5)

    def __call__(self, y, t):
        a_rb = self.a.copy()
        for i in [0, 1]:
            if y[i] <= 0:
                a_rb[i] += self.e[i]
            else:
                a_rb[i] -= self.e[i]
        self.P.a = a_rb
        return self.P(y, t)
```

Note that the `argmax` may be non-convex and therefore difficult to solve. For more complex problems, one way is to use affine arithmetic for help [@liu2007robust].

Multi-parameter Network Problems {#sec:network}
---------------------------------

Given a network represented by a directed graph $G = (V, E)$.
Consider :
$$\begin{array}{ll}
    \text{minimize} & t, \\
    \text{subject to} & u_i - u_j \le h_{ij}(x, t), \; \forall (i, j) \in E,\\
    \text{variables} &x, u,
  \end{array}$$
where $h_{ij}(x, t)$ is the weight function of edge $(i,j)$.

Assume that the network is large but the number of parameters is small. Given $x$ and $t$, the problem has a feasible solution if and only if $G$ contains no negative cycle. Let $\mathcal{C}$ be a set of all cycles of $G$. We can formulate the problem as:
$$\begin{array}{ll}
    \text{minimize} & t, \\
    \text{subject to} & W_k(x, t) \ge 0, \forall C_k \in C ,\\
       \text{variables} & x,
\end{array} $$
where $C_k$ is a cycle of $G$:
$$W_k(x, t) = \sum_{ (i,j)\in C_k} h_{ij}(x, t).$$

### Negative Cycle Detection Algorithm

The negative cycle detection is the most time-consuming part of the proposed method, so it is very important to choose the proper negative cycle detection algorithm. There are lots of methods to detect negative cycles in a weighted graph [@cherkassky1999negative], in which Tarjan’s algorithm [@Tarjan1981negcycle] is one of the fastest algorithms in practice [@alg:dasdan_mcr; @cherkassky1999negative].

The separation oracle only needs to determine:

-   If there exists a negative cycle $C_k$ under $x_0$, then

-   the cut $(g, \beta)$ = $(-\partial W_k(x_0), -W_k(x_0))$

-   If $f_0(x_0) \ge t$, then the cut $(g, \beta)$ = $(\partial f_0(x_0), f_0(x_0) - t)$.

-   Otherwise, $x_0$ is feasible, then
    -   $t := f_0(x_0)$.
    -   The cut $(g, \beta)$ = $(\partial f_0(x_0), 0)$

### Example: symmetric scalings under the min-max-ratio criterion

This example is taken from [@orlin1985computing]. Given a matrix $A \in \mathbb{R}^{N\times N}$. A *symmetric scaling* of $A$ is a matrix $B$ of the form $U A U^{-1}$ where $U$ is a nonnegative diagonal matrix with the same dimension. According to the *min-max criterion*, the aim is to minimize the largest absolute value of $B$'s elements [@orlin1985computing] (Program 3).

Another possible criterion is to minimize the ratio of largest absolute value of the element $B$ to the smallest. One motivation for using this criterion is that high ratios cause difficulties in performing the simplex method. With this *min-max-ratio* criterion, the symmetric scaling problem can be formulated as [@orlin1985computing] (Program 8):
$$\begin{array}{ll}
    \text{minimize}   &  \pi/\psi  \\
    \text{subject to} &  \psi \le u_i |a_{ij}| u_j^{-1} \le \Pi, \; \forall a_{ij} \neq 0 , \\
                      &  \pi, \psi, u_1 \cdot u_N \, \text{positive}. \\
  \end{array}$$
Let $k’$ denotes $\log( | k | )$. By taking the logarithm of the variables, we can transform the above programming into a two-parameter network optimization problem:
$$\begin{array}{ll}
    \text{minimize}   &  \pi’ - \psi’ \\
    \text{subject to} &  u_i’ - u_j’  \le \pi’ - a_{ij}’, \; \forall a_{ij} \neq 0 \,, \\
                      &  u_j’ - u_i’ \le a_{ij}’ - \psi’, \; \forall a_{ij} \neq 0 \,, \\
    \text{variables}  &  \pi’, \psi’, u’ \, .
  \end{array}$$
where $x = (\pi’, \psi’ )^\mathsf{T}$.
The authors of [@orlin1985computing] claim that they have developed an efficient algorithm for solving such a multi-parameter problem, but we could not find any follow-up publication on this. Interestingly, by using the cutting-plane method, one can easily extend the single-parameter network algorithm to a multi-parameter one.

In this application, $h_{ij}(x)$ is:
$${h}_{ij}(x) = \left\{ \begin{array}{cll}
     -\pi’ + a_{ij}’, & \forall a_{ij} \neq 0 \, ,\\
     \psi’ -a_{ji}’,  & \forall a_{ji} \neq 0 \, ,\\
\end{array} \right.$$
We can find fast algorithms for finding a negative cycle in [@dasdan1998faster; @dasdan2004experimental]. More applications to clock skew scheduling can be found in [@zhou2015multi].

Problems Involving Matrix Inequalities {#sec:lmi}
--------------------------------------

Consider the following problem:
$$\begin{array}{ll}
    \text{find}        & x, \\
    \text{subject to}  & F(x) \succeq 0,
  \end{array}
$$
where $F(x)$ is a matrix-valued function, $A \succeq 0$ denotes $A$ is positive semidefinite.
Recall that a matrix $A$ is positive semidefinite if and only if $v^\mathsf{T} A v \ge 0$ for all $v \in \mathbb{R}^N$.
We can transform the problem into:
$$\begin{array}{ll}
        \text{find}          & x, \\
        \text{subject to}    & v^\mathsf{T} F(x) v \ge 0, \; \forall v \in \mathbb{R}^N.
  \end{array}
$$
Consider $v^\mathsf{T} F(x) v$ is concave for all $v \in \mathbb{R}^N$ w.r.t.
$x$, then the above problem is a convex programming.
Reduce to *semidefinite programming* if $F(x)$ is linear w.r.t.
$x$, i.e., $F(x) = F_0 + x_1 F_1 + \cdots + x_n F_n$.

### Cholesky Factorization Algorithm

An alternative form, eliminating the need to take square roots, is the symmetric indefinite factorization:

$$\begin{aligned}
\mathbf{A} = \mathbf{LDL}^\mathsf{T} & =
\begin{pmatrix}   1 & 0 & 0 \\
   L_{21} & 1 & 0 \\
   L_{31} & L_{32} & 1\\
\end{pmatrix}
\begin{pmatrix}   D_1 & 0 & 0 \\
   0 & D_2 & 0 \\
   0 & 0 & D_3\\
\end{pmatrix}
\begin{pmatrix}   1 & L_{21} & L_{31} \\
   0 & 1 & L_{32} \\
   0 & 0 & 1\\
\end{pmatrix} \\
& = \begin{pmatrix}   D_1 &   &(\mathrm{symmetric})   \\
   L_{21}D_1 & L_{21}^2D_1 + D_2& \\
   L_{31}D_1 & L_{31}L_{21}D_{1}+L_{32}D_2 & L_{31}^2D_1 + L_{32}^2D_2+D_3.
\end{pmatrix}.
\end{aligned} $$

If $A$ is real, the following recursive relations apply for the entries of $D$ and $L$:

$$ D_{j} = A_{jj} - \sum_{k=1}^{j-1} L_{jk}L_{jk}^* D_k, $$
$$ L_{ij} = \frac{1}{D_j} \left( A_{ij} - \sum_{k=1}^{j-1} L_{ik} L_{jk}^* D_k \right) \quad \text{for } i>j.
$$

Again, the pattern of access allows the entire computation to be performed in-place if desired.

The following is the algorithm written in Python:

```python
def factor(self, getA):
    T = self.T
    for i in range(self.n):  # from 0 to n-1
        for j in range(i+1): # from 0 to i
            d = getA(i, j) - np.dot(T[:j, i], T[j, :j])
            T[i, j] = d
            if i != j:
                T[j, i] = d / T[j, j]
        if d <= 0.:  # strictly positive
            self.p = i
            return
    self.p = self.n
```

The vector $v$ can be found.
The following is the algorithm written in Python:

```python
def witness(self):
    p = self.p
    n = p + 1
    v = np.zeros(n)
    v[p] = 1
    for i in range(p, 0, -1): # backward substitution
        v[i-1] = -np.dot(self.T[i-1, i:n], v[i:n])
    return v, -self.T[p, p]
```

The oracle only needs to:

-   Perform a *row-based* Cholesky factorization such that
    $F(x_0) = R^\mathsf{T} R$.

-   Let $A_{:p,:p}$ denotes a submatrix
    $A(1:p, 1:p) \in \mathbb{R}^{p\times p}$.

-   If Cholesky factorization fails at row $p$,
    -   there exists a vector
        $e_p = (0, 0, \cdots, 0, 1)^\mathsf{T} \in \mathbb{R}^p$, such that
        -   $v = R_{:p,:p}^{-1} e_p$, and
        -   $v^\mathsf{T} F_{:p,:p}(x_0) v < 0$.
    -   The cut $(g, \beta)$ =
        $(-v^\mathsf{T} \partial F_{:p,:p}(x_0) v, -v^\mathsf{T} F_{:p,:p}(x_0) v)$

### Example: Matrix Norm Minimization

Let $A(x) = A_0 + x_1 A_1 + \cdots + x_n A_n$.
Problem $\min_x \| A(x) \|$ can be reformulated as
$$\begin{array}{ll}
    \text{minimize}      & t, \\
    \text{subject to}    & \begin{pmatrix}
                             t\,I_m   & A(x) \\
                             A^\mathsf{T}(x) & t\,I_n
                            \end{pmatrix} \succeq 0.
  \end{array}
$$
A binary search on $t$ can be used for this problem.

### Example: Estimation of Correlation Function

Random Field [@Schabenberger05]
-------------------------------

*Random field*, also known as *stochastic process*, can be regarded as an indexed family of random variables denoted as {$Z(\mathbf{s}): \mathbf{s}\in D$}, where $D$ is a subset of $d$-dimensional Euclidean space $\mathbb{R}^d$. To specify a stochastic process, the joint probability distribution function of any finite subset $(Z(\mathbf{s}_1), \ldots, Z(\mathbf{s}_n))$ must be given in a consistent way, which is called *distribution* of the process. For ease of analysis, a random field is often assumed to be with *Gaussian* distribution and is called Gaussian random field.

A random field has several key properties useful in practical problems. The field is *stationary* under translations, or *homogeneous*, if the distribution is unchanged when the point set is translated. The field is *isotropic* if the distribution is invariant under any rotation of the whole points in the parameter space. We study the homogeneous isotropic field in this paper.

The *covariance* $C$ and *correlation* $R$ of a stochastic process are defined by:
$$C(\mathbf{s}_i,\mathbf{s}_j) = \mathrm{cov}(Z(\mathbf{s}_i),Z(\mathbf{s}_j)) = \mathrm{E}\lbrack (Z(\mathbf{s}_i)-\mathrm{E}\lbrack Z(\mathbf{s}_i)\rbrack)(Z(\mathbf{s}_j)-\mathrm{E}\lbrack Z(\mathbf{s}_j)\rbrack)\rbrack $$
and
$$R(\mathbf{s}_i,\mathbf{s}_j)=C(\mathbf{s}_i,\mathbf{s}_j)/ \sqrt{C(\mathbf{s}_i,\mathbf{s}_i)C(\mathbf{s}_j,\mathbf{s}_j)} $$
respectively for all $\mathbf{s}_i,\mathbf{s}_j\in D$, where $\mathrm{E}\lbrack Z(\mathbf{s})\rbrack$ denotes the expectation of $Z(\mathbf{s})$. Thus a process is homogeneous if $C$ and $R$ depend
only on the separation vector $\mathbf{h}=\mathbf{s}_i-\mathbf{s}_j$. Furthermore, it is isotropic if $C$ and $R$ depend upon $\mathbf{h}$ only through its length $h$, i.e.,
$$C(\mathbf{s}_i,\mathbf{s}_j)=C(\mathbf{h})=C(h),
$$

$$R(\mathbf{s}_i,\mathbf{s}_j)=R(\mathbf{h})=R(h)=C(h)/C(0).$$ 
{#eq:corr_def}
If we denote $C(0)$, the variance of $Z(\mathbf{s})$, as $\sigma^2$, then the relationship between covariance and correlation is $C(h)=\sigma^2 R(h)$.

When the two components are considered, the measurement data can still be regarded as a Gaussian random field, but the correlation function will have a discontinuity at the origin. We call this phenomenon "nugget effect" [@Diggle07].

$$\begin{array}{ll}
   \min_{\kappa, p}   & \| \Omega(p) + \kappa I - Y \| \\
   \text{s.
t.} & \Omega(p) \succcurlyeq 0,  \kappa \ge 0 \; .\\
  \end{array}
$$
Let $\rho(h) = \sum_i^n p_i \Psi_i(h)$, where $p_i$’s are the unknown coefficients to be fitted $\Psi_i$’s are a family of basis functions. The covariance matrix $\Omega(p)$ can be recast as:
$$\Omega(p) = p_1 F_1 + \cdots + p_n F_n, $$
where $\{F_k\}_{i,j} =\Psi_k( \| s_j - s_i \|_2)$.

Ellipsoid Method Revisited {#sec:ellipsoid}
===========================================

Some History of the Ellipsoid Method [@BGT81]. Introduced by Shor and Yudin and Nemirovskii in 1976. It used to show that linear programming (LP) is polynomial-time solvable (Kachiyan 1979), settled the long-standing problem of determining the theoretical complexity of LP. In practice, however, the simplex method runs much faster than the method, although its worst-case complexity is exponential.

Basic Ellipsoid Method
-----------------------

An ellipsoid $\mathcal{E}_k(x_k, P_k)$ is specified as a set
$$\{x \mid (x-x_k) P^{-1}_k (x - x_k) \le 1 \}, $$
where $x_k \in \mathbb{R}^n$ is the center of the ellipsoid and $P_k \in \mathbb{R}^{n \times n}$ is a positive definite matrix.

\begin{figure}
\centering
\begin{tikzpicture}[scale=0.6]
 \draw[top color=lightgray, bottom color=lightgray] plot[smooth, tension=.7] coordinates {(-3,2) (-5,2) (-6,4) (-5,5) (-3,4) (-3,2)};
 \node at (-5,4) {$\mathcal{K}$};
\draw (0,8) -- (-3,-2);
\draw [fill=qqqqff] (-1,3) circle (1.5pt)
   node [above right] {$x_k$};
\draw  (-1,3) ellipse (7 and 3);
  \node at (5,4) {$\mathcal{E}_k$};
\end{tikzpicture}
\caption{Ellipsoid, feasible region, and cut}
\label{fig:ellipsoid}
\end{figure}

Updating the ellipsoid (deep-cut)

Calculation of minimum volume ellipsoid covering:
$$\mathcal{E}_k \cap \{z \mid g^\mathsf{T} (z - x_k) + \beta \le 0 \}$$
Let $\tilde{g} = P_k\,g$, $\tau^2 = g^\mathsf{T} P_k g$. We can make the following observations:

1. If $n \cdot \beta < -\tau$ (shallow cut), then no smaller ellipsoid can be found.

2. If $\beta > \tau$, then intersection is empty.

3. Otherwise,
   $$x_c^+ = x_c - \frac{\rho}{ \tau^2 } \tilde{g}, \qquad
     P^+ = \delta\cdot\left(P - \frac{\sigma}{ \tau^2 } \tilde{g}\tilde{g}^\mathsf{T}\right)
   $$
   where
   $$\rho = \frac{ \tau+nh}{n+1}, \qquad
     \sigma = \frac{2\rho}{ \tau+\beta}, \qquad
     \delta = \frac{n^2(\tau^2 - \beta^2)}{(n^2 - 1)\tau^2} $$

Even better, split $P$ into two variables $\kappa \cdot Q$. Let $\tilde{g} = Q \cdot g$, $\omega = g^\mathsf{T}\tilde{g}$, $\tau = \sqrt{\kappa\cdot\omega}$.
$$x_c^+ = x_c - \frac{\rho}{\omega} \tilde{g}, \qquad
  Q^+ = Q - \frac{\sigma}{\omega} \tilde{g}\tilde{g}^\mathsf{T}, \qquad
  \kappa^+ =  \delta\cdot\kappa $$
Reduce $n^2$ multiplications per iteration. Note that:

-   The determinant of $Q$ decreases monotonically.

-   The range of $\delta$ is $(0, \frac{n^2}{n^2 - 1})$

Central Cut
-----------

A Special case of  when $\beta = 0$. Deserve a separate implement because it is much simpler. Let $\tilde{g} = Q\,g$, $\tau = \sqrt{\kappa\cdot\omega}$,

$$\rho = \frac{\tau}{n+1}, \qquad
  \sigma = \frac{2}{n+1}, \qquad
  \delta = \frac{n^2}{n^2 - 1}. $$

Parallel Cuts
-------------

Oracle returns a pair of cuts instead of just one. The pair of cuts is given by $g$ and $(\beta_1, \beta_2)$ such that:
$$\begin{array}{l}
    g^\mathsf{T} (x - x_c) + \beta_1 \le 0,  \\
    g^\mathsf{T} (x - x_c) + \beta_2 \ge 0,
  \end{array}$$
for all $x \in \mathcal{K}$.

Only linear inequality constraint can produce such parallel cut:
$$ l \le a^\mathsf{T} x + b \le u, \qquad L \preceq F(x) \preceq U.$$

Usually, provide faster convergence.

![Parallel cuts](ellipsoid.files/parallel_cut.pdf){#fig:parallel_cut width="80%"}

Updating the ellipsoid.

Let $\tilde{g} = Q\,g$, $\tau^2 = \kappa\cdot\omega$.

-   If $\beta_1 > \beta_2$, intersection is empty.

-   If $\beta_1 \beta_2 < -\tau^2/n$, no smaller ellipsoid can be found.

-   If $\beta_2^2 > \tau^2$, it reduces to deep-cut with $\alpha = \alpha_1$.

Otherwise,
$$x_c^+ = x_c - \frac{\rho}{\omega} \tilde{g}, \qquad
    Q^+ = Q - \frac{\sigma}{\omega} \tilde{g}\tilde{g}^\mathsf{T}, \qquad
    \kappa^+ =  \delta \kappa.
$$
where
$$\begin{array}{lll}
      \bar{\beta} &=& (\beta_1 + \beta_2)/2 \\
      \xi^2 &=& (\tau^2 - \beta_1^2)(\tau^2 - \beta_2^2) + (n(\beta_2 - \beta_1)\bar{\beta})^2, \\
      \sigma &=& (n + (\tau^2 - \beta_1\beta_2 - \xi)/(2\bar{\beta}^2)) / (n + 1), \\
      \rho &=& \bar{\beta}\cdot\sigma, \\
      \delta &=& (n^2/(n^2-1)) (\tau^2 - (\beta_1^2 + \beta_2^2)/2 + \xi/n) / \tau^2 .
\end{array}
$$

### Example: FIR filter design

A typical structure of digital Finite Impulse Response (FIR) filter is shown in @fig:fir-strctr, where the coefficients $h[0], h[1], \ldots, h[n-1]$ must be determined to meet given specifications. Usually, they can be manually designed using windowing or frequency-sampling techniques [@oppenheim1989discrete].

However, the experience and knowledge of designers are highly demanded in this kind of design methods. Moreover, there is no guarantee about the design’s quality. Therefore, the optimization-based techniques (e.g. [@wu1999fir], more reference) have attracted tons of research effort. In this kind of method, facilitated with growing computing resources and efficient optimization algorithms, the solution space can be effectively explored.

![A typical structure of an FIR filter\ @mitra2006digital.](ellipsoid.files/fir_strctr.pdf){#fig:fir-strctr width="80%"}

In optimization algorithms, what is particularly interesting is the convex optimization. If a problem is in a convex form, it can be efficiently and optimally solved. Convex optimization techniques are also implementable in designing FIR filters, including the Parks-McClellan algorithm [@park1972chebyshev], METEOR [@steiglitz1992meteor], and peak-constrained least-squares (PCLS) [@selesnick1996constrained; @adams1998peak]. In the mentioned articles, with the help of exchange algorithms (e.g. Remez exchange algorithm), certain FIR filter design problems can be formed as linear or quadratic programs. They are two simple forms of convex optimization problems, which can be optimally solved with existing algorithms, such as the interior-point method [@boyd2009convex]. Tempted by the optimality, more efforts were devoted to forming the problem convex. Particularly, in [@wu1999fir], via spectral factorization [@goodman1997spectral], the problem of designing an FIR filter with magnitude constraints on frequency-domain is formulated as a convex optimization problem. More examples are provided in [@davidson2010enriching].

Its time response is
$$y[t] = \sum_{k=0}^{n-1}{h[k]u[t-k]}$$
{#eq:t_res}
where $\mathbf{h} = (h(0), h(1),..., h(n-1))$ is the filter coefficients. Its frequency response $H: [0,\pi] \rightarrow \mathbb{C}$ is
$$H(\omega) = \sum_{m=0}^{n-1}{h(m)e^{-jm\omega}}$$ 
{#eq:f_res}
where $j = \sqrt{-1}$, $n$ is the order of the filter.
The design of a filter with magnitude constraints is often formulated as a constraint optimization problem as the form
$$\begin{aligned}
  \min            &  \gamma \\
  \mathrm{s.t.}   &  f(\mathbf{x}) \le \gamma \\
                  &  g(\mathbf{x}) \le 0.\end{aligned}$$
{#eq:ori}
where $\mathbf{x}$ is the vector of design variables, $g(\mathbf{x})$ represents the characteristics of the desirable filter and $f(\mathbf{x})$ is the performance metric to be optimized. For example, the magnitude constraints on frequency domain are expressed as
$$L(\omega) \le |H(\omega)| \le U(\omega), \forall \omega\in(-\infty,+\infty)$$
{#eq:mag_cons}
where $L(\omega)$ and $U(\omega)$ are the lower and upper (nonnegative) bounds at frequency $\omega$ respectively. Note that $H(\omega)$ is $2\pi$ periodic and $H(\omega)=\overline{H(-\omega)}$.
Therefore, we can only consider the magnitude constraint on $[0,\pi]$ [@wu1999fir].

Generally, the problem might be difficult to solve, since we can only obtain the global optimal solution with resource-consuming methods, such as branch-and-bound [@davidson2010enriching]. However, the situation is totally different if the problem is convex, where $f(\mathbf{x})$ and $g(\mathbf{x})$ are convex functions. In such a case, the problem can be optimally solved with many efficient algorithms.

Attracted by the benefits, the authors of\ [@wu1999fir] transformed (?), originally non-convex, into a convex form via spectral factorization:

$$L^2(\omega) \le R(\omega) \le U^2(\omega), \forall \omega\in(0,\pi)$$ {#eq:r_con}
where $R(\omega)=\sum_{i=-n+1}^{n-1}{r(t)e^{-j{\omega}t}}=|H(\omega)|^2$ and $\mathbf{r}=(r(-n+1),r(-n+2),\ldots,r(n-1))$ are the autocorrelation coefficients. Especially, $\mathbf{r}$ can be determined by $\mathbf{h}$, with the following equation vice versa\ [@wu1999fir]:

$$r(t) = \sum_{i=-n+1}^{n-1}{h(i)h(i+t)}, t\in\mathbb{Z}.$$ {#eq:h_r}
where $h(t)=0$ for $t<0$ or $t>n-1$.

![Result](ellipsoid.files/lowpass.pdf){#fig:lowpass width="80%"}

### Example: Maximum Likelihood estimation

Consider
$$\begin{array}{ll}
    \min_{\kappa, p}  & \log\det(\Omega(p) + \kappa\cdot I) +
                \mathrm{Tr}((\Omega(p) + \kappa\cdot I)^{-1}Y), \\
    \text{s.t.}       & \Omega(p) \succeq 0, \kappa \ge 0 .
\\
  \end{array}
$$
Note that the first term is concave, the second term is convex. However, if there are enough samples such that $Y$ is a positive definite matrix, then the function is convex within $[0, 2Y]$.
Therefore, the following problem is convex:
$$\begin{array}{ll}
    \min_{\kappa, p}  & \log\det V(p) + \mathrm{Tr}(V(p)^{-1}Y),\\
    \text{s.t.}       & \Omega(p) + \kappa \cdot I = V(p) \\
                      & 0 \preceq V(p) \preceq 2Y, \kappa {>} 0.
  \end{array}$$

Discrete Optimization
---------------------

Many engineering problems can be formulated as a convex/geometric programming, such as digital circuit sizing. However, in ASIC design, there is often only a limited set of cell types to choose from the cell library. In other words, some design variables are discrete. We can express the discrete version as Mixed-Integer Convex programming (MICP) by mapping the design variables to integers.

What are the problems with the existing methods? It is mostly based on relaxation. The relaxed solution is then used as the lower bound and the branch-and-bound method is used to find the discrete optimal solution. Note that the branch-and-bound method does not exploit the convexity of the problem. What if only constraints on discrete data can be evaluated?

A relaxed optimal solution (convex) is usually obtained first. Then the optimized discrete solution is obtained by exhausting the neighborhood search. However, sometimes the constraints are tight so that the relaxed continuous optimal solution is far from the discrete one. Enumeration of the discrete domains is difficult.

Consider:
$$\begin{array}{ll}
        \text{minimize}      & f_0(x), \\
        \text{subject to}    & f_j(x) \le 0, \; \forall j=1,2,\ldots, \\
                             & x \in \mathbb{D},
  \end{array}$$
where $f_0(x)$ and $f_j(x)$ are “convex”. Note that some design variables are discrete. The oracle looks for a nearby discrete solution $x_d$ of $x_c$ with the cutting-plane:
$$ g^\mathsf{T} (x - x_d) + \beta \le 0, \beta \ge 0, g \neq 0. $$
Note that the cut may be a shallow cut.
Suggestion: use as many different cuts as possible for each iteration (e.g. round-robin the evaluation of constraints).

### Example: Multiplierless FIR Filter Design

However, there are still many filter design problems that are non-convex, such as multiplierless FIR filter design problems. Note that in [@fig:fir-strctr], each coefficient associated with a multiplier unit makes the filter power-hungry, especially in *application specific integrated circuits* (ASIC). Fortunately, if each coefficient is quantized and represented as a sum of Singed Power-of-Two (SPT), a multiplierless filter can be implemented. Such coefficients can be uniquely represented by a Canonical Signed-Digit (CSD) code with a minimum number of non-zero digits[@george1960csd]. In this case, it confines the multiplication to addition and shift operations. An example is shown in @fig:multi-shift. The coefficient 0.40625 = 13/32 can be written as $2^{-1} - 2^{-3} + 2^{-5}$. Thus, as shown in @fig:shift, the multiplier can be replaced with three shifters and two adders at a much lower cost. However, the coefficient quantization constraint is non-convex, making the convex optimization algorithm not directly applicable. A similar case is the consideration of the finite word-length effect [@lim1982finite].

Attracted by the benefits of this "multiplier-free" approach, many efforts have been devoted to its design techniques. For its general problems, integer programming (e.g. [@kodek1980design; @lim1982finite; @lim1983fir; @lim1999signed]) can be implemented to achieve the optimal solution. However, it requires excessive computational resources. Other heuristic techniques, such as genetic algorithm [@xu1995design] and dynamic-programming-like method [@chen1999trellis], also have inefficiency. If the quantization constraint is the only non-convex constraint in the design problem, a lower bound can be efficiently obtained by solving the relaxed problem [@davidson2010enriching]. Then to make the solution feasible, it can be rounded to the nearest CSD code or used as a starting point of a local search algorithm to obtain a better solution [@kodek1981comparison]. However, neither method guarantees the feasibility of the final solution. Besides, the local search problem remains non-convex. Therefore, the adopted algorithm may also be inefficient, such as branch-and-bound in [@kodek1981comparison].

![Result](ellipsoid.files/csdlowpass.pdf){#fig:csdlowpass width="80%"}

Concluding Remarks
==================

Should be known to students.
The ellipsoid method is not a competitor but a companion of interior-point methods.

TBD.

References {-}
==========

\ 

