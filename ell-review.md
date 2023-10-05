---
title: Ellipsoid Method and the Amazing Oracles
bibliography: ['ellipsoid.bib', 'fir-ref.bib', 'Geostatistics.bib', 'mpcss1.bib', 'mpcss2.bib']
csl: 'applied-mathematics-letters.csl'
abstract: |
  The ellipsoid method is an optimization technique that offers distinct advantages over interior-point methods because it does not require evaluating all constraint functions. This makes it ideal for problems with many or even infinite constraints. The method utilizes an ellipsoid as a search space and employs a separation oracle to provide a cutting plane for updating the search space. Notably, the importance of the separation oracle is often overlooked. This article evaluates the usage of the ellipsoid method in three specific applications, including robust convex optimization, semidefinite programming, and parametric network optimization. The effectiveness of separation oracles is analyzed for each application. We also discuss implementation issues of the ellipsoid method, such as utilizing parallel cuts to update the search space and enhance computation time. In some instances, parallel cuts can drastically reduce computation time, as observed in FIR filter design. Discrete optimization is also investigated, illustrating how the ellipsoid method can be applied to problems that involve discrete design variables. An oracle implementation is required solely for locating the nearest discrete solutions

---

# Introduction

The ellipsoid method's reputation suffers due to its perceived slowness when solving large-scale convex problems in comparison to interior-point methods. This perception is unfair. Unlike the interior-point method, the ellipsoid method does not need explicit evaluation of all constraint functions. Instead, it utilizes an ellipsoid as a search space and only requires a separation oracle that provides a *cutting plane* (@sec:cutting_plane). The method is well-suited for problems that involve a moderate number of design variables but have many or even infinite constraints. Some criticize the method, claiming it is unable to leverage sparsity. However, although the ellipsoid method cannot take advantage of the sparsity of the problem, the separation oracle is capable of taking advantage of certain structural types.

Despite decades of investigation into the ellipsoid method [@BGT81], the importance of the separation oracle is often overlooked. In this article, we examine three specific applications: robust convex optimization, network optimization, and semidefinite programming. The effectiveness of separation oracles is analyzed for each application. 

Robust optimization incorporates parameter uncertainties into the optimization problem by analyzing the worst-case scenario. The goal is to find a strong solution that performs optimally under various possible parameter values in a given set of uncertainties. A robust counterpart of a convex problem preserves its convexity, although the number of constraints grows infinitely. This makes the ellipsoid method an excellent choice for tackling such problems. This is detailed in @sec:robust.

An example of network optimization where the ellipsoid method can be employed is also presented. The separation oracle involves constructing a cutting plane by identifying a negative cycle in a network graph. There are algorithms available that utilize network locality and other properties, resulting in an effective oracle implementations. This is discussed in more detail in @sec:network.

Meanwhile, @sec:lmi describes concerns surrounding matrix inequalities. Remember that utilizing $LDL^\mathsf{T}$ decomposition can efficiently check the positive definiteness of a symmetric matrix. If a symmetric matrix $A$ with dimensions $m \times m$ encounters a non-positive diagonal entry during decomposition that causes the process to stop at row $p$, $A$ cannot be positive definite. At the same time, a witness vector $v$ can be constructed to certify that $v^\mathsf{T} A v \leq 0$. With the row-based style decomposition and a lazy evaluation technique, the cutting plane can be constructed exactly in $O(p^3)$, enabling its use for efficient oracle implementations.

The implementation of the ellipsoid method is discussed in @sec:ellipsoid. This method is a cutting plane approach where the *search space* is an ellipsoid, conventionally represented as 
  $$\{x \mid (x-x_c)P^{-1}(x-x_c) \le 1\},$$ 
where $x_c \in \mathbb{R}^n$ is the center of the ellipsoid. The matrix $P \in \mathbb{R}^{n \times n}$ is symmetric positive definite. During each iteration, the ellipsoid method updates $x_c$ and $P$. While updating ellipsoids can be straightforward and has been implemented for decades, we show that the cost can be further reduced by $n^2$ floating point operations by multiplying $\kappa$ by $Q$ and splitting $P$, resulting in this form: 
  $$\{ x \mid (x-x_c)Q^{-1}(x-x_c) \le \kappa \}.$$
In addition, @sec:parallel_cut discusses the use of parallel cuts. Two constraints can be used simultaneously to update the ellipsoid when a pair of parallel inequalities occur, one of which is violated. Some articles suggest that this method does not lead to significant improvements. However, we show that in situations where certain constraints have tight upper and lower bounds, such as in some filter designs, the implementation of parallel cuts can significantly speed up the runtime. In addition, we show that if the method is implemented carefully, any update, whether it uses a deep cut or a parallel cut, results in at most one square root operation.

In many practical engineering problems, some design variables may be restricted to discrete forms. Since the cutting-plane method requires only a separation oracle, it can also be used for discrete problems...

Cutting-plane Method Revisited {#sec:cutting_plane}
==============================

Convex Feasibility Problem
--------------------------

Let $\mathcal{K} \subseteq \mathbb{R}^n$ be a convex set. Consider the feasibility problem:

1.  find a point $x^* \in \mathbb{R}^n$ in $\mathcal{K}$, or
2.  determine that $\mathcal{K}$ is empty (i.e., has no feasible solution).

A separation oracle, also known as a cutting-plane oracle, is a concept in the mathematical theory of convex optimization¹. It is a method used to describe a convex set that is given as an input to an optimization algorithm¹. Separation oracles are particularly used as input to ellipsoid methods¹.

The separation oracle operates in the following way: Given a vector `y` in `R^n`, it returns one of the following¹:

1. Asserts that `y` is in `K`, where `K` is a convex and compact set in `R^n`.
2. Finds a hyperplane that separates `y` from `K`: a vector `a` in `R^n`, such that `a.y > a.x` for all `x` in `K`.

This oracle can be implemented in polynomial time as long as the number of constraints is polynomial¹. It's important to note that a strong separation oracle is completely accurate, and thus may be hard to construct. For practical reasons, a weaker version is considered, which allows for small errors in the boundary of `K` and the inequalities¹.

Source: Conversation with Bing, 5/10/2023
(1) Separation oracle - Wikipedia. https://en.wikipedia.org/wiki/Separation_oracle.
(2) Lecture # 8: Separation oracles and the ellipsoid algorithm. https://courses.cs.duke.edu/compsci530/cps230/current/notes/lec8.pdf.
(3) Solving convex programs defined by separation oracles?. https://or.stackexchange.com/questions/2899/solving-convex-programs-defined-by-separation-oracles.
(4) Introduction to Optimization Theory - Stanford University. https://web.stanford.edu/~sidford/courses/20fa_opt_theory/sidford_2020fa_mse213_cs269o_lec18.pdf.
(5) Separation oracle - Wikiwand. https://www.wikiwand.com/en/Separation_oracle.

When a *separation oracle* $\Omega$ is *queried* at $x_0$, it either

1.  asserts that $x_0 \in \mathcal{K}$, or
2.  returns a separating hyperplane between $x_0$ and $\mathcal{K}$:

$$g^\mathsf{T} (x - x_0) + \beta \le 0, \beta \ge 0, g \neq 0, \; \forall x \in \mathcal{K}.$$

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

Here's a basic outline of how the cutting-plane method works:

1. **Initialization**: Start with a search space $\mathcal{S}$ that is guaranteed to contain a minimizer of the function.

2. **Iteration**: Denote the center of the current $\mathcal{S}$ as $x_c$. In each iteration, query the separation oracle at $x_c$. Compute a subgradient of the function at the center of the current search space. This gives us a half-space that is guaranteed to contain a minimizer.

3. **Update**: Compute the smaller $\mathcal{S}^+$ that contains the half-space from step 2. This new search space is guaranteed to contain a minimizer of the function.

4. **Repeat**: Repeat steps 2 and 3 until $\mathcal{S}$ is empty or it is small enough.

For example, consider we are minimizing a convex function `f` in `R^n`. We start with an ellipsoid `E(k)` that is guaranteed to contain a minimizer of `f`. We compute a subgradient `g(k)` of `f` at the center, `x(k)`, of `E(k)`. We then know that the half ellipsoid `E(k) ∩ {z | g(k)T (z − x(k)) ≤ 0}` contains a minimizer of `f`. We compute the ellipsoid of minimum volume that contains this sliced half ellipsoid; this new ellipsoid `E(k+1)` is then guaranteed to contain a minimizer of `f`².


From Feasibility to Optimization
--------------------------------

Consider:
$$\begin{array}{ll}
    \text{minimize}     & f_0(x), \\
    \text{subject to}   & x \in \mathcal{K}.
  \end{array}$$
We treat the optimization problem as a feasibility problem with an additional constraint $f_0(x) \le t$.
Here, $f_0(x)$ can be a convex function or a quasi-convex function.
$t$ is the best-so-far value of $f_0(x)$.
We can reformulate the problem as:
$$\begin{array}{ll}
    \text{minimize}   & t, \\
    \text{subject to} & \Phi(x, t) \le 0, \\
                      & x \in \mathcal{K},
  \end{array} $$
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

```python
def cutting_plane_optim(omega, space, t, options=Options()):
    x_best = None
    for niter in range(options.max_iters):
        cut, t1 = omega.assess_optim(space.xc(), t)
        if t1 is not None:  # better t obtained
            t = t1
            x_best = copy.copy(space.xc())
            status = space.update_central_cut(cut)
        else:
            status = space.update_deep_cut(cut)
        if status != CutStatus.Success or space.tsq() < options.tol:
            return x_best, t, niter
    return x_best, t, options.max_iters
```

Example: Profit Maximization {#sec:profit}
------------------------------------------

This example is taken from [@Aliabadi2013Robust]. Consider the following *short-run* profit maximization problem:
$$\begin{array}{ll}
   \text{maximize}   & p(A x_1^\alpha x_2^\beta) - v_1 x_1 - v_2 x_2, \\
   \text{subject to} & x_1 \le k, \\
                     & x_1 > 0, x_2 > 0,
  \end{array}$$ {#eq:profit-max-in-original-form}
where $p$ is the market price per unit, $A$ is the scale of production, $\alpha$ and $\beta$ are output elasticities, $x_i$ and $v_i$ are the i-th input quantity and output price, $A x_1^\alpha x_2^\beta$ is the Cobb-Douglas production function which is a widely accepted model used to represent the relationship between inputs and outputs in production. $k$ is a constant that limits the quantity of $x_1$. The above formulation is not in convex form. First, we rewrite the problem:

$$\begin{array}{ll}
    \text{maximize}   & t, \\
    \text{subject to} & t  + v_1 x_1  + v_2 x_2 \le p A x_1^{\alpha} x_2^{\beta}, \\
                      & x_1 \le k, \\
                      & x_1 > 0, x_2 > 0.
  \end{array} $$
By the change of variables, we can obtain the following convex form of\ @eq:profit-max-in-original-form:
$$\begin{array}{ll}
    \text{maximize}   & t, \\
    \text{subject to} & \log(t + v_1 e^{y_1} + v_2 e^{y_2}) -
                    (\alpha y_1 + \beta y_2) \le \log(p\,A), \\
                      & y_1 \le \log k,
  \end{array}$$
{#eq:profit-in-cvx-form}
where $y_1 = \log x_1$ and $y_2 = \log x_2$. 


```python
class ProfitOracle(OracleOptim):
    def __init__(self, params, elasticities, price_out):
        unit_price, scale, limit = params
        self.log_pA = math.log(unit_price * scale)
        self.log_k = math.log(limit)
        self.price_out = price_out
        self.elasticities = elasticities

    def assess_optim(self, y, t):
        if (fj := y[0] - self.log_k) > 0.0:  # constraint
            grad = np.array([1.0, 0.0])
            return (grad, fj), None

        log_Cobb = self.log_pA + self.elasticities.dot(y)
        q = self.price_out * np.exp(y)
        vx = q[0] + q[1]
        if (fj := math.log(t + vx) - log_Cobb) >= 0.0:
            grad = q / (t + vx) - self.elasticities
            return (grad, fj), None

        t = np.exp(log_Cobb) - vx
        grad = q / (t + vx) - self.elasticities
        return (grad, 0.0), t
```

Some readers may recognize that we can also write the problem in a geometric program by introducing one additional variable [@Aliabadi2013Robust].




Amazing Oracles {#sec:oracles}
================

-   Robust convex optimization
    -   oracle technique: affine arithmetic

-   Parametric network potential problem
    -   oracle technique: negative cycle detection

-   Semidefinite programming
    -   oracle technique: Cholesky decomposition

Robust Convex Optimization {#sec:robust}
--------------------------

Overall, robust optimization accounts for parameter uncertainties by formulating problems that consider worst-case scenarios. This approach allows for more reliable and robust solutions when dealing with uncertainty. The study presented in this paper addresses profit maximization using a robust geometric programming approach with interval uncertainty. The authors study the well-established Cobb-Douglas production function and introduce an approximate equivalent of the robust counterpart using piecewise convex linear approximations. This approximation takes the form of a geometric programming problem. An example is used to demonstrate the impact of uncertainties.


Interval uncertainties refer to uncertainties in model parameters that are represented as intervals. The authors of this study consider interval uncertainties in the model parameters. Upper and lower piecewise convex linear approximations of the robust counterpart are presented that are efficiently solvable using interior point methods. These approximations are used to incorporate the interval uncertainties into the model.

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

Consider again the profit maximization problem in @sec:profit. Uncertainties in the model parameters over a given interval. Now suppose that the parameters $\{\alpha, \beta, p, v_1, v_2, k\}$ are subject to interval uncertainties:
[@Aliabadi2013Robust].
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

In [@Aliabadi2013Robust], the authors propose the use of piecewise convex linear approximations as a close approximation of the robust counterpart, leading to more solvability using interior-point algorithms. This requires a lot of programming, but the results are imprecise. However, this can be easily solved using the cutting plane method. Note that in this simple example, the worst case happens when:

-   $\hat{p} = p - e_3$, $k = \bar{k} - e_3$
-   $v_1 = \bar{v}_1 + e_3$, $v_2 = \bar{v}_2 + e_3$,
-   if $y_1 > 0$, $\alpha = \bar{\alpha} - e_1$, else
    $\alpha = \bar{\alpha} + e_1$
-   if $y_2 > 0$, $\beta = \bar{\beta} - e_2$, else
    $\beta = \bar{\beta} + e_2$

We can even reuse the original oracle to compose the robust counterpart.

```python
class ProfitRbOracle(OracleOptim):
    def __init__(self, params, elasticities, price_out, vparams):
        e1, e2, e3, e4, e5 = vparams
        self.elasticities = elasticities
        self.e = [e1, e2]
        unit_price, scale, limit = params
        params_rb = unit_price - e3, scale, limit - e4
        self.omega = ProfitOracle(
            params_rb, elasticities, price_out + np.array([e5, e5])
        )

    def assess_optim(self, y, t):
        a_rb = copy.copy(self.elasticities)
        for i in [0, 1]:
            a_rb[i] += -self.e[i] if y[i] > 0.0 else self.e[i]
        self.omega.elasticities = a_rb
        return self.omega.assess_optim(y, t)
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

The minimum cycle ratio (MCR) problem is a fundamental problem in the analysis of directed graphs. Given a directed graph, the MCR problem seeks to find the cycle with the minimum ratio of the sum of the edge weights to the number of edges in the cycle. In other words, the MCR problem tries to find the "tightest" cycle in the graph, where the tightness of a cycle is measured by the ratio of the total weight of the cycle to its length.

The MCR problem has many applications in the analysis of discrete event systems, such as digital circuits and communication networks. It is closely related to other problems in graph theory, such as the shortest path problem and the maximum flow problem. Efficient algorithms for solving the MCR problem are therefore of great practical importance.

### Negative Cycle Detection Algorithm

The negative cycle detection is the most time-consuming part of the proposed method, so it is very important to choose the proper negative cycle detection algorithm. There are lots of methods to detect negative cycles in a weighted graph [@cherkassky1999negative], in which Tarjan’s algorithm [@Tarjan1981negcycle] is one of the fastest algorithms in practice [@alg:dasdan_mcr; @cherkassky1999negative].

Howard's method is a minimum cycle ratio (MCR) algorithm that uses a policy iteration algorithm to find the minimum cycle ratio of a directed graph. The algorithm maintains a set of candidate cycles and iteratively updates the cycle with the minimum ratio until convergence. 

To detect negative cycles, Howard's method uses a cycle detection algorithm based on the Bellman-Ford algorithm. Specifically, the algorithm maintains a predecessor graph of the original graph and performs cycle detection on this graph using the Bellman-Ford algorithm. If a negative cycle is detected, the algorithm stops and returns the cycle.

The separation oracle only needs to determine:

-   If there exists a negative cycle $C_k$ under $x_0$, then
-   the cut $(g, \beta)$ = $(-\partial W_k(x_0), -W_k(x_0))$
-   If $f_0(x_0) \ge t$, then the cut $(g, \beta)$ = $(\partial f_0(x_0), f_0(x_0) - t)$.
-   Otherwise, $x_0$ is feasible, then
    -   $t := f_0(x_0)$.
    -   The cut $(g, \beta)$ = $(\partial f_0(x_0), 0)$

### Example: Optimal matrix scalings under the min-max-ratio criterion

This example is taken from [@orlin1985computing]. According to [@orlin1985computing], optimal matrix scaling has several practical applications. One application is in linear programming, where groups of constraints and groups of variables might express the same physical commodity for which common measurement units are used. Another application is in telecommunication, where matrix scaling can be used to optimize the transmission of signals. Additionally, matrix scaling has been used in approximation theory to approximate functions of several variables by the sum of functions of fewer variables. Finally, matrix scaling has been used in Gaussian elimination, a widely used method for solving systems of linear equations, to improve the numerical stability of the algorithm.

Given a matrix $A \in \mathbb{R}^{N\times N}$. A *symmetric scaling* of $A$ is a matrix $B$ of the form $U A U^{-1}$ where $U$ is a nonnegative diagonal matrix with the same dimension. According to the *min-max criterion*, the aim is to minimize the largest absolute value of $B$'s elements [@orlin1985computing, (Program\ 3)]:
$$\begin{array}{ll}
    \text{minimize}   &  \pi  \\
    \text{subject to} &  1 \le u_i |a_{ij}| u_j^{-1} \le \Pi, \; \forall a_{ij} \neq 0 , \\
                      &  \pi, u_1 \cdot u_N \, \text{positive}. \\
  \end{array}$$

The authors show that the problems of determining the best symmetric scalings under the min-max criterion can be converted into a single parameter network optimization problem, which can be solved efficiently using the parameteric network algorithms.

Another possible criterion is to minimize the ratio of largest absolute value of the element $B$ to the smallest. One motivation for using this criterion is that high ratios cause difficulties in performing the simplex method. With this *min-max-ratio* criterion, the symmetric scaling problem can be formulated as [@orlin1985computing, (Program\ 8)]:
$$\begin{array}{ll}
    \text{minimize}   &  \pi/\psi  \\
    \text{subject to} &  \psi \le u_i |a_{ij}| u_j^{-1} \le \Pi, \; \forall a_{ij} \neq 0 , \\
                      &  \pi, \psi, u_1 \cdot u_N \, \text{positive}. \\
  \end{array}$$
Let $k’$ denotes $\log( | k | )$. By taking the logarithm of the variables, we can transform the above programming into a two-parameter network problem:
$$\begin{array}{ll}
    \text{minimize}   &  \pi’ - \psi’ \\
    \text{subject to} &  u_i’ - u_j’  \le \pi’ - a_{ij}’, \; \forall a_{ij} \neq 0 \,, \\
                      &  u_j’ - u_i’ \le a_{ij}’ - \psi’, \; \forall a_{ij} \neq 0 \,, \\
    \text{variables}  &  \pi’, \psi’, u’ \, .
  \end{array}$$
where $x = (\pi’, \psi’ )^\mathsf{T}$.
The authors of [@orlin1985computing] claim they have devised an algorithm for solving multi-parameter problems. However, we did not uncover any further publications to support this claim. Notably, the cutting-plane method readily extends the single-parameter network algorithm to be multi-parameter.

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

In convex optimization, a **Linear Matrix Inequality (LMI)** is an expression of the form:

$$A(y) = A_0 + y_1 A_1 + y_2 A_2 + \cdot + y_m A_n \succeq 0$$

where $y = [y_i, i = 1, \cdot, n]$ is a real vector, $A_0, A_1, A_2, \cdot, A_n$ are symmetric matrices, and $\succeq 0$ is a generalized inequality meaning $A(y)$ is a positive semidefinite matrix¹.

This linear matrix inequality specifies a convex constraint on $y$. There are efficient numerical methods to determine whether an LMI is feasible (e.g., whether there exists a vector $y$ such that $A(y) \succeq 0$), or to solve a convex optimization problem with LMI constraints¹.

Many optimization problems in control theory, system identification, and signal processing can be formulated using LMIs¹. Also, LMIs find application in Polynomial Sum-Of-Squares¹.

### Cholesky decomposition Algorithm

The Cholesky decomposition, also known as Cholesky factorization, is a method used in linear algebra to decompose a Hermitian, positive-definite matrix into the product of a lower triangular matrix and its conjugate transpose¹. This decomposition is useful for efficient numerical solutions, such as Monte Carlo simulations¹.

The Cholesky decomposition of a Hermitian positive-definite matrix A is a decomposition of the form A = LL*, where L is a lower triangular matrix with real and positive diagonal entries, and L* denotes the conjugate transpose of L¹. Every Hermitian positive-definite matrix (and thus also every real-valued symmetric positive-definite matrix) has a unique Cholesky decomposition¹.

If A is a real matrix (hence symmetric positive-definite), the decomposition may be written as A = LLT, where L is a real lower triangular matrix with positive diagonal entries¹.

The Cholesky decomposition is generally twice as efficient as the LU decomposition for solving systems of linear equations when it is feasible².

The **Cholesky decomposition** and the **LDLT decomposition** are both methods used in linear algebra to decompose a matrix, but they are used in different contexts and have different properties¹².

- **Cholesky decomposition**: This method decomposes a Hermitian, positive-definite matrix into the product of a lower triangular matrix and its conjugate transpose¹. It's generally faster and more numerically stable than the LDLT decomposition³. However, it requires the input matrix to be positive-definite¹.

- **LDLT decomposition**: This is a variant of the LU decomposition that is valid for positive-definite symmetric matrices². The LDLT decomposition can be applied to a broader range of matrices (we don't need a matrix to be positive-definite)¹. It decomposes a matrix into the product of a lower triangular matrix, a diagonal matrix, and the transpose of the lower triangular matrix². The LDLT decomposition is just as fast as Cholesky decomposition, but it avoids performing any square roots and is therefore faster and more numerically stable³.

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

The Cholesky or LDLT decomposition can be computed in different ways depending on the order of operations, which can be categorized into row-based and column-based methods.

- Column-Based: In this method, the computation proceeds by columns. The inner loops compute the current column using a matrix-vector product that accumulates the effects of previous columns.

- Row-Based: In this method, the computation proceeds by rows. The inner loops compute the current row by solving a triangular system involving previous rows.

Each choice of index in the outer loop yields a different Cholesky algorithm, named for the portion of the matrix updated by the basic operation in the inner loops. The choice between row-based and column-based methods depends on the specific requirements of your problem and the properties of your system (like memory layout and access patterns). With the row-based style decomposition and a lazy evaluation technique, the cutting plane can be constructed exactly in $O(p^3)$, enabling its use for efficient oracle implementations.

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

The Cholesky factorization can be used to provide a witness vector that certifies a matrix is not positive definite. If a matrix is not positive definite, the Cholesky factorization will fail¹².

During the Cholesky factorization process, there is a step where you compute the diagonal of the lower diagonal matrix as the square root of a value (let's say x). If x<0 then, this means the matrix is not positive definite¹. This failure provides evidence that the matrix is not positive definite.

When the Cholesky factorization fails due to a negative diagonal element, it means that the leading principal submatrix up to that point is not positive definite. The vector that witnesses this is one of the standard basis vectors¹². This basis vector has a 1 in the position corresponding to the failed diagonal element and zeros elsewhere. When pre-multiplied by the original matrix and post-multiplied by its transpose, it will yield a negative value, thus serving as a witness vector¹².

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

-   Perform a *row-based* Cholesky decomposition such that
    $F(x_0) = R^\mathsf{T} R$.

-   Let $A_{:p,:p}$ denotes a submatrix
    $A(1:p, 1:p) \in \mathbb{R}^{p\times p}$.

-   If Cholesky decomposition fails at row $p$,
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

Parallel Cuts {#sec:parallel_cut}
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

![Parallel cuts](ellipsoid.files/parallel_cut.pdf){width="80%"}

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

In optimization algorithms, what is particularly interesting is the convex optimization. If a problem is in a convex form, it can be efficiently and optimally solved. Convex optimization techniques are also implementable in designing FIR filters, including the Parks-McClellan algorithm [@park1972chebyshev], METEOR [@steiglitz1992meteor], and peak-constrained least-squares (PCLS) [@selesnick1996constrained; @adams1998peak]. In the mentioned articles, with the help of exchange algorithms (e.g. Remez exchange algorithm), certain FIR filter design problems can be formed as linear or quadratic programs. They are two simple forms of convex optimization problems, which can be optimally solved with existing algorithms, such as the interior-point method [@boyd2009convex]. Tempted by the optimality, more efforts were devoted to forming the problem convex. Particularly, in [@wu1999fir], via spectral decomposition [@goodman1997spectral], the problem of designing an FIR filter with magnitude constraints on frequency-domain is formulated as a convex optimization problem. More examples are provided in [@davidson2010enriching].

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

Attracted by the benefits, the authors of\ [@wu1999fir] transformed (?), originally non-convex, into a convex form via spectral decomposition:

$$L^2(\omega) \le R(\omega) \le U^2(\omega), \forall \omega\in(0,\pi)$$ {#eq:r_con}
where $R(\omega)=\sum_{i=-n+1}^{n-1}{r(t)e^{-j{\omega}t}}=|H(\omega)|^2$ and $\mathbf{r}=(r(-n+1),r(-n+2),\ldots,r(n-1))$ are the autocorrelation coefficients. Especially, $\mathbf{r}$ can be determined by $\mathbf{h}$, with the following equation vice versa\ [@wu1999fir]:

$$r(t) = \sum_{i=-n+1}^{n-1}{h(i)h(i+t)}, t\in\mathbb{Z}.$$ {#eq:h_r}
where $h(t)=0$ for $t<0$ or $t>n-1$.

![Result](ellipsoid.files/lowpass.pdf){width="80%"}

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

Discrete Optimization {#sec:discrete}
---------------------

Many engineering problems can be formulated through convex/geometric programming, such as digital circuit sizing. However, in ASIC design, there is frequently a limited number of cell types to select from in the cell library. This means that some design variables are discrete. We can map the design variables to integers to represent the discrete version as Mixed-Integer Convex programming (MICP).

What are the issues with existing methods? They are primarily based on relaxation. The more relaxed solution is used as the lower bound, and then the branch-and-bound method is applied to find the discrete optimal solution. It should be noted that the branch-and-bound method does not utilize the convexity of the issue. What if only constraints regarding discrete data could be evaluated?

Typically, a more relaxed optimal solution (convex) is obtained beforehand. After that, the optimized discrete solution is obtained using an exhaustive neighborhood search. However, tight constraints can cause a significant difference between the discrete and relaxed continuous optimal solutions. Enumerating the discrete domains can be challenging.

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

However, there are still many filter design problems that are non-convex, such as multiplierless FIR filter design problems. Note that in [@fig:fir-strctr], each coefficient associated with a multiplier unit makes the filter power-hungry, especially in *application specific integrated circuits* (ASIC). Fortunately, if each coefficient is quantized and represented as a sum of Singed Power-of-Two (SPT), a multiplierless filter can be implemented. Such coefficients can be uniquely represented by a Canonical Signed-Digit (CSD) code with a minimum number of non-zero digits[@george1960csd]. In this case, it confines the multiplication to addition and shift operations. The coefficient 0.40625 = 13/32 can be written as $2^{-1} - 2^{-3} + 2^{-5}$. Thus, the multiplier can be replaced with three shifters and two adders at a much lower cost. However, the coefficient quantization constraint is non-convex, making the convex optimization algorithm not directly applicable. A similar case is the consideration of the finite word-length effect [@lim1982finite].

Attracted by the benefits of this "multiplier-free" approach, many efforts have been devoted to its design techniques. For its general problems, integer programming (e.g. [@kodek1980design; @lim1982finite; @lim1983fir; @lim1999signed]) can be implemented to achieve the optimal solution. However, it requires excessive computational resources. Other heuristic techniques, such as genetic algorithm [@xu1995design] and dynamic-programming-like method [@chen1999trellis], also have inefficiency. If the quantization constraint is the only non-convex constraint in the design problem, a lower bound can be efficiently obtained by solving the relaxed problem [@davidson2010enriching]. Then to make the solution feasible, it can be rounded to the nearest CSD code or used as a starting point of a local search algorithm to obtain a better solution [@kodek1981comparison]. However, neither method guarantees the feasibility of the final solution. Besides, the local search problem remains non-convex. Therefore, the adopted algorithm may also be inefficient, such as branch-and-bound in [@kodek1981comparison].

![Result](ellipsoid.files/csdlowpass.pdf){width="80%"}

Concluding Remarks
==================

Should be known to students.
The ellipsoid method is not a competitor but a companion of interior-point methods.

TBD.

References {-}
==========

\ 

