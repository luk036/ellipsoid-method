---
title: |
  Intra-die Spatial Correlation Extraction with Maximum Likelihood
  Estimation Method
author:
  - Qiang Fu
  - Wai-Shing Luk
  - Xuan Zeng
bibliography:
  - "Geostatistics.bib"
  - "statistics.bib"
  - "ref-spatial.bib"
csl: "siam-numeric.csl"
abstract: |
  In this paper, a novel intra-die spatial correlation extraction method is
  proposed. The proposed method is based on maximum likelihood estimation of
  multivariate normal distribution for multiple samples. The obtained
  likelihood function represents the underlying statistical relationships of
  all sample chips. The proposed method can account for the nugget effect due
  to white noise, and deals with deviations among samples caused by inter-die
  variation. Experimental results have shown that the proposed method is
  efficient and practical.
keywords:
  - intra-die variations
  - spatial correlation
  - maximum likelihood estimation
---

**Key words**: intra-die variations, spatial correlation, maximum likelihood estimation.

## Introduction

As the minimum feature size of semiconductor device continues scaling
down, integrated circuits suffer from increasing variations in the
manufacturing process. These process variations lead to the geometric
variations in the devices and interconnects, and greatly affect their
electrical parameters. As a result, the performances of the fabricated
circuits are degraded from the design specifications, and the
manufacturing yield is lost. Therefore, it is more desirable to develop
statistical analysis and design methodologies to tackle variation
problems in the design stages [@Nassif00].

Process variations can be classified into two categories according to
the spatial scales. _Inter-die variations_ affect device parameters with
the same value on a die but with different values beyond the die scope,
while _intra-die variations_ cause device parameter values to vary
across different locations within a single die. As technology generation
marches, intra-die variations exceed inter-die variations and become
predominant in total process variations. Intra-die variations often show
spatial correlated patterns, which means devices that are closely placed
tend to possess similar characteristics than those further apart.

Spatial correlation is defined to describe the degree to which the
values of device characteristics are related to each other as a function
of spatial distances of the devices [@Friedberg05]. Assumed to be a
given function or matrix form, spatial correlation has been widely used
in variation aware circuit analysis and design techniques, such as
statistical timing analysis [@Chang05; @Zhang06], power/leakage
minimization [@Bhardwaj06; @Heloue07]. Recently how to model spatial
correlation from silicon measurement data has also attracted a lot of
attention. The task of spatial correlation
modeling [@Friedberg05; @Doh05; @Xiong07; @Liu07; @Hargreaves08; @Fu08]
aims to extract the characteristic parameters of spatial correlation
function provided with a large amount of silicon measurement data. It
consists of two essential issues. Firstly, an appropriate kind of
function form should be chosen to represent spatial correlation. The
_positive semidefiniteness_ property possessed by a valid spatial
correlation function should be satisfied, which means that any
correlation matrix generated from the correlation function is positive
semidefinite. Secondly, the unknown parameters in the function should be
extracted (or estimated) efficiently from measurement data, which
contain unpredictable measurement error. In [@Liu07], the exponential
function was used as spatial correlation function, and the parameter was
extracted using the curve-fitting procedure. However, how to deal with
measurement error in the extraction was not mentioned in [@Liu07].
In [@Xiong07], the Matérn function was chosen to represent spatial
correlation, and a least squares estimation (LSE) based method was
presented to extract parameters of the correlation function. Moreover,
measurement error was modeled as white noise. The method is robust
enough to recover the correlation function even if the data are
distorted by significant white noise. In [@Fu08], the unknown parameters
of Matérn function were estimated by fitting the spectral domain
counterpart of Matérn correlation function to the measurement data in
spectral domain. The method can not only handle white noise measurement
error, but also get rid of the influences of measurement error occurring
in high frequency range. In [@Hargreaves08], both exponential and Matérn
functions were chosen as spatial correlation models. The maximum
likelihood estimation (MLE) technique for a single test chip was
proposed to find the required parameters to fit spatial correlation
function to the measurement data of one chip. Therefore, a set of test
chips will result in a set of correlation functions with different
parameter values [@Hargreaves08], which are hard to use in statistical
analysis and design since only one correlation function is required. In
addition, the method did not consider the purely random component and
measurement error in measurement data, and the extraction results would
be highly distorted.

In this paper, MLE-M (MLE for multiple test chips) method is proposed to
obtain one spatial correlation function with a unique group of parameter
values from a set of test chip measurement data. The main idea is to
derive a single likelihood function for multiple chips by multiplying
the set of likelihood functions for all test chips. The unknown
parameters of the spatial correlation function are estimated from the
single likelihood function for multiple chips. Moreover, in order to
deal with the purely random component and measurement error in
measurement data which are modeled as white noise, the spatial
correlation function combined with the correlation of white noise is
used to compute the likelihood function for multiple chips. By this way,
the spatial correlation function is recovered with high accuracy.
Compared with the LSE based algorithm [@Xiong07], the proposed MLE-M
method has higher accuracy and efficiency as demonstrated by the
experimental results.

The rest of this paper is organized as follows. In Section 2, we
introduce the background about spatial correlation modeling based on
random field theory. In Section 3, the MLE-M method for intra-die
spatial correlation extraction is proposed. In Section 4, experimental
results are presented. The paper is concluded in Section 5.

## Background Material

According to [@Pitchumani05], the intra-die variation $Z$ can be further
decomposed into three components:

-1ex

- a _deterministic_ component $Z_\mathrm{det}$, which is determined by
  layout context and can be modeled by deliberatively exploring layout
  patterns;

- a _correlated random_ (or _spatially correlated_) component
  $Z_\mathrm{cor}$, which is random but shows correlated patterns due
  to proximity effects;

- a _purely random_ component $Z_\mathrm{rnd}$, which is spatially
  uncorrelated and can be treated as statistically random.

In this paper, for the sake of simplicity, the deterministic component
is assumed to be well modeled and taken away from the whole process
variation. We concentrate on the spatially correlated component,
together with the purely random component and measurement error
as [@Xiong07; @Fu08]. The spatially correlated component is generally
modeled as random field in the literature. In this section, fundamental
concepts and theories of random field are reviewed, and the Matérn
correlation function is given. The purely random component and
measurement error cause a discontinuity at the origin of the correlation
function, which is called nugget effect. We will describe this
phenomenon and give the modification form of the Matérn correlation
function in this section.

### Random Field [@Schabenberger05]

_Random field_, also known as _stochastic process_, can be regarded as
an indexed family of random variables denoted as
{$Z(\mathbf{s}): \mathbf{s}\in D$}, where $D$ is a subset of
$d$-dimensional Euclidean space $\mathbb{R}^d$. To specify a stochastic
process, the joint probability distribution function of any finite
subset $(Z(\mathbf{s}_1), \ldots, Z(\mathbf{s}_n))$ must be given in a
consistent way, which is called _distribution_ of the process. For ease
of analysis, a random field is often assumed to be with _Gaussian_
distribution, and is called Gaussian random field.

A random field has several key properties that are useful in practical
problems. The field is _stationary_ under translations, or
_homogeneous_, if the distribution is unchanged when the point set is
translated. The field is _isotropic_ if the distribution is invariant
under any rotation of the whole points in the parameter space. We study
homogeneous isotropic field in this paper.

The _covariance_ $C$ and _correlation_ $R$ of a stochastic process are
defined by
$$C(\mathbf{s}_i,\mathbf{s}_j) = \mathrm{cov}(Z(\mathbf{s}_i),Z(\mathbf{s}_j)) = \mathrm{E}\lbrack (Z(\mathbf{s}_i)-\mathrm{E}\lbrack Z(\mathbf{s}_i)\rbrack)(Z(\mathbf{s}_j)-\mathrm{E}\lbrack Z(\mathbf{s}_j)\rbrack)\rbrack \nonumber$$
and
$$R(\mathbf{s}_i,\mathbf{s}_j)=C(\mathbf{s}_i,\mathbf{s}_j)/ \sqrt{C(\mathbf{s}_i,\mathbf{s}_i)C(\mathbf{s}_j,\mathbf{s}_j)}$$
respectively for all $\mathbf{s}_i,\mathbf{s}_j\in D$, where
$\mathrm{E}\lbrack Z(\mathbf{s})\rbrack$ denotes the expectation of
$Z(\mathbf{s})$. Thus a process is homogeneous if $C$ and $R$ depend
only on the separation vector $\mathbf{h}=\mathbf{s}_i-\mathbf{s}_j$.
Furthermore, it is isotropic if $C$ and $R$ depend upon $\mathbf{h}$
only through its length $h$, i.e.,
$$C(\mathbf{s}_i,\mathbf{s}_j)=C(\mathbf{h})=C(h),$$

$$R(\mathbf{s}_i,\mathbf{s}_j)=R(\mathbf{h})=R(h)=C(h)/C(0).$$ {#eq:corr_def}
If we denote $C(0)$, the variance of $Z(\mathbf{s})$, as $\sigma^2$, then the
relationship between covariance and correlation is $C(h)=\sigma^2 R(h)$.

### Correlation Function for The Spatially Correlated Component

-----------------------------------------------------------

The spatially correlated component is modeled as random field with
variance $\sigma^2$ and correlation function $R(h)$ to be extracted.
Positive semidefiniteness is the necessary condition for a parametric family
of functions to define a legitimate class of correlation functions.
However, this is not an easy condition to check directly. In practice,
this is usually ensured by working within one of several standard
classes of parametric model for $R(h)$. Although no certain criteria
exist for spatial correlation representation, we prefer to select the
Matérn function as it is more flexible than other correlation functions.
However, the proposed extraction methods are also adaptive to other
correlation functions. The Matérn function for a homogeneous isotropic
field is in the form as
$$R(h)=\frac{\displaystyle 1}{\displaystyle 2^{\nu-1}\Gamma(\nu)}(\alpha h)^\nu K_\nu(\alpha h)$$
where $K_\nu$ is the second kind modified Bessel function of order
$\nu$, and $\Gamma$ is the gamma function.

In Matérn function, $\alpha$ is the range parameter which measures how
quickly correlation decays with distance, and $\nu$ is the shape
parameter which controls smoothness of the process. When $\nu$ is small,
it implies the field is rough. When $\nu$ is large, it means the field
is smooth. Matérn function is a generalization of several other widely
used correlation functions. For $\nu=1/2$, Matérn function reduces to
exponential function; when $\nu$ tends to infinity, it degenerates as
Gaussian function.
@fig:Matern_funs shows Matérn functions with different
parameter groups.

::: {#fig:Matern_funs}
Matérn functions with different parameter groups.
:::

### Correlation Function Considering Nugget Effect

-----------------------------------------------------------

In the measurement data, apart from the spatially correlated component,
the purely random component and the unavoidable measurement error also
exist. Typically, the two components are modeled as random variables
with independent identical Gaussian distribution, or in short, Gaussian
white noise. The correlation function of white noise and the
corresponding spectral representation (named as "spectral
density" [@Fu08]) are shown in
@fig:white_noise.

::: {#fig:white_noise}
Correlation function and spectral density of white noise.
:::

When the two components are considered, the measurement data can still
be regarded as a Gaussian random field, but the correlation function
will have a discontinuity at the origin. This phenomenon is called
"nugget effect\" [@Diggle07]. Denoting the sum variance of the above two
components as $\tau^2$, the original Matérn correlation function $R(h)$
should be modified as $$\tilde{R}(h)=\left\{ \begin{array}{l l}
 1 & \textrm{if $h=0$} \\
 \frac{\displaystyle \sigma^2}{\displaystyle \sigma^2+\displaystyle \tau^2}\cdot R(h) = \frac{\displaystyle \sigma^2}{\displaystyle \sigma^2+\displaystyle \tau^2}\cdot \frac{\displaystyle 1}{\displaystyle 2^{\nu-1}\Gamma(\nu)}(\alpha h)^\nu K_\nu(\alpha h) & \textrm{if $h>0$}
\end{array} \right.
$$

@fig:modi_Matern_funs shows the modified Matérn correlation
functions considering the nugget effect for the original Matérn
correlation functions in
@fig:Matern_funs.

::: {#fig:modi_Matern_funs}
Modified Matérn functions with different parameter groups.
:::

In [@Hargreaves08], $R(h)$ was simply used in the extraction. In this
paper, the proposed method uses the modified form $\tilde{R}(h)$ to
account for the nugget effect and the extraction results are more
accurate.

## Intra-die Spatial Correlation Extraction Based on MLE Method

In this section, we will first formulate the problem of spatial
correlation function extraction. We further deduce the likelihood
function for multiple samples, and present the extraction method in
details.

### Problem Formulation

In the measurement process, we sample to gather measurement data over a
batch of $M$ chips, with each chip comprising $N$ measurement sites. The
spatial locations of the measurement sites on each chip should be
identical.

Assume that the spatially correlated component of intra-die variation is
modeled as a Gaussian random field with Matérn correlation function
$R(\vec{\psi};h)$ (denoted as $R(\vec{\psi})$ for short). The problem of
spatial correlation function extraction is defined to be: given
measurement data representing process variations, extract the unknown
parameter vector $\vec{\psi}=(\alpha,\nu)$ that completely specifies the
correlation function, and the variance of the spatially correlated
component $\sigma^2$, while the interference of nugget effect and
inter-die variation should be considered. The objective is to recover
the spatial correlation function $R(\vec{\psi})$ as accurately as
possible after plugging the estimated $\vec{\psi}$ into Matérn function.

### Spatial Correlation Extraction Based on MLE for Multiple Samples

For a sample $\vec{z}=(z(\mathbf{s}_1),\ldots,z(\mathbf{s}_N))^T$ from
an $N$-variate Gaussian random field with zero mean vector and
$N \times N$ covariance matrix $C=\sigma^2 R(\vec{\psi})$, the
distribution is
$$p(\vec{z}|\sigma^2,\vec{\psi})=\frac{1}{(2\pi)^{N/2}(\mathrm{det}(\sigma^2 R))^{1/2}} \cdot \mathrm{exp}\left(-\frac{1}{2\sigma^2}\vec{z}^T R^{-1}\vec{z} \right)$$
_Likelihood function_ is defined by reversing the roles of the data
vector and the parameter vector in the distribution, as
$L(\sigma^2,\vec{\psi}|\vec{z})$, to represent the likelihood of the
parameters given the observed data, and $L(\sigma^2,\vec{\psi})$ is
often used with the understanding that the observed data $\vec{z}$ is
implicit. In practice, for computational convenience the logarithm of
the likelihood is taken and estimates are obtained by maximizing the
log-likelihood function. The log-likelihood function of a sample is

$$
\log L(\sigma^2,\vec{\psi}) = -\frac{N}{2}\log 2\pi-\frac{N}{2}\log \sigma^2-\frac{1}{2}\log \mathrm{det}~R-\frac{1}{2\sigma^2}\vec{z}^T R^{-1}\vec{z}
$$ {#eq:likhood_uni}

As there are $M$ samples obtained in the measurement process,
in [@Hargreaves08],
@eq:likhood_uni is used $M$ times to estimate the unknown
parameters of each sample separately, resulting in $M$ groups of
different parameter values unrelated to each other. However, designers
do not simply care for the parameters of each individual chip. Except
for the particular characteristic of each specific chip in these
samples, there is no confidence for which group of parameter values to
choose for inference. The method in [@Hargreaves08] provides no guidance
for robust circuit design, where we actually need the unique group of
parameters that capture the characteristics buried in all samples.

In fact since we have $M$ samples of observations
$\vec{z}_1,\ldots,\vec{z}_M$, with the $m$-th sample containing $N$
variates as $\vec{z}_m=(z_m(\mathbf{s}_1),\ldots,z_m(\mathbf{s}_N))^T$,
it is more natural to calculate the likelihood based upon all samples,
i.e., we want to search for the parameters that maximize the likelihood
of all the samples. Before we do that, two processing steps should be
performed to coincide with the actual extraction process. Firstly, as
the $M$ samples of observations are from $M$ different chips, they are
inevitably affected by inter-die variations. Although the percentage of
inter-die variations in total process variations becomes comparatively
smaller in the sub-micron regime, the extraction results will be
distorted if the influence is neglected. To alleviate the sample
deviations caused by inter-die variations, we remove the mean effect
through a simple averaging process as
$$z_m^*(\mathbf{s}_i)=z_m(\mathbf{s}_i) - \frac{1}{M}\sum_{k=1}^M z_k(\mathbf{s}_i)$$
and the $m$-th sample of data vector becomes
$\vec{z}_m^*=(z_m^*(\mathbf{s}_1),\ldots, z_m^*(\mathbf{s}_N))^T$.
Secondly, as mentioned in Sect. 2.3, the purely random component and
measurement error result in the nugget effect, which should be taken
into consideration in the extraction. Replacing $R(h)$ with
$\tilde{R}(h)$, the modified covariance matrix of measurement data
becomes $$\tilde{C}=\sigma^2 R(\vec{\psi})+\tau^2 I$$ where
$\sigma^2 R(\vec{\psi})$ is the covariance of the spatially correlated
component, while $\tau^2 I$ is the covariance caused by nugget effect.
Parameterizing $$\kappa=\tau^2/\sigma^2,$$ the modified correlation
matrix of measurement data is
$$\tilde{R}(\kappa,\vec{\psi})=R(\vec{\psi})+\kappa I.$$ The $N$-variate
Gaussian measurement process can be denoted as
$Z\sim \mathrm{N}(\mathbf{0},\sigma^2 \tilde{R}(\kappa,\vec{\psi}))$.

After the processings, the likelihood function for all the $M$ samples
is [@Anderson03]
$$L(\sigma^2,\kappa,\vec{\psi}) = \prod_{m=1}^{M}L_m(\sigma^2,\kappa,\vec{\psi})=\frac{1}{(2\pi)^{MN/2}(\mathrm{det}(\sigma^2 \tilde{R}))^{M/2}} \cdot \mathrm{exp}\left(-\frac{1}{2\sigma^2}\sum_{m=1}^{M}\vec{z}_m^{*T} \tilde{R}^{-1}\vec{z}_m^{*} \right)$$
and the log-likelihood function is
$$\log L(\sigma^2,\kappa,\vec{\psi}) = \log \left(\prod_{m=1}^{M}L_m(\sigma^2,\kappa,\vec{\psi})\right)=-\frac{MN}{2}\log 2\pi-\frac{MN}{2}\log\sigma^2-\frac{M}{2}\log \mathrm{det}~\tilde{R}-\frac{1}{2\sigma^2}\sum_{m=1}^{M}\vec{z}_m^{*T} \tilde{R}^{-1}\vec{z}_m^*$$ {#eq:loglik}
By setting $\frac{\partial \log L}{\partial \sigma^2}=0$, we can get the
estimation of $\sigma^2$ as
$$\hat{\sigma}^2(\tilde{R})=\frac{1}{MN}\sum_{m=1}^{M}\vec{z}_m^{*T} \tilde{R}^{-1}\vec{z}_m^*$$ {#eq:sigma2_est}
Substituting
@eq:sigma2_est back into
@eq:loglik and ignoring the constants, we can obtain the
concentrated log-likelihood function
$$\log L_0(\kappa,\vec{\psi})=-\log \mathrm{det}~\tilde{R} - N\log \hat{\sigma}^2(\tilde{R})$$ {#eq:con_loglik}

Regard the scalar $\vec{z}^{*T} \tilde{R}^{-1}\vec{z}^*$ as the trace of
a $1\times 1$ matrix. Using the property of the trace of a matrix that
$\mathrm{tr}(AB)=\mathrm{tr}(BA)$ whenever $A$ and $B$ are matrices so
shaped that both products exist, we can reformulate $\hat{\sigma}^2$ as

$$
\hat{\sigma}^2 = \frac{1}{MN}\sum_{m=1}^{M}\mathrm{tr}(\vec{z}_m^{*T} \tilde{R}^{-1}\vec{z}_m^*)=\frac{1}{MN}\sum_{m=1}^{M}\mathrm{tr}(\vec{z}_m^*\vec{z}_m^{*T} \tilde{R}^{-1})=\frac{1}{N}~\mathrm{tr}\left(\frac{1}{M}\sum_{m=1}^{M}\vec{z}_m^*\vec{z}_m^{*T} \tilde{R}^{-1}\right)=\frac{1}{N}~\mathrm{tr}(Y\tilde{R}^{-1})
$$ {#eq:sigma2_est2}

where
$$Y=\frac{1}{M}\sum_{m=1}^{M}\vec{z}_m^*\vec{z}_m^{*T}$$
Substituting
@eq:sigma2_est2 back into
@eq:con_loglik and ignoring the constant, we obtain the
log-likelihood function of MLE for multiple samples as
$$\log L_0(\kappa,\vec{\psi})=-\log \mathrm{det}~\tilde{R} - N\log (\mathrm{tr}(Y\tilde{R}^{-1}) )$$
This must be optimized numerically with respect to $\kappa$ and
$\vec{\psi}$, followed by back substitution to obtain $\hat{\sigma}^2$.

In the numerical evaluation steps of log-likelihood function, we find
that direct computation of the log determinant of $\tilde{R}$ often
incurs "log of zero\" error, as the determinant of $\tilde{R}$ may
approach to zero too closely. However, if there is no numerical errors,
a positive value close to zero will result in a negative value in the
acceptable range after taking the logarithm. To overcome this problem,
we first decompose $\tilde{R}$ into two matrices $L$ and $U$ using LU
factorization, where $U$ is an upper triangular matrix, and $L$ is the
permutation of a lower triangular matrix whose diagonal entries are all
ones. Since $\tilde{R}$ should be positive definite,
$\mathrm{det}~\tilde{R}>0$. Using the property of the determinant that
$\mathrm{det}~\tilde{R}=|\mathrm{det}~L| \cdot |\mathrm{det}~U|$, and
$|\mathrm{det}~L|=1$, the log determinant of $\tilde{R}$ can be obtained
by summing up the absolute values of the diagonal elements of the $U$
matrix, i.e.,
$$\log \mathrm{det}~\tilde{R}=\log |\mathrm{det}~U|=\log |\prod_{j=1}^N u_{jj}|=\sum_{j=1}^N \log |u_{jj}|$$
where $u_{mm}$ is the diagonal element of $U$ and $|\bullet|$ means
taking the absolute value. The finally obtained log-likelihood function
is in the form as
$$\log L_0(\kappa,\vec{\psi})=-\sum_{j=1}^N \log |u_{jj}| - N\log (\mathrm{tr}(Y\tilde{R}^{-1}))$$
This can be solved by any standard nonlinear optimization technique. In
our implementation, we use the _fmincon_ function in MATLAB which is
based on a sequential quadratic programming method.

## Experimental Results

The proposed method was implemented in MATLAB on an Intel machine with
3.0 GHz XEON CPU. Without real silicon measurement data, we synthesized
a pseudo measurement process and attained data from a batch of $M$ chips
with $N$ sites on each chip. The area of the chip is set to be
$10\,\mathrm{mm} \times 10\,\mathrm{mm}$. We conducted two groups of experiments according to
different sampling schemes. In the first group, the $N$ sites are evenly
distributed on the chip in horizontal and vertical directions
respectively, which is called "uniform gridding sampling scheme". In the
second group, the $N$ sites reside on the two-dimensional plane in a
completely random way, which is called "Monte Carlo sampling scheme".
The synthetic data consist of four components, i.e., the spatially
correlated component of intra-die variation, the purely random component
of intra-die variation, inter-die variation, and measurement error. The
spatially correlated component is the most important part, which is
generated as a homogeneous isotropic field with Gaussian distribution
and Matérn correlation function. The "Cholesky decomposition" method
described in [@Cressie91] and used in [@Hargreaves08; @Fu08] was also
utilized for the generation.
@fig:fig2_a shows the surface plot of a generated random
field. The purely random component and measurement error were added as
Gaussian white noise, whose variances were chosen with different
percentages (labeled as $\kappa$) w.r.t. the variance of the spatially
correlated component. For inter-die variation, a global variation
component was added to the characteristic values of all sites on a same
chip. Different chips were added with different global variation values
from a Gaussian distribution. The above three variation components were
generated with a set of different variance values.
@fig:fig2_b shows the surface plot of the finally generated
measurement data from a sample chip.

::: {#fig:fig2_a}
Random field.
:::

::: {#fig:fig2_b}
Measurement data.
:::

We implemented the proposed method in two branches for comparison. One
does not consider nugget effect and use the original $R(h)$ in
likelihood function for estimation which is called MLEsim, while the
other accounts for nugget effect and use the modified form
$\tilde{R}(h)$ in likelihood function which is called MLEnug. We also
implemented the LSE based extraction algorithm in [@Xiong07] named as
RESCF. We performed multi-runs (30 times) in each experiment and present
the average results in
Tables~\ref{tab:1_uni_grid}
and~\ref{tab:1_MC}.

```{=latex}
\begin{table*}[tb]
  \caption{Results of spatial correlation extraction with uniform gridding sampling scheme}\label{tab:1_uni_grid}
  \centering
  \resizebox{\linewidth}{!}{%
      \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|} \hline
               & & & \multicolumn{3}{|c|}{RESCF} & \multicolumn{3}{|c|}{MLEsim} & \multicolumn{3}{|c|}{MLEnug} \\ \cline{4-12}
         $M$  & $N$ & $\kappa$ & err($\sigma^2$) & err($R(h)$) & $t$(sec.) & err($\sigma^2$) & err($R(h)$) & $t$(sec.) & err($\sigma^2$) & err($R(h)$) & $t$(sec.) \\ \hline
        500 & 11$\times$11 & 10\% & 2.14\% & 1.69\% & 10.43 & 10.98\% & 5.50\% & 2.73 & 0.95\% & 0.58\% & 7.63 \\ \cline{3-12}
              & & 50\% & 6.05\% & 3.81\% & 12.33 & 50.21\% & 23.83\% & 2.04 & 1.95\% & 1.17\% & 5.48 \\ \cline{3-12}
              & & 100\% & 7.31\% & 4.26\% & 10.31 & 98.76\% & 38.69\% & 2.80 & 4.12\% & 2.27\% & 5.89 \\ \cline{3-12} \hline
        1000 & 11$\times$11 & 10\% & 2.45\% & 1.58\% & 9.62 & 11.66\% & 5.44\% & 2.60 & 0.84\% & 0.45\% & 7.85 \\ \cline{3-12}
              & & 50\% & 4.94\% & 3.02\% & 15.23 & 50.64\% & 23.71\% & 2.14 & 2.17\% & 1.02\% & 5.91 \\ \cline{3-12}
              & & 100\% & 6.58\% & 3.90\% & 12.24 & 99.18\% & 38.57\% & 2.78 & 3.18\% & 1.75\% & 5.97 \\ \cline{3-12} \hline
        500 & 21$\times$21 & 10\% & 1.19\% & 1.09\% & 137.77 & 13.16\% & 13.90\% & 30.72 & 0.38\% & 0.30\% & 131.28 \\ \cline{3-12}
              & & 50\% & 1.65\% & 1.66\% & 231.66 & 25.93\% & 42.41\% & 34.15 & 0.99\% & 0.63\% & 87.84 \\ \cline{3-12}
              & & 100\% & 2.75\% & 1.98\% & 227.44 & 79.64\% & 55.72\% & 54.56 & 1.79\% & 0.77\% & 70.38 \\ \cline{3-12} \hline
        1000 & 21$\times$21 & 10\% & 0.95\% & 0.75\% & 185.31 & 12.90\% & 13.88\% & 30.72 & 0.30\% & 0.27\% & 132.76 \\ \cline{3-12}
              & & 50\% & 1.69\% & 1.32\% & 238.66 & 26.09\% & 42.42\% & 34.82 & 1.07\% & 0.48\% & 90.48 \\ \cline{3-12}
              & & 100\% & 2.46\% & 1.60\% & 222.42 & 80.04\% & 55.66\% & 55.32 & 1.37\% & 0.68\% & 68.15 \\ \cline{3-12} \hline
      \end{tabular}%
  }
\end{table*}
```

```{=latex}
\begin{table*}[tb]
  \caption{Results of spatial correlation extraction with Monte Carlo sampling scheme}\label{tab:1_MC}
  \centering
  \resizebox{\linewidth}{!}{%
      \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|} \hline
               & & & \multicolumn{3}{|c|}{RESCF} & \multicolumn{3}{|c|}{MLEsim} & \multicolumn{3}{|c|}{MLEnug} \\ \cline{4-12}
         $M$  & $N$ & $\kappa$ & err($\sigma^2$) & err($R(h)$) & $t$(sec.) & err($\sigma^2$) & err($R(h)$) & $t$(sec.) & err($\sigma^2$) & err($R(h)$) & $t$(sec.) \\ \hline
        500 & 11$\times$11 & 10\% & 1.85\% & 2.20\% & 15.93 & 25.53\% & 45.34\% & 2.27 & 0.63\% & 0.70\% & 7.55 \\ \cline{3-12}
              & & 50\% & 2.01\% & 2.33\% & 16.65 & 27.80\% & 63.16\% & 2.57 & 1.59\% & 1.01\% & 5.14 \\ \cline{3-12}
              & & 100\% & 2.32\% & 2.18\% & 14.98 & 84.23\% & 71.75\% & 2.50 & 1.67\% & 1.39\% & 5.00 \\ \cline{3-12} \hline
        1000 & 11$\times$11 & 10\% & 1.63\% & 1.65\% & 15.04 & 23.69\% & 42.48\% & 2.17 & 0.35\% & 0.45\% & 7.41 \\ \cline{3-12}
              & & 50\% & 1.57\% & 1.30\% & 15.50 & 27.48\% & 63.49\% & 2.50 & 1.29\% & 0.91\% & 5.19 \\ \cline{3-12}
              & & 100\% & 1.93\% & 1.52\% & 15.03 & 86.13\% & 72.93\% & 2.19 & 1.49\% & 1.15\% & 4.93 \\ \cline{3-12} \hline
        500 & 21$\times$21 & 10\% & 1.83\% & 1.27\% & 200.43 & 47.15\% & 56.20\% & 32.88 & 0.39\% & 0.32\% & 109.78 \\ \cline{3-12}
              & & 50\% & 1.71\% & 1.33\% & 217.19 & 7.58\% & 72.74\% & 32.81 & 1.34\% & 0.60\% & 84.32 \\ \cline{3-12}
              & & 100\% & 1.92\% & 1.36\% & 177.16 & 70.16\% & 77.92\% & 25.16 & 1.32\% & 0.67\% & 79.22 \\ \cline{3-12} \hline
        1000 & 21$\times$21 & 10\% & 1.33\% & 1.20\% & 209.48 & 46.04\% & 53.49\% & 34.27 & 0.32\% & 0.29\% & 107.84 \\ \cline{3-12}
              & & 50\% & 1.12\% & 1.06\% & 223.88 & 8.70\% & 71.51\% & 33.83 & 0.80\% & 0.54\% & 84.05 \\ \cline{3-12}
              & & 100\% & 0.91\% & 1.14\% & 204.29 & 69.44\% & 77.61\% & 27.86 & 0.81\% & 0.47\% & 79.96 \\ \cline{3-12} \hline
      \end{tabular}%
  }
\end{table*}
```

In the experiments, the relative error of the variance of the spatially
correlated component was given by
$$\mathrm{err}(\sigma^2)=|\frac{\hat{\sigma}^2-\sigma^2}{\sigma^2}|,$$
and the relative error of the spatial correlation function was computed
by $$\mathrm{err}(R(h))=\frac{||\hat{R}(h)-R(h)||}{||R(h)||},$$ in which
$\sigma^2$ and $R(h)$ denote the actual variance and spatial correlation
function used to generate the random field, while $\hat{\sigma}^2$ and
$\hat{R}(h)$ refer to the extracted variance and spatial correlation
function after plugging in the estimated parameters.

From the two tables, we can see that when nugget effect is ignored in
the extraction, MLEsim works badly. For most of the test cases, the
relative errors of the extracted spatial correlation function are above
30%, or even higher to 70% which are totally unacceptable. When nugget
effect is considered, MLEnug works well, and achieves better results
than RESCF for all test cases with less runtime for both two sampling
schemes. This shows the accuracy and efficiency of the proposed method.
@fig:corr_funs illustrates the correlation functions
extracted by the three methods compared with the actual correlation
function in one experiment.

::: {#fig:corr_funs}
Comparison of correlation function extraction.
:::

## Conclusion Remarks

Intra-die spatial correlation extraction has been a prevailing subject
these years. In this paper, we propose a novel extraction method using
MLE technique. The proposed method is based on MLE of multivariate
normal distribution for multiple samples. After removing the sample
deviations caused by inter-die variations, the proposed method takes all
sample chips as a whole to formulate the likelihood function, in which
nugget effect is integrated to account for the purely random component
of intra-die variation with measurement error. Experimental results show
that the proposed method outperforms the LSE based extraction
algorithm [@Xiong07] for all test cases with less runtime. The proposed
method is more efficient and stable for spatial correlation extraction
provided with real silicon data.

## Acknowledgment {#acknowledgment .unnumbered .unnumbered}

This research is supported partly by NSFC research project

## References {-}

\


