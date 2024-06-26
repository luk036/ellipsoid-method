---
标题：椭圆体法与神奇的天体 椭圆体法和神奇的天体组织
参考书目。['ellipsoid.bib'、'fir-ref.bib'、'Geostatistics.bib'、'mpcss1.bib'、'mpcss2.bib']
摘要: |
  '椭圆体法被重新审视。除此之外，还研究了三种分离神谕的应用。它们分别是鲁棒优化、半定义编程和网络优化。讨论了稳定性问题。最后，介绍了并行切割'。
---

# 介绍

椭圆法的名声很不好。与内点法相比，该方法在处理大规模凸问题时往往被认为很慢[@未知]。而这是不公平的。首先，与内点法不同，椭圆法不需要在每次迭代时对所有约束函数进行明确的评估。它所需要的只是一个提供切平面的分离神谕（@sec:cutting-plane）。这可以使该方法对于具有大量甚至无限多约束的问题具有吸引力。其次，虽然椭圆方法本身不能利用问题的稀疏性，但分离神谕可以利用问题的某些类型结构......

在@sec:robust中，鲁棒优化...

在@sec:network中，我们展示了对于网络参数问题，可以通过寻找一个有向图的负循环来获得切割平面。存在高效的算法，其中可以利用网络的位置性以及其他特性。

在@sec:lmi中，讨论了涉及矩阵不等式的问题。回顾一下，对称矩阵的正定性可以使用Cholesky，或者更准确的说是$LDL^\mathsf{T}$，因式化来有效检查。让一个对称矩阵$A \in mathbb{R}^{m x m}$。回想一下，如果因数分解过程由于遇到$A$的非正对角线条目而停在$p$行，那么$A$就不是正定性的。此外，通过使用懒惰评估技术，可以在$O(p^3)$而不是$O(m^3)$中构造一个切割平面。因此，它可以用于高效的神谕实现。

椭圆方法的实现问题在@sec:椭圆中讨论。这种方法是一种切平面方法，搜索空间是一个椭圆体，通常用以下方法表示。
$${ x \mid (x-x_c)P^{-1}(x-x_c) \le 1 \}，$$
矩阵$P \in \mathbb{R}^{n \times n}$是对称正定的。在每次迭代时，$x_c$和$P$都会根据神谕进行更新。虽然椭圆的更新很简单，而且几十年来已经找到了实现方法，但我们表明，通过将$P$拆分为$alpha$乘以$Q$，即。
$${ x \mid (x-x_c)Q^{-1}(x-x_c) \le \alpha \},$$
更新成本可以通过$n^2$ flops（浮点运算）进一步降低。

此外，在4.2节中还讨论了并行切割的使用。当一对平行不等式，其中一个不等式被违反时，可以同时使用两个约束来更新椭圆。一些论文报道，这种技术并没有提供显著的改进。然而，我们表明，对于一些约束的上界和下界很紧的情况下，例如在一些滤波器设计的情况下，使用并行切割可以极大地改善运行时间（什么的？）此外，我们还表明，如果方法实现正确，每次更新，无论是使用深切还是平行切，最多只需要一次平方根操作。

在许多实践工程问题中，一些设计变量可能被限制为离散形式。由于切割平面法只需要一个分离神谕，所以它对离散问题也能...
