ell-review.pdf

# Abstract:
The main topic of this document is the ellipsoid method and its application in various optimization problems. The document begins by addressing the reputation of the ellipsoid method, which is commonly believed to be slow compared to interior point methods. However, the author argues that the ellipsoid method has advantages, such as not requiring the evaluation of all constraint functions and being able to exploit certain types of problem structure.

The document then delves into three specific applications of the ellipsoid method: robust optimization, semidefinite programming, and network optimization. In each case, the author discusses the use of separation oracles and their efficiency in solving these types of problems.

The implementation issues of the ellipsoid method are also explored, including the use of parallel cuts to update the search space. The author presents evidence that, in certain cases, parallel cuts can significantly improve runtime.

Additionally, the document mentions that the cutting-plane method can be applied to discrete problems as well, where design variables are restricted to discrete forms.

Overall, this document presents a comprehensive exploration of the ellipsoid method and its applications in optimization problems, highlighting its advantages and discussing implementation issues and techniques.


# How can multiplier-less FIR filter design problems be implemented?

Multiplier-less FIR filter design problems can be implemented by quantizing the filter coefficients and representing them as a sum of Signed Power-of-Two (SPT) terms. Each coefficient can be uniquely represented by a Canonic Signed-Digit (CSD) code with the smallest number of non-zero digits. This allows the multiplications to be replaced with add and shift operations, resulting in lower cost implementation. The coefficient quantization constraint, which is non-convex, makes it challenging to directly apply convex optimization algorithms. Integer programming techniques can be used to solve these problems, such as applying a trellis search algorithm or using a parallel genetic algorithm .
