## Minimum cycle ratio problem

The minimum cycle ratio (MCR) problem is a fundamental problem in the analysis of directed graphs. Given a directed graph, the MCR problem seeks to find the cycle with the minimum ratio of the sum of edge weights to the number of edges in the cycle. In other words, the MCR problem seeks to find the "tightest" cycle in the graph, where the tightness of a cycle is measured by the ratio of the total weight of the cycle to its length.

The MCR problem has many applications in the analysis of discrete event systems, such as digital circuits and communication networks. It is closely related to other problems in graph theory, such as the shortest path problem and the maximum flow problem. Efficient algorithms for solving the MCR problem are therefore of great practical importance.


## Single Negative Cycle

According to the PDF file, there are several negative cycle detection algorithms that are considered to be among the fastest in practice. These include Tarjan's and Szymanski's negative cycle detection algorithms, as well as the Goldberg-Radzik algorithm. The choice of algorithm may depend on the specific application and the requirements for running time, space complexity, and other factors.

## Howard's method

Howard's method is a minimum cycle ratio (MCR) algorithm that uses a policy iteration algorithm to find the minimum cycle ratio of a directed graph. The algorithm maintains a set of candidate cycles and iteratively updates the cycle with the minimum ratio until convergence. 

To detect negative cycles, Howard's method uses a cycle detection algorithm that is based on the Bellman-Ford algorithm. Specifically, the algorithm maintains a predecessor graph of the original graph and performs cycle detection on this graph using the Bellman-Ford algorithm. If a negative cycle is detected, the algorithm terminates and returns the cycle.
