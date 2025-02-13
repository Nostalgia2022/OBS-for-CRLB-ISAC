# Optimal Beamforming Structure and Efficient Designs for CRLB Optimization in ISAC Systems
Simulation code for our paper: "Optimal Beamforming Structure and Efficient Designs for CRLB Optimization in ISAC Systems". 

# Content of Code Package
Here is a detailed description of the package:
The code in all packages are implemented in Matlab environment, and part of them is assisted by CVX toolbox.

This script provides a basic implementation of the simulation code for our paper. 

- By adjusting the weight coefficients, we obtain:
 - Figure 2: Convergence behavior
 - Figure 3: Tradeoff region

- By modifying the number of sensing streams, we obtain Figure 4.
- By changing the number of transmit and receive antennas, we obtain Figure 5 and Table 1.
- By varying the number of users, we obtain Figure 6.
- By adjusting the transmit power, we obtain Figure 7.

# Obstract of the Article
Integrated sensing and communications (ISAC) has emerged as a promising paradigm to unify wireless communications and radar sensing, enabling efficient spectrum and hardware utilization. A core challenge with realizing the gains of ISAC stems from the unique challenges of dual purpose beamforming design due to the highly non-convex nature of key performance metrics such as sum rate for communications and the Cramér–Rao lower bound (CRLB) for sensing. In this paper, we propose a low-complexity structured approach to ISAC beamforming optimization to simultaneously enhance spectral efficiency and estimation accuracy. Specifically, we develop a successive convex approximation (SCA) based algorithm which transforms the original non-convex problem into a sequence of  convex subproblems ensuring convergence to a locally optimal solution. Furthermore, leveraging the proposed SCA framework and the Lagrange duality, we derive the optimal beamforming structure for CRLB optimization in ISAC systems. Our findings characterize the reduction in radar streams one can employ without affecting performance. This enables a dimensionality reduction that enhances computational efficiency. Numerical simulations validate that our approach achieves comparable or superior performance to the considered benchmarks while requiring much lower computational costs.

# License and Referencing
This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.

# Acknowledgements
Comming soon
