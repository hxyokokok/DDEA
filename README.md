# DDEA

This code is to reproduce the experiments in the paper: X. He, Y. Zhou,  Z. Chen, and Q. Zhang, “Evolutionary Many-objective Optimization based  on Dynamical Decomposition,” IEEE Transactions on Evolutionary  Computation, pp. 1–1, 2018. 



## File structure

Three files are contained in the root directory: 

1. DDEA.m: The proposed DDEA algorithm. It relies only on the dynamic decomposition based ranking scheme (DDR) to sort the solutions.
2. DDEANS.m: The proposed DDEA+NS algorithm. It is similar to the NSGA-II except that the crowding distances are replaced with the DDR ranking results in each non-dominated front.
3. binaryTournamentSelection.m: A simple procedure for environmental selection or mating selection.



## Usage

This code is implemented on the open source software PlatEMO (http://bimk.ahu.edu.cn/index.php?s=/Index/Software/index.html).

To run the proposed two algorithms within PlatEMO, first create a new directory named DDEA in <PlatEMO_root>/Algorithms, and then put all three files in it. Launch PlatEMO, it will automatically detect the two executable algorithms and show them in the algorithm list.

For more details about how to configure the environmental settings, please see the manual of PlatEMO.
