# CP-KERN

This repository contains the code implementing the CP-KERN algorithm described in the following paper:
- Singh, Abhishek. (2023). Cutting-plane algorithms for preemptive uniprocessor scheduling problems. Real-Time Systems. 10. DOI: [1007/s11241-023-09408-y](https://doi.org/10.1007/s11241-023-09408-y). 

CP-KERN checks the unschedulability condition for FP-systems. The counterpart for EDF-systems is still under development.

### Input file
The input file has the format: 
- [number of tasks]

and for each task: 
- [period] [deadline] [jitter] [wcet]


## Assumptions

### Integer values
$C_i$ and $T_i$ are strictly positive integer values.


### Utilization 
The sum of the utilization $U_j$ of all tasks is less than 1, that is : 
 $\sum_{j\in [n]} U_j <= 1$



#### Example systems

##### Example 1 (from page 28 of Singh's paper)

The following file describes a system that doesn't pass the schedulability test, because of an unsatisfied deadline at t= 10: 

```
3
1e-100
17 10 0 6
13 10 0 5
20 31 0 1
```

$\hat{D}_min = 10$ and $L = 13$, while $p=2$ and $q=3$. The solution is $t=-10$.


#### Example 2 (from page 12 of Singh's paper)

An FP task system with implicit deadlines and zero release jitter:

```
3
1e-100
40 40 0 20
50 50 0 10
150 150 0 33
```

The system is unschedulable and CP-KERN detects a missed deadline at $t=143$ after 3 iterations.


