# Scheduling-and-Routing-of-Roaming-Conductors-of-American-Railways

"Revised model.docx" is the revised model;

"Zezhou's code in Python3.ipynb" contains the work event generation, construction heuristics and metaheuristics from Zezhou's work, and I made them runable in python 3 environment

"vrp_10nodes.ipynb" and "vrp_5nodes.gms" trails solving the revised model with smaller size and in simpler situation, where all work events start at the same time and have the same time window. And both of them work out.

"vrp_20nodes.ipynb" is the trail solving the work events generated by "work event generation" of Zezhou's code. But the problem is infeasible.


"vrp_20nodes-deleting unserved nodes.ipynb" is the trail checking whether the solutions solved by heuristic method (after deleting the unserved points) satisfy the constraints of our revised model. But there are two problems: 1. Some work events are visited several times by different conductors in the solutions of heuristic method. 2. The solutions derived by heuristic method don't satisfy the constraints we constructed.
