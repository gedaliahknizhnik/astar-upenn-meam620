MEAM620 Project 1.3 
========================================

###### Author: Gedaliah Knizhnik

These files contain Gedaliah Knizhnik's work on project 1.3 for MEAM620 Advanced Robotics.  The goal of this phase of the project is to combine the results of projects 1.1 and 1.2,. The Dijkstra/A* implementation from 1.2 is used to find a path through a pre-defined map, which is then passed to a modified trajectory generator from 1.1, which produces a feasible trajectory that approximates the Dijkstra path. A simulated quadrotor then attempts to follow the path, using the geometric controller from 1.1, without colliding with any obstacles. The ideal trajectory is fast, efficient, feasible, and results in no collisions.

The following files contain work done by Gedaliah in addition to the code skeleton provided by the teaching staff:

1. dijkstra.m

2. trajectory_generator.m

3. crazyflie.m

4. controller.m

   

## DIJKSTRA.m

This file contains the implementation of both graph search algorithms, with the astar flag indicating whether A* or Dijkstra is used. The function takes in a map object and (x,y,z) coordinates for the start and goal positions, and returns a path, formatted as an Nx3 matrix where each row is an (x,y,z) coordinate, as well as the number of nodes expanded to find the goal.

MATLAB does not perform loops well, so in order to improve efficiency, the code is vectorized so that the only loop required is the while() loop that searches until the goal is found. 

1. Data is stored as a series of Nx1 vectors where the index corresponds to the single-value indexing of a 3D matrix.
2. The set of valid neighbors is found by exploiting how MATLAB adds a row vector to a column vector (credit to Luca Scheuer).
3. The set of nodes to explore next (i.e. the minimum cost unexplored nodes) is found using a min() function, called on a list that has infinity for all nodes but those currently open.
4. Although it is possible to return the list of all nodes with this cost, finding them takes time. I determined experimentally that it is more efficient to find the minima one at a time, since actually processing them takes next to no time and finding them this way is much faster.

A number  of helper functions are used to clean up the main loop. Detailed comments on their purpose, inputs, and outputs are located in the m-file.

Several updates were implemented to improve the function of the 1.2 Dijkstra implementation.

1. Since the voxel grid representation returns the minimum corner of the voxel, rather than its center, the output was updated so that all nodes except the start and end (which are arbitrary and not-voxelized) are moved to the center of their voxel. This gives more padding for safe trajectories.

2. I explored using the 1-norm as a heuristic for A* rather than the 2-norm. Although the 1-norm is not admissible for this 26-connected implementation, it produces much straighter paths (smoother trajectories) and improves computation time significantly.

   1. When run locally, I recommend using the 1-norm. The post-processing on the path done by the trajectory generator means the sub-optimality is largely irrelevant, and the computation time savings are highly useful during repeated testing.

   2. For submission and grading I use the 2-norm. Although it takes longer to compute a path, the paths are slightly more optimal and result in trajectories just a tad shorter.

      

## TRAJECTORY_GENERATOR.m

This file contains the trajectory generator, which takes the Dijkstra output path and transforms it into a feasible trajectory that the quadrotor can follow. This is performed as follows:

1. The path is scanned for straight line bridges between points on the original path. The farthest connectible point that doesn't result in a collision is then connected, and the intermediate original nodes removed. One intermediate point along the ray is added. This is repeated until the goal is reached, resulting in a much shorter and simpler path. This search is performed forwards and backwards, and the shorter path is used.

2. A maximum acceleration is assigned (based on experiments), and this value is used to calculate the time required for each step, assuming acceleration from 0 and deceleration to 0. That is:
   $$
   \Delta t = 2\sqrt{\frac{d}{a}}
   $$
   

   where we have used the kinematic equation 
   $$
   d = \frac{1}{2} a t^2
   $$
   and assumed that we want to be at the half-way point in distance at the half-way point in time (credit to Luca Scheuer for this method). Although the quadrotor is not accelerating  from and decelerating to zero on every segment, this approaches gives sufficient time to accelerate and decelerate overall.

3. Matrix algebra is used to calculate the coefficient for N cubic splines which will represent the N segments. The waypoints are each constrained via position, the beginning and end are constrained to have 0 velocity, and each intermediate waypoint is constrained to have continuous velocity and acceleration.

4. These coefficients are stored and used to present the trajectory when queried.

I determined experimentally that a single intermediate waypoint was sufficient. Multiple intermediate waypoints made the trajectory much slower, while no intermediate waypoints gave the cubic splines a lot of freedom to deviate from the straight line, meaning collisions would be more likely. Cubic splines were used largely for their simplicity.



## CREDITS

I originally used a set of 3D matrices rather than column vectors for my Dijkstra data structure. I shifted to a set of 1D matrices based on advice from **Mathew Halm** and **Luca Scheuer**. The method I used to calculate the set of neighbors in this representation was described to me by **Luca Scheuer**.  

The function used to return values was adapted from one written by **Jessica McWilliams**, although the particular matrix multiplication methodology used is mine. 

The acceleration based timing approach I used was shown to me by **Luca Scheuer**, and the recommendation to use continuous splines was his as well, though the implementation is my own (I had originally been using individual cubic splines).