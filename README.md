# 2D Lattice Boltzmann (BGK) fluid model tweaks

Originally written by Ian Haslam, commented by Greg Tucker.

I have been playing around with the puffed out Matlab code for the 2D BGK fluid model that Greg gave us. It is fun however to play with the variables to see how long it takes to run, hopefully less than 4000 time steps each time (reaching equilibrium).

I show that in the model, the minimum amount of time for the model to run falls around 0.91 for the ratio of space to obstacles (see final output of the statisitics).

In my code it is easy to change the spacing of obstacles, domain size, and flow velocities before equilibrium is reached. multiple outputs are given in .png files for various obstacle densities within the domian.

-- Trying now to figure out how to add a wall instead of random points --
