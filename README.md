# Comp_phys
Some computational physics projects/code.

* Active homogeneous simulates an active fluid in shear flow. Newton's Method is used to solve for the dynamics of a nematic order parameter, and stress vs strain curves are generated. Requires scipy, matplotlib, and numpy.

* Crank solves for the inhomogeneous case, using a Crank-Nicholson finite difference method, and assuming grad(Q) = 0 at the boundaries. Requires numpy and matplotlib.

* Scalar wave solves the scalar wave equation in box, using spectral methods, and generates a animation of this. Requires fftw.

* Finite methods is a python program that uses finite element analysis to determine the final temperature on a bar
attached to a wall, given an initial temperature, a function defining the bar over its length, and other parameters. Requires scipy and numpy.


* Finance predicts the future of a particular stock via Monte Carlo methods, given its volatility, strike price, and other parameters. Requires matplotlib.

* Thermrand creates a random grid of cells with some integer from 1-10, and for each step in time, the highest number moves in a random direction into another cell. The grid is eventually plotted with red representing a higher sum of numbers in the cell, blue being a lower sum. Requires Tkinter.

