files of type parameters_*.m create .mat file in MAT_files with model parameters.

equilibria.m calculates all possible equilibria of the model described in equation.

operating_diagram.m calls a .mat file and uses equilibria.m for different pairs (s_in,D), it registers the Operating_diagram_*.mat file that contains the different zones.

plot_operating_diagram.m uses a Operating_diagram_*.mat files and plot the different zones.

test_equilibria.m compares the steady state of a simulation with the numerical stability analysis.

d_monod: monod growth function derivative.
inv_growth_monod: 1/monod function.
dynamic: model differential equations.
dynamic_jacobian: returns the jacobian of the system at a certain point.
