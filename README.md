# Robust time-inconsistent LQ control: A stochastic differential game approach

## Usage

The Matlab code relies on the `ode45` function to solve the ODE systems in the paper.

* `compare_mu.m` generates the results when varying risk aversion, which is Figure 1 in the paper. 

* `compare_xi.m` is for varying ambiguity aversion in Figure 2.

* `compare_linear.m` considers a linear relationship between ambiguity and risk aversions. It gives Figure 3.

* The ODE system for state-dependent ambiguity aversion is solved in `sol_state.m` and `ODE_state.m`.

* The ODE system for control-dependent ambiguity aversion is solved in `sol_ctrl.m` and `ODE_ctrl.m`.

* `closedloop.m` is for the closed-loop control derived in Han, Pun, and Wong (Finance and Stochastics, 2021).
