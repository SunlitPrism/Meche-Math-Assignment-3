clear;
t = 0;
XA = 0;
h = 0.2;
rate_func_in = @rate_func01;

% G = @(X_in) XA + h * rate_func01(t+h, X_in) - X_in; 
% [XB, ~, num_evals] = multi_newton_asst3(G, XA, struct());

[XB, num_evals] = backward_euler_step(rate_func_in, t, XA, h);

h_ref = 0.1; % initial h_ref value for iterating through plots at different time step sizes
t0 = 0;
tf = 6;
X0 = 1;
tspan = [t0, tf];

[t_input_list,X_list,h_avg, tot_num_evals] = backward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref);
[t_list,x_list,h_avg, num_evals] = implicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);