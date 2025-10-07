%Runs numerical integration using forward Euler approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration2(rate_func_in,tspan,X0,h_ref)
    
   
    t0 = tspan(1);
    tf = tspan(2);

    N = ceil((tf-t0) / h_ref);
    h = (tf-t0) / N;

    global t_input_list;
    t_input_list = linspace(t0,tf,N+1);
    global X_list;
    X_list = [X0];
    tot_num_evals = 0;
   
    for i  = 1:N
        [XB, num_evals] = forward_euler_step(rate_func_in,t_input_list(i),X0,h);
        X_list(end+1) = XB;
        tot_num_evals = tot_num_evals + num_evals;
        X0 = XB;

    end
    figure;
    plot(t_input_list,X_list)

end