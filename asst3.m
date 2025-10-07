%This function computes the value of X at the next time step
%using the Forward Euler approximation
function asst3()
    %define constants
    %t = linspace(0,5,6);
    t = 0;
    h = 0.5; %timestep
    t0 = 0;
    tf = 6;
    X0 = 1;
    tspan = [t0, tf];

    longX = [];
    longdX = [];

    for t_iter = 1:length(t)-1
        X = solution01(t_iter)
        dXdt = rate_func01(t_iter,X)
    end
    
    

    [XB,~] = forward_euler_step(@rate_func01, t, X0, h);
    [t_list, x_list, h_avg, num_evals] = forward_euler_fixed_step_integration2(@rate_func01, tspan, X0, h);
    figure()

    plot(t_list, x_list);
    X = solution01(linspace(tspan(1), tspan(2), h));
    plot(linspace(tspan(1), tspan(2), h), X);

end

