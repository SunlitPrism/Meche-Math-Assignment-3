%This function computes the value of X at the next time step
%using the Forward Euler approximation
function asst3()
    %define constants
    t = linspace(0,5,6);
    %t = 0;
    t_ref = 0.492;
    h_ref = 0.2; %reference timestep
    t0 = 0;
    tf = 6;
    X0 = 1;
    tspan = [t0, tf];
    hspan = [1e-8, 1e-1, 15];

    for t_iter = 1:length(t)-1
        X = solution01(t_iter);
        dXdt = rate_func01(t_iter,X);
    end
    
    %[XB, ~] = forward_euler_step(@rate_func01, t, X0, h);
    [t_list, x_list, h_avg, num_evals] = forward_euler_fixed_step_integration2(@rate_func01, tspan, X0, h_ref);
    [t_mid_list,x_list_expmid,h_avg_expmid, num_evals_expmid] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);  
    
    figure()
    hold on;
    plot(t_list, x_list,'gO'); % Plot Euler step integration over time
    plot(t_mid_list, x_list_expmid,'rO'); % Plot Explict Midpoint step integration over time

    X = solution01(linspace(tspan(1), tspan(2), (tspan(2)/h_avg))); %Calc Analytical
    plot(linspace(tspan(1), tspan(2), (tspan(2)/h_avg)), X,'b'); % Plot analytical solution over time
    legend("Forward Euler","Explicit Mid","Analytical")
    
    % [error_list] = truncation_error(t_ref, hspan);
    % plot(t_list, error_list);
end

