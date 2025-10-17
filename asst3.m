%This function computes the value of X at the next time step
%using the Forward Euler approximation
function asst3()
    close all
    %define constants
    t = linspace(0,5,6);
    %t = 0;
    t_ref = 0.1;
    %h_ref = 0.2; % Initially chosen h_ref value
    %h_ref = 0.45; % Need plots for h_ref = 0.38 and 0.45
    h_ref = 0.1; % initial h_ref value for iterating through plots at different time step sizes
    t0 = 0;
    tf = 6;
    X0 = 1;
    tspan = [t0, tf];
    hspan = [-5, -1, 100];

    for t_iter = 1:length(t)-1
        X = solution01(t_iter);
        dXdt = rate_func01(t_iter,X);
    end
    
    %[XB, ~] = forward_euler_step(@rate_func01, t, X0, h);
    [t_list, x_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
    [t_list, x_list, h_avg, num_evals] = backward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
    [t_list,x_list,h_avg, num_evals] = implicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);  
    
   
    % Plot Euler vs. Analytical throughout different time step sizes
    figure()
    X = solution01(linspace(tspan(1), tspan(2), (tspan(2)/h_avg))); %Calc Analytical
    plot(linspace(tspan(1), tspan(2), (tspan(2)/h_avg)), X,'b'); % Plot analytical solution over time
    hold on
    for i = 1:4
        plot(t_list, x_list, '.', 'MarkerSize', 20); % Plot Euler step integration over time
        h_ref = h_ref + 0.3;
        [t_list, x_list, h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
        hold on
    end
    ylim([-4, 4]);
    title("Euler step integration over time for different time step sizes");
    xlabel("Time");
    ylabel("X(t)");
    legend("Analytical Solution","0.1 timestep", "0.4 timestep", "0.7 timestep", "1.0 timestep");
    hold off
    
    % Plot Midpoint vs. Analytical throughout different time step sizes
    figure()
    h_ref = 0.1; %reset back to 0.1
    [t_list, x_list, h_avg, num_evals] = backward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref); % re-claim parameters with reset h_ref
    plot(linspace(tspan(1), tspan(2), (tspan(2)/h_avg)), X,'b'); % Plot analytical solution over time
    hold on
    for i = 1:4
        [t_mid_list,x_list_expmid,h_avg_expmid, num_evals_expmid] = explicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
        plot(t_mid_list, x_list_expmid(:, 1), '.', 'MarkerSize', 20); % Plot Explict Midpoint step integration over time
        h_ref = h_ref + 0.2;
        hold on;
    end
    ylim([-4, 4]);
    title("Explicit midpoint step integration over time for different time step sizes");
    xlabel("Time");
    ylabel("X(t)");
    legend("Analytical Solution","0.1 timestep", "0.3 timestep", "0.5 timestep", "0.7 timestep");
    hold off

    % Plot Euler + Midpoint vs. Analytical for h_ref = 0.38 and 0.45
    figure()
    X = solution01(linspace(tspan(1), tspan(2), (tspan(2)/h_avg))); %Calc Analytical
    plot(linspace(tspan(1), tspan(2), (tspan(2)/h_avg)), X,'b'); % Plot analytical solution over time
    hold on
    plot(t_list, x_list,'gO'); % Plot Euler step integration over time
    hold on
    plot(t_mid_list, x_list_expmid,'rO'); % Plot Explict Midpoint step integration over time
    title("Euler and Midpoint step integration over time for href = 0.45");
    xlabel("Time");
    ylabel("X(t)");
    legend("Analytical Solution", "Euler", "Explicit Midpoint");
    hold off
   
    
    % Find local truncation error for Forward Euler, Explicit Midpoint, and Analytical Difference for test 1
    [h_list, analytical_difference, fel_error_list, expmid_error_list, bel_error_list, impmid_error_list] = truncation_error(t_ref, hspan,@rate_func01);
    
    % Print out values for table of local truncation errors with step size
    % (values of p). p values chosen are the first and last values of
    % h_span, or -5 and -1.
    % fprintf("analytical_difference: ");
    % analytical_difference(1)
    % analytical_difference(end)
    % fprintf("fel_error_list: ")
    % fel_error_list(1)
    % fel_error_list(end)
    % fprintf("expmid_error_list: ")
    % expmid_error_list(1)
    % expmid_error_list(end)

    % Find line of best fit coefficients and estimate p-value for test 1
    [ana_p, ~] = loglog_fit(h_list, analytical_difference);
    fprintf("Local analytical p-value: ");
    disp(ana_p);

    [fel_p,fel_k] = loglog_fit(h_list,fel_error_list);
    fel_y_data = fel_k.*((h_list.^fel_p));
    fprintf("Local forward Euler p-value: ");
    disp(fel_p);

    [expmid_p,expmid_k] = loglog_fit(h_list,expmid_error_list);
    expmid_y_data = expmid_k.*((h_list.^expmid_p));
    fprintf("Local explicit Midpoint p-value: ");
    disp(expmid_p);

    [bel_p, ~] = loglog_fit(h_list,bel_error_list);
    fprintf("Local backward Euler p value: ");
    disp(bel_p);

    [impmid_p, ~] = loglog_fit(h_list,impmid_error_list);
    fprintf("Local implicit Midpoint p-value: ");
    disp(impmid_p);

    % % Find local truncation error for Forward Euler, Explicit Midpoint, and Analytical Difference for test 2
    % [h_list2, analytical_difference, fel_error_list2, expmid_error_list2, bel_error_list2, impmid_error_list2] = truncation_error(t_ref, hspan,@rate_func02);
    % 
    % % Estimate p-value for test 2
    % [fel_p2, ~] = loglog_fit(h_list2,fel_error_list2);
    % fprintf("Forward Euler p value: ");
    % disp(fel_p2);
    % 
    % [expmid_p2,~] = loglog_fit(h_list2,expmid_error_list2);
    % fprintf("Explicit Midpoint p-value: ");
    % disp(expmid_p2);
    % 
    % [bel_p2, ~] = loglog_fit(h_list2,bel_error_list2);
    % fprintf("Backward Euler p value: ");
    % disp(bel_p2);
    % 
    % [impmid_p,~] = loglog_fit(h_list2,impmid_error_list2);
    % fprintf("Implicit Midpoint p-value: ");
    % disp(impmid_p);


    % Plot log scale graph of errors vs h_list (all the different step sizes used)
    figure;

    loglog(h_list,fel_error_list,'ro','MarkerFaceColor','r'); hold on;
    %loglog(h_list,fel_y_data,'k','LineWidth',2);

    loglog(h_list,expmid_error_list,'go','MarkerFaceColor','g');
    %loglog(h_list,expmid_y_data,'k','LineWidth',2); hold on;

    loglog(h_list,bel_error_list,'bo');
    %loglog(h_list,bel_y_data,'k','LineWidth',2); hold on;

    loglog(h_list,impmid_error_list,'mo');
    %loglog(h_list,impmid_y_data,'k','LineWidth',2); hold off;

    legend('Forward Euler', 'Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint')
    %loglog(h_list,analytical_difference, 'bo', 'MarkerFaceColor', 'b'); hold off;
    %legend('Forward Euler','Fit line','Explicit Midpoint', '', 'Analytical Solution')
    % Set axes and legend
    xlim([1e-5 1e0])
    title("Local error vs. h list for all four methods")
    xlabel("Time step sizes");
    ylabel("Error");

    % Find Global truncation error for Forward Euler, Explicit Midpoint
    [global_h_list,global_fel_error_list,global_expmid_error_list, global_bel_error_list,global_impmid_error_list,global_analytical_list, tot_evals_fel,tot_evals_expmid, tot_evals_bel, tot_evals_impmid] = global_truncation_error(tspan, hspan, @rate_func01);
 
 
    % Find line of best fit coefficients and estimate p-values for test 1
    [g_fel_p,g_fel_k] = loglog_fit(global_h_list,global_fel_error_list);
    g_fel_y_data = g_fel_k.*((global_h_list.^g_fel_p));
    fprintf("Global forward Euler p-value: ");
    disp(g_fel_p);

    [g_expmid_p,g_expmid_k] = loglog_fit(global_h_list,global_expmid_error_list);
    g_expmid_y_data = g_expmid_k.*((global_h_list.^g_expmid_p));
    fprintf("Global backward Euler p-value: ");
    disp(g_expmid_p);

    [g_bel_p, g_bel_k] = loglog_fit(h_list,bel_error_list);
    g_bel_y_data = g_bel_k.*((global_h_list.^g_bel_p));
    fprintf("GLobal explicit midpoint p-value: ");
    disp(g_bel_p);

    [g_impmid_p, g_impmid_k] = loglog_fit(h_list,impmid_error_list);
    g_impmid_y_data = g_impmid_k.*((global_h_list.^g_impmid_p));
    fprintf("Global implicit Midpoint p-value: ");
    disp(g_impmid_p);


    % Plot log scale graph of errors vs h_list (all the different step sizes used)
    figure;

    loglog(global_h_list,global_fel_error_list,'ro','MarkerFaceColor','r'); hold on;
    loglog(global_h_list,g_fel_y_data,'k','LineWidth',2);

    loglog(global_h_list,global_expmid_error_list,'go','MarkerFaceColor','g');
    loglog(global_h_list,g_expmid_y_data,'k','LineWidth',2); hold on; 

    %loglog(global_h_list,analytical_difference, 'bo', 'MarkerFaceColor', 'b'); hold off;
    legend('Forward Euler','Fit line','Explicit Midpoint', '')
    title("Global error vs. h list")
    xlabel("Time step sizes");
    ylabel("Error");
    

    % Plot log scale graph of errors vs # of calls
    figure;

    loglog(tot_evals_fel,global_fel_error_list,'ro','MarkerFaceColor','r'); hold on;
    loglog(tot_evals_fel,g_fel_y_data,'k','LineWidth',2);

    loglog(tot_evals_expmid,global_expmid_error_list,'go','MarkerFaceColor','g');
    loglog(tot_evals_expmid,g_expmid_y_data,'k','LineWidth',2); hold on; 

    %loglog(tot_evals_bel,global_bel_error_list,'ro','MarkerFaceColor','r'); hold on;
    %loglog(tot_evals_bel,g_bel_y_data,'k','LineWidth',2);

    %loglog(tot_evals_impmid,global_impmid_error_list,'go','MarkerFaceColor','g');
    %loglog(tot_evals_impmid,g_impmid_y_data,'k','LineWidth',2); hold on; 

    %loglog(tot_evals_fel,analytical_difference, 'bo', 'MarkerFaceColor', 'b'); hold off;
    %legend('Forward Euler','Fit line','Explicit Midpoint', '', 'Backward Euler','','Implicit Midpoint', '')
    legend('Forward Euler','Fit line','Explicit Midpoint')
    title('Global Error vs # of Calls')
    xlabel("# of calls");
    ylabel("Global error")

    % Estimate p-values for global truncation error vs. rate function calls for test 1
    [gcall_fel_p,~] = loglog_fit(tot_evals_fel, global_fel_error_list);
    fprintf("Global/rate function calls forward Euler p-value: ");
    disp(gcall_fel_p);

    [gcall_expmid_p,~] = loglog_fit(tot_evals_expmid,global_expmid_error_list);
    fprintf("Global/rate function calls backward Euler p-value: ");
    disp(gcall_expmid_p);

    [gcall_bel_p, ~] = loglog_fit(tot_evals_bel,bel_error_list);
    fprintf("Global/rate function calls explicit midpoint p-value: ");
    disp(gcall_bel_p);

    [gcall_impmid_p, ~] = loglog_fit(tot_evals_impmid,impmid_error_list);
    fprintf("Global/rate function calls implicit Midpoint p-value: ");
    disp(gcall_impmid_p);

    % Plot log scale graph of global truncation error scaling for all four
    % methods
    figure;

    loglog(global_h_list,global_fel_error_list,'ro','MarkerFaceColor','r'); hold on;

    loglog(global_h_list,global_expmid_error_list,'go','MarkerFaceColor','g'); hold on;

    loglog(global_h_list,global_bel_error_list,'co'); hold on;

    loglog(global_h_list,global_impmid_error_list,'mo'); hold on;

    legend('Forward Euler','Explicit Midpoint','Backwards Euler','Implicit Midpoint')
    title("Global error vs. h list")
    xlabel("Time step sizes");
    ylabel("Error");
    hold off; 

end

