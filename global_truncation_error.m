% function [eror_list] = truncation_error(analytical_soln, t_ref, hspan)
% 
%     h_list = logspace(hspan(1),hspan(2),hspan(3));
%     x_approx_list = [];
% 
%     XA = analytical_soln(t_ref);
%     for i = 1:hspan(3)
% 
%         [x_approx,~] = forward_euler_step(@analytical_soln,t_ref,XA,hspan(i));
%         x_approx_list(end+1) = x_approx;
%     end
%     error_list = x_approx_list-XA;
% 
% end

function [h_list,g_fel_error_list,g_expmid_error_list, g_bel_error_list,g_impmid_error_list , g_analytical_difference, tot_evals_fel,tot_evals_expmid, tot_evals_bel] = global_truncation_error(tspan, hspan, test_function)

    h_list = logspace(hspan(1), hspan(2), hspan(3));
    t0 = tspan(1);
    tf = tspan(2);

    x_approx_fel_list = []; % Forward Euler Local
    x_approx_expmid_list = []; % Explicit Midpoint
    x_analytical_list = []; % Analytical Solution
    x_approx_bel_list = []; % Backward Euler 
    x_approx_impmid_list = []; % Implicit Midpoint
    global tot_evals_fel;
    global tot_evals_expmid;
    global tot_evals_bel;
    global tot_evals_impmid;
    %global tot_evals_analytical;
    tot_evals_fel = [];
    tot_evals_expmid = [];
    %tot_evals_analytical = [];
    tot_evals_bel = [];
    tot_evals_impmid = [];

    X0 = solution01(t0);
    for i = 1:(length(h_list)) % subtracting by two to make error_list the same size as t_list in asst3. Fix later.
        
        [~,x_approx_fel, ~, fel_evals] = forward_euler_fixed_step_integration(test_function, tspan, X0, h_list(i));
        [~,x_approx_expmid,~, expmid_evals] = explicit_midpoint_fixed_step_integration(test_function, tspan, X0, h_list(i));
        [~,x_approx_bel, ~, bel_evals] = backward_euler_fixed_step_integration(test_function, tspan, X0, h_list(i));
        [~,x_approx_impmid, ~, impmid_evals] = implicit_midpoint_fixed_step_integration(test_function, tspan, X0, h_list(i));
        x_approx_fel_list(end+1) = x_approx_fel(end);
        x_approx_expmid_list(end+1) = x_approx_expmid(end);
        x_analytical_list(end+1) = solution01(tf);
        x_approx_bel_list(end+1) = x_approx_bel(end);
        x_approx_impmid_list(end+1) = x_approx_impmid(end);
        tot_evals_fel(end+1) = fel_evals;
        tot_evals_expmid(end+1) = expmid_evals;
        tot_evals_bel(end+1) = bel_evals;
        tot_evals_impmid(end+1) = impmid_evals;
    end
    g_fel_error_list = abs(x_approx_fel_list - x_analytical_list);
    g_expmid_error_list = abs(x_approx_expmid_list - x_analytical_list);
    g_analytical_difference = abs(x_analytical_list - x_analytical_list);
    g_bel_error_list = abs(x_approx_bel_list - x_analytical_list);
    g_impmid_error_list = abs(x_approx_impmid_list - x_analytical_list);
end