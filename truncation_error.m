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

function [h_list,analytical_difference,fel_error_list,expmid_error_list, bel_error_list, impmid_error_list] = truncation_error(t_ref, hspan, test_function)

    h_list = logspace(hspan(1), hspan(2), hspan(3));
    x_approx_fel_list = []; % Forward Euler Local
    x_approx_expmid_list = []; % Explicit Midpoint
    x_analytical_list = []; % Analytical Solution
    x_approx_bel_list = []; % Backward Euler Local
    x_approx_impmid_list = []; % Implicit Midpoint

    XA = solution01(t_ref);
    for i = 1:(length(h_list)) % subtracting by two to make error_list the same size as t_list in asst3. Fix later.
        
        [x_approx_fel,~] = forward_euler_step(test_function,t_ref,XA,h_list(i));
        [x_approx_expmid,~] = explicit_midpoint_step(test_function,t_ref,XA,h_list(i));
        [x_approx_bel, ~] = backward_euler_step(test_function,t_ref,XA,h_list(i));
        [x_approx_impmid, ~] = implicit_midpoint_step(test_function,t_ref,XA,h_list(i));

        x_approx_fel_list(end+1) = x_approx_fel;
        x_approx_expmid_list(end+1) = x_approx_expmid;
        x_analytical_list(end+1) = solution01(t_ref+h_list(i));
        x_approx_bel_list(end+1) = x_approx_bel;
        x_approx_impmid_list(end+1) = x_approx_impmid;
    end
    analytical_difference = abs(x_analytical_list-XA); % |X(t + h) − X(t)|
    fel_error_list = abs(x_approx_fel_list - x_analytical_list);
    expmid_error_list = abs(x_approx_expmid_list - x_analytical_list);
    bel_error_list = abs(x_approx_bel_list - x_analytical_list);
    impmid_error_list = abs(x_approx_impmid_list - x_analytical_list);

end