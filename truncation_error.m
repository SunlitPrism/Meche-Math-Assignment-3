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

function [error_list] = truncation_error(t_ref, hspan)

    h_list = logspace(hspan(1), hspan(2), hspan(3));
    x_approx_list = [];
    x_analytical_list = [];

    XA = solution01(t_ref);
    for i = 1:(length(h_list) - 2) % subtracting by two to make error_list the same size as t_list in asst3. Fix later.
        
        [x_approx,~] = forward_euler_step(@rate_func01,t_ref,XA,h_list(i));
        x_approx_list(end+1) = x_approx;
        x_analytical_list(end+1) = solution01(t_ref+h_list(i));
    end
    error_list = abs(x_approx_list - x_analytical_list);
    error_list
    analytical_difference = abs(x_analytical_list-XA); % |X(t + h) − X(t)|
    % next step: plot as a function of hspan

end