function [eror_list] = truncation_error(analytical_soln, t_ref, hspan)

    h_list = logspace(hspan(1),hspan(2),hspan(3));
    x_approx_list = [];

    XA = analytical_soln(t_ref);
    for i = 1:hspan(3)
        
        [x_approx,~] = forward_euler_step(@analytical_soln,t_ref,XA,hspan(i));
        x_approx_list(end+1) = x_approx;
    end
    error_list = x_approx_list-XA;

end
