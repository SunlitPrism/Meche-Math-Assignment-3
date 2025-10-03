%This function computes the value of X at the next time step
%using the Forward Euler approximation
function asst3()
    %define constants
    t = linspace(0,5,6);
    h = 0.5; %timestep

    longX = [];
    longdX = [];

    for t_iter = 1:length(t)-1
        X = solution01(t_iter)
        dXdt = rate_func01(t_iter,X)
    end
    
    

    [XB,~] = forward_euler_step(@rate_func01,t,X,h)
end

