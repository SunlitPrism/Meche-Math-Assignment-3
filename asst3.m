%This function computes the value of X at the next time step
%using the Forward Euler approximation
function asst3()
    %define constants
    t = 1;
    h = 1; %timestep
    X = [0,1]; %inital X
    
    dXdt = rate_func01(t,X);
    XA = solution01(t);

    disp(dXdt)
    

    [XB,~] = forward_euler_step(rate_func01,t,X,h);
end

