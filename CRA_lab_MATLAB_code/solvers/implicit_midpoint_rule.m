function [sol, time] = implicit_midpoint_rule(f,t0,t_end,h,initial_conditions)
    n = round((t_end-t0)/h)+1;
    time = zeros(1,n);
    sol = zeros(length(initial_conditions),n);
    sol(:,1) = initial_conditions;
    time(1,1) = t0;
    k = zeros(length(initial_conditions),1);
    options = optimoptions('fsolve','Display','none');
    for i = 2:n
        k1 = fsolve(@(x) f(time(1,i-1)+h/2,sol(:,i-1)+h/2*x)-x,k, options);
        time(1,i) = time(1,i-1) + h;
        sol(:,i) = sol(:,i-1) + h*k1;
    end
end

