function [sol, time] = gauss_order_4(f,t0,t_end,h,initial_conditions)
    number_of_iterations = round((t_end-t0)/h)+1;
    s = length(initial_conditions);
    time = zeros(1,number_of_iterations);
    sol = zeros(s,number_of_iterations);
    sol(:,1) = initial_conditions;
    time(1,1) = t0;
    k = zeros(2*s,1);
    options = optimoptions('fsolve','Display','none');
    c = [1/2-sqrt(3)/6;1/2+sqrt(3)/6];
    b = [1/2 1/2];
    a = [1/4 1/4-sqrt(3)/6;1/4+sqrt(3)/6 1/4];
        options = optimoptions('fsolve','Display','none');
    for i = 2:number_of_iterations
        k = fsolve(@(x) [f(time(1,i-1)+h*c(1), ...
                         sol(:,i-1)+h*(a(1,1)*x(1:s)+a(1,2)*x(s+1:2*s)))-x(1:s);
                         f(time(1,i-1)+h*c(2), ...
                         sol(:,i-1)+h*(a(2,1)*x(1:s)+a(2,2)*x(s+1:2*s)))-x(s+1:2*s)], ...
                         k, options);
        sol(:,i) = sol(:,i-1) + h*(b(1)*k(1:s) + b(2)*k(s+1:2*s));
        time(1,i) = time(1,i-1) + h;
    end
end