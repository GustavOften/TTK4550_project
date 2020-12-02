function [sol, time] = matrix_initial_value_problem(A,shape,t_0,T,time_step,initial_conditions)
    f = @(t,x) reshape(A(t)*reshape(x,shape(1),shape(2)),[shape(1)*shape(2),1]);
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    %sol = gls(f,JAC,tspan,reshape(initial_conditions,shape(1)*shape(2),1),2);
    [time,sol] = ode113(f,[t_0 T], reshape(initial_conditions,shape(1)*shape(2),1),options);
    sol = sol';
    %[sol, time] = gauss_order_4(f,t_0,T,time_step,reshape(initial_conditions,shape(1)*shape(2),1));
end