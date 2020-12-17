sdpvar a b c d e
constraints = [a<=0];
JRm = bf.J_s/bf.R_b/bf.m_b;
theta_eq = @(x) a*atan(b*sin(2*x)+c*sin(4*x)+d*sin(6*x)+e*sin(8*x)) + x;
dtheta_eq = @(x) a*(2*b*cos(2*x)+4*c*cos(4*x)+6*d*cos(6*x)+e*8*cos(8*x))/((b*sin(2*x)+c*sin(4*x)+d*sin(6*x)+e*sin(8*x))*(b*sin(2*x)+c*sin(4*x)+d*sin(6*x)+e*sin(8*x))+1)+1;
alpha = @(x) (bf.get_ds(x)*[0 0 1]*cross(bf.get_rho(x),bf.get_tau(x))-bf.get_dsf(x)*bf.J_s/bf.R_b/bf.m_b)*dtheta_eq(x)+bf.get_ds(x)^2+bf.get_dsf(x)^2*bf.J_s/bf.R_b^2/bf.m_b;
constraint_func = @(x) (bf.get_ds(x)^2+bf.get_dsf(x)^2*JRm/bf.R_b)/ ...
                (bf.get_ds(x)*([0 0 1]*cross(bf.get_rho(x),bf.get_tau(x))) - bf.get_dsf(x)*JRm);
diff_gamma = @(x) [0 1 0]*(bf.diff_R(theta_eq(x))*dtheta_eq(x)*bf.get_tau(x)*bf.get_ds(x)+...
    bf.R(theta_eq(x))*bf.get_kappa(x)*bf.get_ds(x)^2+bf.R(theta_eq(x))*bf.get_tau(x)*bf.get_dds(x));
objective = 0;
for i = 1:100
    x = (i-1)*pi/100;
    i
    constraints = [constraints, alpha(x) >= 0.001];
    %constraints = [constraints, -35<=diff_gamma(x)/alpha(x)*bf.g <= 35];
end

%diff_gamma = @(x) [cos(theta_eq(x)) -sin(theta_eq(x)) 0]*dtheta_eq(x)*bf.get_tau(x)*bf.get_ds(x)+[sin(theta_eq(x)) cos(theta_eq(x)) 0]*bf.get_kappa(x)*bf.get_ds(x)+[sin(theta_eq(x)) cos(theta_eq(x)) 0]*bf.get_tau(x)*bf.get_dds(x);
%constraint_gamma_zeros(0)

eq1 = pi/2-pi/8; eq2 = pi/8; eq3 = pi/2-pi/8; eq4 = pi/2-pi/16;
%% Making additional eq points be in the desired points:
% constraints = [constraints, -0.001 <= theta_eq(eq1)-atan2(-[0 1 0]*bf.get_tau(eq1),[1 0 0]*bf.get_tau(eq1)) <= 0.001];
% constraints = [constraints, -0.001 <= theta_eq(eq2)-atan2(-[0 1 0]*bf.get_tau(eq2),[1 0 0]*bf.get_tau(eq2)) <= 0.001];
% constraints = [constraints, -0.001 <= theta_eq(eq3)-atan2(-[0 1 0]*bf.get_tau(eq3),[1 0 0]*bf.get_tau(eq3)) <= 0.001];
% constraints = [constraints, -0.001 <= theta_eq(eq4)-atan2(-[0 1 0]*bf.get_tau(eq4),[1 0 0]*bf.get_tau(eq4)) <= 0.001];


% constraints = [constraints, diff_gamma(0) <= 0];
constraints = [constraints, diff_gamma(pi/2) >= 0];

constraints = [constraints, diff_gamma(eq1) <= 0];
% constraints = [constraints, diff_gamma(eq2) <= 0];
% constraints = [constraints, diff_gamma(eq3) >= 0];
% constraints = [constraints, diff_gamma(eq4) <= 0];
%objective = (diff_gamma(0)-diff_gamma(pi/2))+diff_gamma(eq1)-diff_gamma(eq2)+diff_gamma(eq3)-diff_gamma(eq4)
objective = (-diff_gamma(pi/2)-diff_gamma(eq1));%+diff_gamma(eq1)-diff_gamma(eq2)+diff_gamma(eq3)-diff_gamma(eq4)
options = sdpsettings('solver', 'fmincon','fmincon.MaxIter',1000);
sol = optimize(constraints,objective,options);
fprintf("a:= %d: b:= %d: c:= %d: d:= %d:\n",value(a),value(b),value(c),value(d));
sol.info
fprintf("a:= %d: b:= %d: c:= %d: d:= %d: e:= %d:\n",value(a),value(b),value(c),value(d),value(e));
fprintf("a= %d; b= %d; c= %d; d= %d; e= %d;;\n",value(a),value(b),value(c),value(d),value(e));