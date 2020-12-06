sdpvar a b c d
constraints = [-0.01<=a>=0.01, -0.01<=b>=0.01, -0.01<=c>=0.01 -0.01<=d>=0.01];
objective;
JRm = bf.J_s/bf.R_b/bf.m_b;
theta_eq = @(x) a*atan(b*sin(2*x)+c*sin(4*x)+d*sin(6*x)) + x;
dtheta_eq = @(x)((a*(2*b*cos(2*x)+4*c*cos(4*x)+6*d*cos(6*x)))/((b*sin(2*x)+c*sin(4*x)+d*sin(6*x))*(b*sin(2*x)+c*sin(4*x)+d*sin(6*x))+1)+1);
alpha = @(x) (bf.get_ds(x)*[0 0 1]*cross(bf.get_rho(x),bf.get_tau(x))-bf.get_dsf(x)*bf.J_s/bf.R_b/bf.m_b)*dtheta_eq(x)+bf.get_ds(x)^2+bf.get_dsf(x)^2*bf.J_s/bf.R_b^2/bf.m_b;
constraint_func = @(x) (bf.get_ds(x)^2+bf.get_dsf(x)^2*JRm/bf.R_b)/ ...
                (bf.get_ds(x)*([0 0 1]*cross(bf.get_rho(x),bf.get_tau(x))) - bf.get_dsf(x)*JRm);
diff_gamma = @(x) [0 1 0]*(bf.diff_R(theta_eq(x))*dtheta_eq(x)*bf.get_tau(x)*bf.get_ds(x)+...
    bf.R(theta_eq(x))*bf.get_kappa(x)*bf.get_ds(x)^2+bf.R(theta_eq(x))*bf.get_tau(x)*bf.get_dds(x));

for i = 1:50
    x = (i-1)*pi/100;
    i
    constraints = [constraints, alpha(x) >= 0.01];
    %constraints = [constraints, -30 <=bf.g*diff_gamma(x)/alpha(x) <= 30];
end

%diff_gamma = @(x) [cos(theta_eq(x)) -sin(theta_eq(x)) 0]*dtheta_eq(x)*bf.get_tau(x)*bf.get_ds(x)+[sin(theta_eq(x)) cos(theta_eq(x)) 0]*bf.get_kappa(x)*bf.get_ds(x)+[sin(theta_eq(x)) cos(theta_eq(x)) 0]*bf.get_tau(x)*bf.get_dds(x);
%constraint_gamma_zeros(0)


eq1 = pi/8; eq2 = pi-pi/8; eq3 = pi/2 - pi/4 ; eq4 = pi/2 + pi/4;


%% Making eq. point 2 and three be in eq1 and eq2:
constraints = [constraints, -0.01 <= theta_eq(eq1)-atan2(-[0 1 0]*bf.get_tau(eq1),[1 0 0]*bf.get_tau(eq1)) <= 0.01];
constraints = [constraints, -0.01 <= theta_eq(eq2)-mod(atan2(-[0 1 0]*bf.get_tau(eq2),[1 0 0]*bf.get_tau(eq2)),2*pi)  <= 0.01];
constraints = [constraints, -0.01 <= theta_eq(eq3)-atan2(-[0 1 0]*bf.get_tau(eq3),[1 0 0]*bf.get_tau(eq3)) <= 0.01];
constraints = [constraints, -0.01 <= theta_eq(eq4)-mod(atan2(-[0 1 0]*bf.get_tau(eq4),[1 0 0]*bf.get_tau(eq4)),2*pi)  <= 0.01];


constraints = [constraints, diff_gamma(0) <= 0 ];
constraints = [constraints, diff_gamma(pi/2) >= 0];

% %% Making eq. points 2 and 3 sadles.
constraints = [constraints, diff_gamma(eq1) >= 0];
constraints = [constraints, diff_gamma(eq3) <= 0];
constraints = [constraints, diff_gamma(eq2) >= 0];
constraints = [constraints, diff_gamma(eq4) <= 0];

%ops = sdpsettings('solver','sdpt3');
objective =(-diff_gamma(0)-diff_gamma(pi/2));%+diff_gamma(eq1)+diff_gamma(eq2)+diff_gamma(eq3)+diff_gamma(eq4))^2;
sol = optimize(constraints);
fprintf("a:= %d: b:= %d: c:= %d: d:= %d:\n",value(a),value(b),value(c),value(d));
fprintf("a= %d; b= %d; c= %d; d= %d;",value(a),value(b),value(c),value(d));