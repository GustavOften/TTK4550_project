function y = jalla_lqr(x,x_dot,y,y_dot,theta,theta_dot)
    psi = atan2(y, x);
    phi = atan2(y_dot,x_dot);
    omega_c = 0.5;
    R_d = 1;
%     A = [0 1 0 0 0;
%         0 cos(psi)*sin(psi)*omega_c 0 -cos(psi)^2*omega_c -R_d*omega_c*cos(psi);
%         0 0 0 1 0;
%         0 sin(psi)^2*omega_c 0 -cos(psi)*sin(psi)*omega_c -R_d*omega_c*sin(psi);
%         0 0 0 0 0]
%     B = [0;0;0;0;1];
%     rank([B A*B A^2*B A^3*B A^4*B])
%     Q = eye(5)*10;
%     r = 1;
%     K = lqr(A,B,Q,r);
    states = [x-R_d*cos(psi);x_dot+R_d*sin(psi);y-R_d*sin(psi);y_dot-R_d*cos(psi)*omega_c;theta_dot-omega_c];
y = states;
