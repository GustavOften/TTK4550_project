syms t 
o = 1;
R_d = 1;
psi_0 = 0;
psi(t) = o*t+psi_0;
x(t) = R_d*cos(psi(t));
y(t) = R_d*sin(psi(t));
theta(t) = psi(t)+pi/2;
x_dot(t) = diff(x,t);
y_dot(t) = diff(y,t);
theta_dot(t) = diff(theta,t);
A = [0 1 0 0 0;
    0 cos(psi(t))*sin(psi(t))*o 0 -cos(psi(t))^2*o -R_d*o*cos(psi(t));
    0 0 0 1 0;
    0 sin(psi(t))^2*o 0 -cos(psi(t))*sin(psi(t))*o -R_d*o*sin(psi(t));
    0 0 0 0 0];
B = [0;0;0;0;1];
Q = eye(5);
r = 1;
[K, S] = lqr(A,B,Q,r);

