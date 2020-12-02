clear;
syms psi;
omega_c = 0.5;
psi_0 = 0;
R_d = 1;
psi = @(t) omega_c*t+psi_0; 
g = @(t) 
th = @(t) psi(t) + pi/2;
k1 = 
A = [-k1 0 -k2 -(k3-g*w*sin(th)) -k4;
    0 1 0 0 0;
    k1 0 k2 k3+g*w*cos(th) k4;
    0 0 0 0 1;
    0 0 0 0 0];
B = [0;0;0;0;1];
rank([B A*B A^2*B A^3*B A^4*B])
Q = eye(5)*10;
r = 1;
[X,time] = one_shot_generator(A,B, Q,r,0,2*pi);

ts = append(timeseries(X,time),timeseries(X,time+time(length(time))));