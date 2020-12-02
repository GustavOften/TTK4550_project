clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Creating LTI with periodic solution.
syms t;
omega = 1;
A = [0 1 0 0;0 0 1 0;0 0 0 1; 0 0 0 0];
A = [0 1;0 0];
B = [0;0;0;1];
B = [0;1];
Q = eye(4);
Q = eye(2);
R = 1;
[K,S] = lqr(A,B,Q,R);
transform = [cos(t*omega) -sin(t*omega);sin(t*omega) cos(t*omega)];
P = [transform zeros(2);zeros(2) transform];
P = transform;
new_A = diff(P,t)*P^-1 + P*A*P^-1;
new_B = P*B;
new_Q = transpose(P)^-1*Q*P^-1;
%2. Solving.
[X,time] = one_shot_generator(new_A,new_B, new_Q,R,0,pi);
%plot(t,y);
sim_A = matlabFunction(new_A);
sim_B = matlabFunction(new_B);

ts = append(timeseries(X,time));
Phi = multi_shot(new_A,new_B,Q,R,0,pi,10);
x11 = X(1,1,:)';
x12 = X(1,2,:)';
x21 = X(2,1,:)';
x22 = X(2,2,:)';




