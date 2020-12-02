%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
syms theta
a = pi/2-pi/8; g = 9.81;
f1 = -(2*g*(cos(a) + 1))/sin(a);
h = sqrt(-2*(g*(1+cos(theta))+1/2*f1*sin(theta))/(1-1/2*cos(theta)^2));
fun = matlabFunction(1/h);
T1 = abs(integral(fun, pi, a));
x_time_T1 = 1/2*(1/2*f1*T1^2-2-sin(a));
x_dot_time_T1 = f1*T1/2;
test = true;
theta_2_dot = -5;
f2 = (1/2*theta_2_dot^2-g*cos(a))/(sin(a)+1)*2;
h2 = sqrt(-2*(-g*cos(a)-1/2*f2*sin(a)+g*cos(theta)+1/2*f2*sin(theta))/(1-1/2*cos(theta)^2));
fun2 = matlabFunction(1/h2);
T2 = abs(integral(fun2, a, pi/2));
k = 1/2*(1/2*f2*T2^2+1+2*x_time_T1+sin(a)+2*x_dot_time_T1*T2);
f2 = 500;
% while test
%     %f2 = (1/2*theta_2_dot^2-g*cos(a))/(sin(a)+1)*2;
%     f2 = f2 + 4;
%     h2 = sqrt(-2*(-g*cos(a)-1/2*f2*sin(a)+g*cos(theta)+1/2*f2*sin(theta))/(1-1/2*cos(theta)^2));
%     fun2 = matlabFunction(1/h2);
%     T2 = abs(integral(fun2, a, pi/2));
%     k = 1/2*(1/2*f2*T2^2+1+2*x_time_T1+sin(a)+2*x_dot_time_T1*T2)
%     if k > 0
%         test = false;
%     end
% end
omega = 2*pi/5; g = 9.81;
%phi(theta) = -(1+g/omega^2)*log((1+sin(theta))/cos(theta));


Kd = 2;
A = [0 1 0 0;
     0 0 -g 0;
     0 0 0 1;
     Kd 0 -2*g 0];
B = [0;1;0;1];
Q = [100 0 0 0;
     0 1 0 0;
     0 0 10 0;
     0 0 0 1];
r = 0.1;
K = lqr(A,B,Q,r);
P = -pinv((A-B*K)^-1*B);
%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%
%out = sim('pendulum_on_cart', 10);

%%%%%%%%%%%%%%%%%%%%%%% DRAWING STUFF %%%%%%%%%%%%%%
x_box = [-0.3, 0.3, 0.3, -0.3];
y_box = [-0.2, -0.2, 0.2, 0.2];
g = hgtransform;
patch('XData',x_box,'YData',y_box,'FaceColor','yellow','Parent',g)
x_pendulum = [0.1, -0.1, -0.1, 0.1];
y_pendulum = [-1, -1, 0, 0];
p = hgtransform;
patch('XData',x_pendulum,'YData',y_pendulum,'FaceColor','red','Parent',p)
axis equal
xlim([-10 4])
ylim([-3 3])
pt1 = [-1 0 0];
r1 = pi;
r2 = 2*pi;
x_obstacle = [-0.1,0.1,0.1,-0.1];
y_obstacle = [-1, -1, 0, 0];
A =[0 1 0 0;0 0 -9.81 0;0 0 0 1;0 0 -2*9.81 0];
B = [0;3/2;0;2];
Q = eye(4).*[10;1;10;1];
r = 1;
[K, S] = lqr(A,B,Q,r);
P = -pinv((((A-B*K)^-1)*B));
o = hgtransform;
patch('XData',x_obstacle,'YData',y_obstacle,'FaceColor','black','Parent',o)

for t=1:2:length(out.tout)
  g.Matrix = makehgtform('translate',[out.simout(t,1),0,0]);
  p.Matrix = makehgtform('translate',[out.simout(t,1),0,0], ...
                         'zrotate',pi-out.simout(t,2));
  drawnow
end