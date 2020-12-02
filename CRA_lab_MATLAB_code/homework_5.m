%% Øving 5.
close all;
L = -1.3;
g = 9.81;
dx_0 = 0.2;

%% Defining system
sys = @(theta, theta_dot, u) [(u-cos((1/2*sin(theta)*theta_dot^2+cos(theta)*u/2-g*sin(theta))/(1+cos(theta)^2/2))+sin(theta)*theta_dot^2)/2;
                               (1/2*sin(theta)*theta_dot^2+cos(theta)*u/2-g*sin(theta))/(1+cos(theta)^2/2)];

%% Defining motion generator 
m_g = @(x) x + L*sin(x);
der_m_g = @(x) 1 + L*cos(x);
sec_der_m_g = @(x) -L*sin(x);

%% Defining alpha beta gamma
alpha = @(x) cos(m_g(x)) + der_m_g(x);
beta = @(x) sec_der_m_g(x);
gamma = @(x) -g*sin(m_g(x));
xdd = @(x,xd) -1/alpha(x)*(-beta(x)*xd^2-gamma(x));

%% Defining invariance controller
epsilon = @(x,y,x_dot,y_dot) (2+der_m_g(x)*cos(y+m_g(x)))/(cos(y+m_g(x))+der_m_g(x))*(g*sin(y+m_g(x))-sec_der_m_g(x)*x_dot^2)+cos(y+m_g(x))*sec_der_m_g(x)*x_dot^2-sin(y+m_g(x))*(y_dot+der_m_g(x)*x_dot)^2;
psi = @(x,y) (cos(y+m_g(x))^2-2)/(cos(y+m_g(x))+der_m_g(x));


%% Calculating trajectory
f = @(t,x) [x(2);-1/alpha(x(1))*(beta(x(1))*x(2)^2+gamma(x(1)))];
%phase_plot_2_interactive(f,[-pi pi;-1 1],10,'',[100,100],0.1);
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[ts,ys] = ode45(f,[0,30],[0;dx_0], options);
figure
plot(ys(:,1),ys(:,2))
figure
plot(ys(:,1),m_g(ys(:,1)))
i = 50;
length(ys)
while norm(ys(i,:)-ys(1,:)) > 0.005
    i = i+1;
end
phi_dot = ys(1:i+1,:);
diff_diff_x = zeros(1,length(phi_dot));
for i = 1:length(phi_dot)
    diff_diff_x(1,i) = xdd(phi_dot(i,1),phi_dot(i,2));
end
time = ts(1:i);
f = polyfit(time(:,1), phi_dot(:,1),20);
f_x = matlabFunction(poly2sym(f));
f = polyfit(time(:,1), phi_dot(:,2),25);
f_dx = matlabFunction(poly2sym(f));
f = polyfit(time(:,1), diff_diff_x(1,:)',18);
f_ddx = matlabFunction(poly2sym(f));
figure
hold on;
plot(time,phi_dot(:,1));
plot(time,f_x(time));
figure
hold on;
plot(time,phi_dot(:,2));
plot(time,f_dx(time));
figure
hold on;
plot(time,diff_diff_x);
plot(time,f_ddx(time));


%% Defining linearization  
A = @(x,x_dot, x_dot_dot)[-2*x_dot*beta(x)/alpha(x) 2*x_dot*(g*cos(m_g(x))+sin(m_g(x))*x_dot_dot)/alpha(x) 0;
                            0   0   1;
                            0   0   0];
B = @(x,xd) [-2*xd/alpha(x);0;1];

A_plot = zeros(3,3,length(time));
B_plot = zeros(3,1,length(time));
for i = 1:length(time)
    t = time(i);
    A_plot(:,:,i) = A(f_x(t),f_dx(t),f_ddx(t));
    B_plot(:,:,i) = B(f_x(t),f_dx(t));
end


k = 1;
for i = 1:3
    for j = 1:3
        subplot(3,3,k)
        hold on;
        plot(time(:),reshape(A_plot(i,j,:),1,length(time)));
        k = k+1;
    end
end
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('A(\phi)');
k = 1;

figure
for i = 1:3
    subplot(3,1,k)
    hold on;
    plot(time(:),reshape(B_plot(i,1,:),1,length(time)));
    k = k+1;
end
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('B(\phi)')


%% Calculating Riccati.
func_A = @(t) A(f_x(t),f_dx(t),f_ddx(t));
func_B = @(t) B(f_x(t),f_dx(t));
Q =@(t) [1 0 0;0 1 0;0 0 1];
Gamma = 1;
%[X,time] = multi_shot(func_A,func_B,Q,Gamma,0,time(end),256);
%[X,time] = one_shot_generator(func_A,func_B,Q,Gamma,0,time(end));
time(end)
X = sdp_riccati(func_A,func_B,Q,1,0,3,505,5,3);
interpolation_for_X = @(x)[interp1(time,reshape(X(1,1,:),1,length(X)),x) interp1(time,reshape(X(1,2,:),1,length(X)),x) interp1(time,reshape(X(1,3,:),1,length(X)),x);
                           interp1(time,reshape(X(2,1,:),1,length(X)),x) interp1(time,reshape(X(2,2,:),1,length(X)),x) interp1(time,reshape(X(2,3,:),1,length(X)),x);
                           interp1(time,reshape(X(3,1,:),1,length(X)),x) interp1(time,reshape(X(3,2,:),1,length(X)),x) interp1(time,reshape(X(3,3,:),1,length(X)),x)];

g_y_lin= @(x,ddx) sin(m_g(x))+cos(m_g(x))*ddx-sin(m_g(x))-cos(m_g(x));
%diff_I = @(x,y,dy,dx,ddx,v,theta,dtheta,I) -beta(x)/alpha(x)*I+2*dx*g_y_lin(x,ddx)/alpha(x)*y-2*dx/alpha(x)*v;%dx*(2/alpha(x)*(sin(y+m_g(x))+cos(y+m_g(x))*ddx-sin(m_g(x))-cos(m_g(x))*ddx)*y-v-2*beta(x)/alpha(x)*I);                       
v =@(x,dx,time,transverse) -(1/Gamma)*B(x,dx)'*interpolation_for_X(mod(time,time(end)))*transverse;

controller = @(x,dx,x_star,dx_star,time,transverse) epsilon(x_star,transverse(2),dx_star,transverse(3)) + psi(x_star,transverse(2))*v(x,dx,time,transverse);
%invariance_controller = @(x,dx,ddx) 2*ddx+cos(m_g(x))*(sec_der_m_g(x)*dx^2+der_m_g(x))-sin(m_g(x))*(der_m_g(x)*dx)^2 + invariance_controller(x,dx,ddx);

W = @(x,y,ddx,ddx_star,v)-g*(sin(m_g(x))+sin(y+m_g(x)))-(-cos(m_g(x))+cos(y+m_g(x)))*ddx-v;

X = zeros(3,3,length(time));
for i = 1:length(time)
    X(:,:,i) = interpolation_for_X(mod(time(i),time(end)));
end
not_positive = [];
for i = 1:length(X(:,:,:))
    [~,r] = chol(X(:,:,i));
    j = 0;
    if r
        fprintf("X not positive definite");
        not_positive = [not_positive i];
    else
    end
end
not_controllable = [];
for i = 1:length(time)
    %[A,B] = bf.get_linearization(phi(i),bf.function_for_dphi(phi(i)));
    A_ = func_A(time(i));
    B_ = func_B(time(i));
    control = [B_ A_*B_ A_^2*B_];
    if rank(control) ~= 3 
        fprintf("A B not controllable")
        not_controllable = [not_controllable i];
    else
    end
end

figure
k = 1;
X = zeros(3,3,length(time));
for i = 1:length(time)
    X(:,:,i) = interpolation_for_X(mod(time(i),time(end)));
end
for i = 1:3
    for j = 1:3
        subplot(3,3,k)
        hold on;
        plot(time,reshape(X(i,j,:),1,length(time)));
        k = k+1;
    end
end
sgtitle('Riccati sol X');
y = @(theta, x) theta - m_g(x);
y_dot = @(theta_dot , x, dx) theta_dot - der_m_g(x)*dx;


A_plot = zeros(3,3,length(time));
B_plot = zeros(3,1,length(time));
for i = 1:length(time)
    t = time(i);
    A_plot(:,:,i) = func_A(t);
    B_plot(:,:,i) = func_B(t);
end


k = 1;
figure
for i = 1:3
    for j = 1:3
        subplot(3,3,k)
        hold on;
        plot(time(:),reshape(A_plot(i,j,:),1,length(time)));
        k = k+1;
    end
end
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('A(\phi)');
k = 1;

figure
for i = 1:3
    subplot(3,1,k)
    hold on;
    plot(time(:),reshape(B_plot(i,1,:),1,length(time)));
    k = k+1;
end

set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('B(\phi)')