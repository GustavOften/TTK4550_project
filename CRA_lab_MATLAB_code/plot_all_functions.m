close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
n = 4;

phi = linspace(0,2*pi,1000);
X = zeros(4,4,length(phi));
for i = 1:length(phi)
    X(:,:,i) = model.function_for_X(mod(phi(i),pi));
end
not_positive = [];
for i = 1:length(model.X(:,:,:))
    d = eig(model.X(:,:,i))
    if all(d >= 0)
    else
        fprintf("X %d not positive semi definite ", i)
        not_positive = [not_positive i];
    end
end


figure
k = 1;
phi = linspace(0,2*pi,1000);
X = zeros(4,4,length(phi));
for i = 1:length(phi)
    X(:,:,i) = model.function_for_X(mod(phi(i),pi));
end
for i = 1:4
    for j = 1:4
        subplot(4,4,k)
        hold on;
        plot(phi,reshape(X(i,j,:),1,1000));
        k = k+1;
    end
end
sgtitle('Riccati sol X');
figure 

not_controllable = [];
for i = 1:length(phi)
    %[A,B] = model.get_linearization(phi(i),model.function_for_dphi(phi(i)),true);
    A = model.A(phi(i));
    B = model.B(phi(i));
    control = [B A*B A^2*B];
    if rank(control) ~= 4 
        fprintf("A B not controllable")
        not_controllable = [not_controllable i];
    else
    end
end

plot(phi,model.function_for_dphi(mod(phi,pi)));
title('$\frac{d\phi}{dt}$','Interpreter','latex')
A = zeros(4,4,length(phi));
B = zeros(4,1,length(phi));
gwy = zeros(3,1,length(phi));
abg = zeros(3,1,length(phi));
delta = zeros(3,1,length(phi));
rho = zeros(3,1,length(phi));
tau = zeros(3,1,length(phi));
bad_rho = zeros(3,1,length(phi));
diff_s = zeros(1,length(phi));
diff_diff_s = zeros(1,length(phi));
kappa = zeros(3,1,length(phi));
theta = model.theta(phi);
u = zeros(length(phi),1);
my_u = zeros(length(phi),1);

for i = 1:length(phi)
    A(:,:,i) = model.A(phi(i));
    B(:,:,i) = model.B(phi(i));
    gwy(:,1,i) = model.get_g_w_y_on_trajectory(phi(i),model.function_for_dphi(mod(phi(i),pi)));
    abg(:,1,i) = model.alpha_beta_gamma(phi(i));
    tau(:,1,i) = model.tau(phi(i));
    rho(:,1,i) = model.rho(phi(i));
    delta(:,1,i) = model.delta(phi(i));
    bad_rho(:,1,i) = model.bad_rho(phi(i));
    diff_s(1,i) = model.ds(phi(i));
    diff_diff_s(1,i) = model.dds(phi(i));
    kappa(:,1,i) = model.kappa(phi(i));
    u(i) = model.get_u([model.theta(phi(i));phi(i)],[model.diff_theta(phi(i))*model.function_for_dphi(phi(i));model.function_for_dphi(phi(i))],[0;0;0;0]);
    my_u(i) = model.get_my_u([model.theta(phi(i));phi(i)],[model.diff_theta(phi(i))*model.function_for_dphi(phi(i));model.function_for_dphi(phi(i))],[0;0;0;0]);
end
figure
hold on;
plot(phi,u);
plot(phi,my_u);
legend("u","my u")
title("U");

figure

k = 1;
for i = 1:3
    for j = 1:3
        subplot(3,3,k)
        hold on;
        plot(phi,reshape(A(i,j,:),1,length(phi)));
        k = k+1;
    end
end
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('A(\phi)');
k = 1;

figure
for i = 1:3
    subplot(3,1,k)
    hold on;
    plot(phi,reshape(B(i,1,:),1,length(phi)));
    k = k+1;
end
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('B(\phi)')

figure
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(phi,reshape(gwy(i,1,:),1,length(phi)));
    k = k+1;
end
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
sgtitle('$g_w, g_y, g_{\dot{y}}$ of $\phi$', 'Interpreter','latex');

figure
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(phi,reshape(abg(i,1,:),1,length(phi)));
    k = k+1;
end
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
grid on
sgtitle('$\alpha, \beta, \gamma$ of $\phi$', 'Interpreter','latex');

figure
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(phi,reshape(delta(i,1,:),1,length(phi)));
    plot(phi,reshape(rho(i,1,:),1,length(phi)));
    plot(phi,reshape(bad_rho(i,1,:),1,length(phi)));
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    legend('$\delta$', '$\rho$','$\rho_{bad}$','Interpreter','latex');
    grid on;
    k = k+1;
end
sgtitle('\delta and \rho')


figure
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(phi,reshape(tau(i,1,:),1,length(phi)));
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    grid on;
    k = k+1;
end
sgtitle('\tau')


figure
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(phi,reshape(kappa(i,1,:),1,length(phi)));
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    grid on;
    k = k+1;
end

sgtitle('\kappa')

figure
subplot(2,1,1)
hold on;
plot(phi,diff_s(:))
title('diff_s')
subplot(2,1,2)
plot(phi,diff_diff_s(:))
title('diff_diff_s')
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
figure
plot(phi,theta(1,:));
title('\Theta');

k = 1
figure
for i = 1:3
    for j = 1:3
        subplot(3,3,k)
        hold on;
        plot(reshape(model.X(i,j,:),1,length(model.X)));
        k = k+1;
        set(gca,'XTick',0:pi/2:2*pi)    
        set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
        grid on;
    end
end

sgtitle('X without any interpolation');

% Check that diff Theta is good.
val = zeros(1,length(phi));
for i = 1:length(phi)
    rhoxtau = cross(model.rho(phi(i)),model.tau(phi(i)));
    val = -(1+model.J_s/model.m_b/model.R_b^2)/(rhoxtau(3)-model.J_s/model.m_b/model.R_b)*model.diff_s(phi(i));
end
figure
plot(phi,val)
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})


% Plot alpha from phi = 0..2*pi
figure
plot(phi,atan2(-tau(2,:),tau(1,:)));
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
figure
plot(phi,sin(atan2(-tau(2,:),tau(1,:))))
figure
plot(phi,cos(atan2(-tau(2,:),tau(1,:))))


figure
plot(phi,reshape(X(3,3,:),1000,1).*reshape(B(3,1,:),1000,1)+reshape(X(2,3,:),1000,1))
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})

n = zeros(2,length(phi));
for i = 1:length(phi)
    n(:,i) = [sin(model.alpha(phi(i)));cos(model.alpha(phi(i)))];
end
figure
for i = 1:2
    subplot(2,1,i)
    hold on;
    plot(phi,n(i,:));
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
end

