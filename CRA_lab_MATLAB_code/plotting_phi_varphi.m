close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
n = 3;

phi = linspace(0,2*pi,1000);
X = zeros(3,3,length(phi));
for i = 1:length(phi)
    X(:,:,i) = bf.function_for_X(mod(phi(i),pi));
end
not_positive = [];
for i = 1:length(bf.X(:,:,:))
    d = eig(bf.X(:,:,i))
    if all(d >= 0)
    else
        fprintf("X %d not positive semi definite ", i)
        not_positive = [not_positive i];
    end
end


figure
k = 1;
phi = linspace(0,2*pi,1000);
X = zeros(3,3,length(phi));
for i = 1:length(phi)
    X(:,:,i) = bf.function_for_X(mod(phi(i),pi));
end
for i = 1:3
    for j = 1:3
        subplot(3,3,k)
        hold on;
        plot(phi,reshape(X(i,j,:),1,1000));
        k = k+1;
    end
end
sgtitle('Riccati sol X');
figure 

not_controllable = [];
for i = 1:length(phi)
    %[A,B] = bf.get_linearization(phi(i),bf.function_for_dphi(phi(i)),true);
    A = bf.A(phi(i));
    B = bf.B(phi(i));
    control = [B A*B A^2*B];
    if rank(control) ~= 3
        fprintf("A B not controllable")
        not_controllable = [not_controllable i];
    else
    end
end

plot(phi,bf.function_for_dphi(mod(phi,pi)));
title('$\frac{d\phi}{dt}$','Interpreter','latex')
A = zeros(3,3,length(phi));
B = zeros(3,1,length(phi));
gwy = zeros(3,1,length(phi));
abg = zeros(3,1,length(phi));
delta = zeros(3,1,length(phi));
rho = zeros(3,1,length(phi));
tau = zeros(3,1,length(phi));
bad_rho = zeros(3,1,length(phi));
diff_s = zeros(1,length(phi));
diff_diff_s = zeros(1,length(phi));
kappa = zeros(3,1,length(phi));
theta = bf.get_theta(phi);
u = zeros(length(phi),1);
my_u = zeros(length(phi),1);
p = zeros(length(phi),1);
dp = zeros(length(phi),1);
normal_vector = zeros(3,1,length(phi));
for i = 1:length(phi)
    A(:,:,i) = bf.A(phi(i));
    B(:,:,i) = bf.B(phi(i));
    gwy(:,1,i) = bf.get_g_w_y_on_trajectory(phi(i),bf.function_for_dphi(mod(phi(i),pi)));
    abg(:,1,i) = bf.alpha_beta_gamma(phi(i));
    tau(:,1,i) = bf.get_tau(phi(i));
    rho(:,1,i) = bf.get_rho(phi(i));
    delta(:,1,i) = bf.get_delta(phi(i));
%     bad_rho(:,1,i) = bf.bad_rho(phi(i));
    diff_s(1,i) = bf.get_ds(phi(i));
    diff_diff_s(1,i) = bf.get_dds(phi(i));
    kappa(:,1,i) = bf.get_kappa(phi(i));
%    u(i) = bf.get_u([bf.theta(phi(i));phi(i)],[bf.diff_theta(phi(i))*bf.function_for_dphi(phi(i));bf.function_for_dphi(phi(i))],[0;0;0;0]);
 %   my_u(i) = bf.get_my_u([bf.theta(phi(i));phi(i)],[bf.diff_theta(phi(i))*bf.function_for_dphi(phi(i));bf.function_for_dphi(phi(i))],[0;0;0;0]);
    p(i) = bf.get_p(phi(i));
    dp(i) = bf.get_dp(phi(i));
    normal_vector(:,1,i) = bf.get_normal_vector(phi(i));
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
        plot(reshape(bf.X(i,j,:),1,length(bf.X)));
        k = k+1;
        set(gca,'XTick',0:pi/2:2*pi)    
        set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
        grid on;
    end
end
sgtitle('X without any interpolation');

figure
plot(phi,p);
title("p($\varphi$)");

figure
plot(phi,dp);
title("dp($\varphi$)");

figure
for i = 1:3
    subplot(3,1,i)
    hold on;
    plot(phi,reshape(normal_vector(i,1,:),1,length(phi)));
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    grid on;
    k = k+1;
end

sgtitle('Normal vector to frame', 'Interpreter','latex');

