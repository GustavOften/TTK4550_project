k = linspace(0,2*pi,500);

rho = zeros(3,500);
delta = zeros(3,500);
rho_real = zeros(3,500);
alpha = zeros(1,500);
tau = zeros(3,500);
for i = 1:500
   rho(:,i) = bf.rho(k(i));
   delta(:,i) = bf.delta(k(i));
   rho_real(:,i) = bf.delta(k(i))+bf.R(pi/2)*bf.tau_delta(k(i))*bf.R_b;
   tau(:,i) = bf.tau_delta(k(i));
   alpha(1,i) = acos(rho_real(:,i)'*delta(:,i)/(norm(rho_real(:,i))*norm(delta(:,i))));
end
figure
hold on
plot(k, rho(1,:));
plot(k, rho_real(1,:))

figure
hold on
plot(k, rho(2,:));
plot(k, rho_real(2,:))

figure
plot(k,alpha)

figure
plot(k,tau)

alpha = @(phi) atan2(bf.delta)
