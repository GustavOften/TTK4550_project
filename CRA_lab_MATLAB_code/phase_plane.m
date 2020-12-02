y1 = linspace(-pi,2*pi,100);
y2 = linspace(0,2,100);
[x,y] = meshgrid(y1,y2);
syms t phi(t) phi_dot(t)
%f = @(t,Y) [Y(2);1];%-1/((3*(abs(sin(Y(1))*(4860*sin(Y(1))^2 - 1529))^2 + abs(cos(Y(1))*(4860*cos(Y(1))^2 - 6571))^2)^(1/2)*((7507068131429461*(abs(sin(Y(1))*(4860*sin(Y(1))^2 - 1529))^2 + abs(cos(Y(1))*(4860*cos(Y(1))^2 - 6571))^2)^(1/2))/90071992547409920000 + ((43433820*cos(Y(1))^4 - 169740441*cos(Y(1))^2 + 219200*(abs(31*sin(Y(1)) - 81*sin(Y(1))^3)^2 + abs(104*cos(Y(1)) - 81*cos(Y(1))^3)^2)^(1/2) + 165384150)*(113419434*cos(2*Y(1)) + 24355332*cos(2*Y(1))^2 + 1205766*cos(2*Y(1))^3 - 19683*cos(2*Y(1))^4 - 111710449))/(19860000*(abs(104*cos(Y(1)) - 81*cos(Y(1))^3)^2 + abs(31*sin(Y(1)) - 81*sin(Y(1))^3)^2)^(1/2)*(35478*cos(2*Y(1))^3 - 21963888*cos(2*Y(1))^2 - 39455478*cos(2*Y(1)) + 19683*cos(2*Y(1))^4 + 82524205))))/20000000)*(Y(2)^2*((3*(abs(sin(Y(1))*(4860*sin(Y(1))^2 - 1529))^2 + abs(cos(Y(1))*(4860*cos(Y(1))^2 - 6571))^2)^(1/2)*((7507068131429461*((204201*sin(2*Y(1)))/10000000 + (19683*sin(4*Y(1)))/2000000))/(9007199254740992*((sin(Y(1))^2*(4860*sin(Y(1))^2 - 1529)^2)/(400000000*sign(1529*sin(Y(1)) - 4860*sin(Y(1))^3)^2) + (cos(Y(1))^2*(4860*cos(Y(1))^2 - 6571)^2)/(400000000*sign(6571*cos(Y(1)) - 4860*cos(Y(1))^3)^2))^(1/2)) + ((43433820*cos(Y(1))^4 - 169740441*cos(Y(1))^2 + 219200*(abs(31*sin(Y(1)) - 81*sin(Y(1))^3)^2 + abs(104*cos(Y(1)) - 81*cos(Y(1))^3)^2)^(1/2) + 165384150)*((625240800692505*sin(4*Y(1)))/16 - (3605626278600401*sin(2*Y(1)))/8 - (2422380005784549*sin(6*Y(1)))/64 + (4208007806379*sin(8*Y(1)))/4 + (5417885822602753*sin(10*Y(1)))/32768 + (3922566021*sin(12*Y(1)))/16 + (2035950471*sin(14*Y(1)))/64))/(827500*(abs(104*cos(Y(1)) - 81*cos(Y(1))^3)^2 + abs(31*sin(Y(1)) - 81*sin(Y(1))^3)^2)^(1/2)*(19683*cos(2*Y(1))^4 + 35478*cos(2*Y(1))^3 - 21963888*cos(2*Y(1))^2 - 39455478*cos(2*Y(1)) + 82524205)^2) + ((1620*cos(Y(1))^2 - 3331)*((73*sin(Y(1) - conj(Y(1))))/2 - (81*sin(Y(1) - 3*conj(Y(1))))/4 + (27*sin(2*real(Y(1))))/4)*(113419434*cos(2*Y(1)) + 24355332*cos(2*Y(1))^2 + 1205766*cos(2*Y(1))^3 - 19683*cos(2*Y(1))^4 - 111710449)^2)/(20000*(abs(cos(Y(1))*(81*cos(Y(1))^2 - 104))^2 + abs(sin(Y(1))*(81*sin(Y(1))^2 - 31))^2)^(1/2)*(19683*cos(2*Y(1))^4 + 35478*cos(2*Y(1))^3 - 21963888*cos(2*Y(1))^2 - 39455478*cos(2*Y(1)) + 82524205)^2)))/20000000)+((79461*sin(2*Y(1))*(abs(sin(Y(1))*(4860*sin(Y(1))^2 - 1529))^2 + abs(cos(Y(1))*(4860*cos(Y(1))^2 - 6571))^2)^(1/2))/(2000000000*(abs(cos(Y(1))*(81*cos(Y(1))^2 - 104))^2 + abs(sin(Y(1))*(81*sin(Y(1))^2 - 31))^2)^(1/2))))];    
%f = @(t, Y) [Y(2);0];
a = bf.alpha_beta_gamma(phi(t));
alpha = a(1); beta = a(2); gamma = a(3);
vars = [phi(t),phi_dot(t)];
%%f = [phi_dot(t);-1/alpha*(beta*phi_dot(t)^2+gamma)];
f = [phi_dot(t);-1/alpha*(beta*phi_dot(t)^2+gamma)];
syms Phi Phi_dot
a = bf.alpha_beta_gamma(Phi);
alpha = a(1); beta = a(2); gamma = a(3);
f_plane = [Phi_dot;-1/alpha*(beta*Phi_dot^2+gamma)];
der = matlabFunction(f_plane);

%f = [phi_dot(t);0];
odefun = odeFunction(f,vars);
u = zeros(size(x));
v = zeros(size(x));
t = 0;
for i = 1:numel(x)
    derivatives = der(x(i),y(i));
    u(i) = derivatives(1);
    v(i) = derivatives(2);
end
figure
h = quiver(x,y,u,v,'r'); figure(gcf)
set(h,'AutoScale','on', 'AutoScaleFactor', 1)
xlabel('y_1')
ylabel('y_2')
%axis tight equal;
hold on;
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
for y20 = 1.5
    [ts,ys] = ode45(odefun,[0,10],[0;y20], options);
    plot(ys(:,1),ys(:,2))
    plot(ys(1,1),ys(1,2),'bo') % starting point
    plot(ys(end,1),ys(end,2),'ks') % ending point
end
i = 1;
phi_dot = [0 0];
length(ys)
while ys(i,1) <= 3*pi
    phi_dot(i,:) = ys(i,:);
    i = i+1;
   
    if i >= length(ys)
        break
    end
end
phi_dot(i,:) = ys(i,:);
fprintf("Number of elements in phi_dot: ");
length(phi_dot)
%%%%%%%%%%%% Polynomial fitted %%%%%%%%%%%%%%%%
p = polyfit(phi_dot(:,1),phi_dot(:,2), 20);
function_for_dphi = matlabFunction(poly2sym(p));
f = fit(phi_dot(:,1), phi_dot(:,2), 'sin9');
sin_fitted_function = @(phi) f.a1*sin(f.b1*phi+f.c1)+f.a2*sin(f.b2*phi+f.c2)+f.a3*sin(f.b3*phi+f.c3)+f.a4*sin(f.b4*phi+f.c4)+f.a5*sin(f.b5*phi+f.c5)+f.a6*sin(f.b6*phi+f.c6)+f.a7*sin(f.b7*phi+f.c7)+f.a8*sin(f.b8*phi+f.c8)+f.a9*sin(f.b9*phi+f.c9);
%%%%%%%%%%%% Interpolated %%%%%%%%%%%%%%%%%%%%%
interpolation_for_dphi = @(phi) interp1(phi_dot(:,1),phi_dot(:,2), phi);
%plot(ys(:,1),function_for_dphi(mod(ys(:,1),pi)));
plot(ys(:,1),sin_fitted_function(mod(ys(:,1),pi)));
syms phi dphi
% Using the curvfitted function since ode45 has problems with
% interpolations.
[A, B] = bf.get_linearization(phi, sin_fitted_function(phi));
%dphi = @(phi) phi_dot(mod(round(phi/pi*length(phi_dot)),length(phi_dot))+1,2),phi;
Q = eye(3); R = 10;
[X,phi] = one_shot_generator(A,B,Q,R,0,pi);
% poly_number = 20;
% px = zeros(3,3,poly_number+1);
% 
% for i = 1:3
%     for j = 1:3
%         px(i,j,:) = polyfit(phi,reshape(X(i,j,:),1,length(X)),poly_number);
%     end
% end
% 
% %%%%%%%%%%%% Curvfitted polynomial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function_for_X  = matlabFunction([poly2sym(px(1,1,:)) poly2sym(px(1,2,:)) poly2sym(px(1,3,:));
%                                   poly2sym(px(2,1,:)) poly2sym(px(2,2,:)) poly2sym(px(2,3,:));
%                                   poly2sym(px(3,1,:)) poly2sym(px(3,2,:)) poly2sym(px(3,3,:))]);
% 
%                               
% %%%%%%%%%%%%%%%% Interpolated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation_for_X = @(x)[interp1(phi,reshape(X(1,1,:),1,length(X)),x) interp1(phi,reshape(X(1,2,:),1,length(X)),x) interp1(phi,reshape(X(1,3,:),1,length(X)),x);
%                            interp1(phi,reshape(X(2,1,:),1,length(X)),x) interp1(phi,reshape(X(2,2,:),1,length(X)),x) interp1(phi,reshape(X(2,3,:),1,length(X)),x);
%                            interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x)];
% k = linspace(0,pi,length(X));
% figure;
% plot(phi, reshape(X(1,1,:), 1,length(X)));
% hold on;
% poly_X = zeros(length(k));
% interp_X = zeros(3,3,length(k));
% for i = 1:length(k)
%     X_poly_temp = function_for_X(k(i));
%     interp_X(:,:,i) = interpolation_for_X(k(i));
%     poly_X(i) = X_poly_temp(1,1);
% end
% %plot(k,poly_X);
% plot(k,reshape(interp_X(1,1,:),1, length(k)));
% hold off;




