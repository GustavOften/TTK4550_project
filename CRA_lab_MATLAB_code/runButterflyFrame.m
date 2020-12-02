clear all;
%% Butterfly Frame simulation script 
%  This code contains the some basic building blocks for simulating the
%  Butterfly frame system. See also the ButterflyFrame class for how it works. 
%  Change what you will, and please report any bugs you find to christian.f.satre@ntnu.no
%% Define the object bf of class ButterflyFrame
bf = ButterflyFrame;
%% Initial conditions and simulation time
theta0  = 0;                    % Initial angle of the frame
Dtheta0 = 0;                    % Initial angular velocity of the frame
t0      = 0;                    % Intial time
tEnd    = 10;                   % Final time (simulation runs from t0 to tEnd)
x0      = [theta0;Dtheta0];     % Initial conditions 
%% Controller and reference options
% Reference
bf.ref = 2;    % 1: Step response; 2: Sinusoid
% Controller    
bf.cont = 1;   % 1: PD ; 2: Sliding mode; 3: LQR
%----- PD gains ----%
bf.kp = 0.5;
bf.kd = 0.1;
%---- Sliding mode gain-----%
bf.SMCt = 2;         % Sliding mode type. 1: sign(s); 2: sat(s/eps)      
bf.ks1  = 0.1;
bf.ks2  = 0.1;
bf.eps  = 0.1;
%---- LQR ----%
A=[0,1;0,0];B=[0;bf.K/bf.J];            % State matrices: \dot e= Ae+Bu
Q = [0.1, 0.0 ; 0.0, 0.1];              % Weigth for the states
S = 1;                                  % Weigth for the actuator input
[Klqr,~,~] = lqr(A,B,Q,S);bf.Klqr=Klqr; % Calculate the gain matrix
%% ---- Additional controller settings ------%
% Add friction compensation (true/false)
bf.fricComp = 0;
% Add feed-forward of angular acceleration reference (true/false)
bf.ffwd     = 0;
% Set power of noise on measurements for simulink simulation
noisepw     = 0.0;
% Give a specific sample time of the controller (0.001 is default)
sampleTime  = 0.001;
%% Simulate the system
tic
% [t,x,u] = bf.simEOM(t0,tEnd,x0);     % Simulate using MATLAB
sim simulinkFrame                      % Simulate using Simulink 
[~,r,Dr]=bf.get_u_r_Dr(t,x);
toc
%% Animate the system?
animate   = 1;        % Set to 1 in order to run the animation
plotspeed = 1.0;      % Speed of plotting the animation (not very accurate)
%% Make plots
% Coordinates and actutator output vs time
figure(1);clf(1);
hold on;
subplot(4,2,1);
plot(t,r-x(:,1));
xlim([0 3]);
ylabel('error $e$ [rad]','Interpreter','latex');
subplot(4,2,2);
plot(t,r-x(:,1));
xlim([3 10]);
ylabel('error $e$ [rad]','Interpreter','latex');
subplot(4,1,2);
hold on;
plot(t,x(:,1));
plot(t,r);
ylabel(' [rad]','Interpreter','latex');
legend({'$\theta$','$r(t)$'},'Location','Best','Interpreter','latex');
title('Coordinates and actuator output vs time');
subplot(4,1,3);
hold on;
plot(t,x(:,2));
plot(t,Dr);
legend({'$\dot \theta$','$\dot r(t)$'},'Location','Best','Interpreter','latex');
ylabel('[rad/s]','Interpreter','latex');
subplot(4,1,4);
plot(t,u);
ylabel('$u$ [Nm]','Interpreter','latex');
xlabel('$t$ [s]','Interpreter','latex');
%Phase portraits
figure(3);clf(3);
hold on;
plot(x(:,1),x(:,2));
tt = linspace(t0,tEnd,100);
rt = bf.a*sin(tt);Drt=bf.a*bf.w*cos(tt);
plot(rt,Drt,'g--');
ylabel('$\dot{\theta}$ [rad/s]','Interpreter','latex');
xlabel('$\theta$ [rad]','Interpreter','latex');
title('Phase portrait');
legend({'Actual','Nominal'},'Interpreter','latex','Location','Best');
hold off;
%% Maximum value and average of the error over last n steps
n  = 50;
emax= max(abs(r(end-n+1:end)-x(end-n+1:end,1)));
eavg= sum(abs(r(end-n+1:end)-x(end-n+1:end,1)))/n;
disp(strcat('Maximum value of the error over the last',{' '}, num2str(n),' steps: ',{' '}, num2str(emax)));
disp(strcat('Average value of the error over the last',{' '}, num2str(n),' steps: ',{' '}, num2str(eavg)));
%% Animation
if animate
    bf.animate(t,x,plotspeed);
end