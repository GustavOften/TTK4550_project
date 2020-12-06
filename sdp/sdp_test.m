syms t 
f = 1/(1.5+cos(2*pi*t))*[cos(2*pi*t) sin(2*pi*t);-sin(2*pi*t) cos(2*pi*t)];
P = @(t) 1/(1.5+cos(2*pi*t))*[cos(2*pi*t) sin(2*pi*t);-sin(2*pi*t) cos(2*pi*t)];
dP = matlabFunction(diff(f,t));
A_c = [4 3;-4.5 -3.5];
B_c = [1;-1];
R_c = 1;
Q_c = [10 6;6 4];
A = @(t) inv(P(t))*A_c*P(t)-inv(P(t))*dP(t);
B = @(t) inv(P(t))*B_c;
Q = @(t) P(t)*Q_c*P(t)';
X = sdp_riccati(A,B,Q,R_c,0,1,100,50,2);

time = linspace(0,1,100);
figure
k = 1;
for i = 1:2
    for j = 1:2
        subplot(2,2,k)
        hold on;
        plot(time,reshape(X(i,j,:),1,100));
        k = k+1;
    end
end
sgtitle('Riccati sol X'); 