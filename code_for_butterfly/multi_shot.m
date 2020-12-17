function [X, phi] =  multi_shot(A,B,Q,R,t_0, periode, N)
    n = 3;
    delta = periode/N;
    Phi = zeros(2*n,2*n,N);
    H = @(t)([A(t) -B(t)*R^-1*transpose(B(t));
                        -Q -transpose(A(t))]);
%     figure
%     k = 1;
%     time = linspace(0, periode, 100);
%     H_plot = zeros(6,6,length(time));
%     for i = 1:length(time)
%         H_plot(:,:,i) = H(time(i));
%     end
%     for i = 1:6
%         for j = 1:6
%             subplot(6,6,k)
%             hold on;
%             plot(time,reshape(H_plot(i,j,:),1,length(time)));
%             k = k+1;
%         end
%     end                
%     
%                     
    parfor i = 1:N
        [sol,~] = matrix_initial_value_problem(H,[2*n,2*n],delta*(i-1),delta*(i),delta/100,eye(2*n));
        Phi(:,:,i) = reshape(sol(:,end),2*n,2*n);
        fprintf("Number: %d",i);
    end
    %%%%%%%%% Do fast algorithm to compute X.
    %1. Creat matrix S-zT
    S_zT = zeros(2*n*N,2*n*N);
    z = 0;
    S_zT(N*2*n-(2*n-1):N*2*n,1:2*n) = -z*eye(2*n);
    for i = 1:N
        S_zT((i)*(2*n)-(2*n-1):(i)*2*n,(i)*(2*n)-(2*n-1):(i)*2*n) = Phi(:,:,i);
    end
    for i = 2:N
        S_zT((i-1)*(2*n)-(2*n-1):(i-1)*2*n,(i)*(2*n)-(2*n-1):(i)*2*n) = -eye(2*n);
    end

%     k = eig(S_zT);
%     figure
%     hold on;
%     f = linspace(0,2*pi,200);
%     scatter(cos(f),sin(f),'r','.')
%     scatter(real(k),imag(k));
    %U = zeros(4*n,4*n,N);
    for i = 2:N
        I_phi = [S_zT((i-1)*(2*n)-(2*n-1):(i-1)*2*n,(i)*(2*n)-(2*n-1):(i)*2*n);
                 S_zT((i)*(2*n)-(2*n-1):(i)*2*n,(i)*(2*n)-(2*n-1):(i)*2*n)];
        [Q, ~] = qr(I_phi);
        U = Q^-1;
        %U*S_zT((i-1)*(2*n)-(2*n-1):i*2*n,:);
        S_zT((i-1)*(2*n)-(2*n-1):i*2*n,:) = U*S_zT((i-1)*(2*n)-(2*n-1):i*2*n,:);
    end
    %S_zT(N*2*n-(2*n-1):N*2*n,1:4*n);
   	[U,tau] = schur(S_zT(N*2*n-(2*n-1):N*2*n,1:2*n), 'real');
    [US, TS] = ordschur(U,tau, 'udi');
    X = zeros(n,n,N);
    X(:,:,1) = US(n+1:2*n,1:n)*inv(US(1:n,1:n));
    Phi_11 = Phi(1:n,1:n,N);
    Phi_12 = Phi(1:n,n+1:2*n,N);
    Phi_21 = Phi(n+1:2*n,1:n,N);
    Phi_22 = Phi(n+1:2*n,n+1:2*n,N);
    X(:,:,N) = inv(X(:,:,1)*Phi_12-Phi_22)*(Phi_21-X(:,:,1)*Phi_11);
    %%%%% Using X_k = (X_k+1*Phi_12-Phi_22)^-1*(Phi_21-X_k+1*Phi_11)
    for i = N-1:-1:2
         Phi_11 = Phi(1:n,1:n,i);
         Phi_12 = Phi(1:n,n+1:2*n,i);
         Phi_21 = Phi(n+1:2*n,1:n,i);
         Phi_22 = Phi(n+1:2*n,n+1:2*n,i);
         X(:,:,i) = inv(X(:,:,i+1)*Phi_12-Phi_22)*(Phi_21-X(:,:,i+1)*Phi_11);
    end
    phi = linspace(0,periode,N);
end