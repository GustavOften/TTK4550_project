function [X,time]= one_shot_generator(A,B,Q,R,t_0,T)
    % Step one: Compute monodromy Psi(t_0) = Phi(t_0+T,t_0)
    % This is done by solving the initial value problem 
    % d/dt(Phi(t,t_0) = H(t)*Phi(t,t_0);
    % H(t) = [A(t) -B(t)*R(t)^-1*B(t)^T;
    %         -Q(t) -A(t)^T]
    n = 3;
    H = @(t)([A(t) -B(t)*R^-1*transpose(B(t));
                        -Q -transpose(A(t))]);
    [Psi, phi] = matrix_initial_value_problem(H,[2*n,2*n],t_0,t_0+T,0.005,eye(2*n));
    % Step two: Compute the ordered real Schur form of Psi(t_0)
%     Psi_plot = zeros(6,6,length(phi));
%     for i = 1:length(phi)
%         temp_Psi = reshape(Psi(:,i),2*n,2*n);
%         Psi_plot(:,:,i) = temp_Psi;
%     end
%     k = 1;
%     for i = 1:6
%         for j = 1:6
%             %subplot(6,6,k)
%             hold on;
%             plot(phi,reshape(Psi_plot(i,j,:),1,length(phi)));
%             k = k+1;
%         end
%     end
%     sgtitle('Psi from first solution');
    Psi_t0_T = reshape(Psi(:,end),2*n,2*n);
    [U,tau] = schur(Psi_t0_T, 'real');
    [US, TS] = ordschur(U,tau, 'udi');
    %eig(TS(1:3,1:3));
    %US
    stable_subspace = US(1:2*n,1:n)
    %matrix_initial_value_problem(H,[2*n,2*n],t_0,t_0+T,[stable_subspace zeros(4,2)]);
    % Step three: Solve the matrix differential equation: 
    % d/dt(Y(t) = H(t)*Y(t)
    [Y, time] = matrix_initial_value_problem(H,[2*n,n],t_0,t_0+T,0.05,stable_subspace);
    length(Y)
    % Step four: Partition solution from three into [Y1(t);Y2(t)]
    % Solution of PRDE is then X(t) = Y2(t)*Y1(t)^-1;    
    X = zeros(n,n,length(Y));
    for i = 1:length(Y)
        temp_Y = reshape(Y(:,i),2*n,n);
        Y1 = temp_Y(1:n,1:n);
        Y2 = temp_Y(n+1:2*n,1:n);
        X(:,:,i) = Y2*Y1^-1;
    end
end
