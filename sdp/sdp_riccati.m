function [X,time] = sdp_riccati(A,B,Q,R,t_0,T,N,M,n)
    d = 10e10;
    omega = 2*pi/T;
    time_step = T/N;
    F_zero = sdpvar(n,n,'symmetric','real');
    F_complex = sdpvar(n,n,M,'symmetric','complex');
    constraints = [];
    for i = 1:N
        H = F_zero;
        H_dot = 0;
        t = t_0+(i-1)*time_step;
        for j = 1:M
            H = H + exp(1i*j*omega*t)*F_complex(:,:,j);
            H_dot = H_dot + 1i*omega*j*exp(1i*j*omega*t)*F_complex(:,:,j);
            H = H + exp(-1i*j*omega*t)*conj(F_complex(:,:,j));
            H_dot = H_dot - 1i*omega*j*exp(-1i*j*omega*t)*conj(F_complex(:,:,j));
        end
        L = [H_dot+H*A(t)+A(t)'*H+Q(t) H*B(t);
             B(t)'*H R];
        constraints = [constraints,L >= 0,-d*eye(n) <= H <= d*eye(n)];
    end
    options = sdpsettings('solver','sdpt3','debug',1);
    objective = T*trace(F_zero);
    sol = optimize(constraints,-objective,options);
    if sol.problem == 0
            % Extract and display value
            F_0 = value(F_zero);
            complex_sol = value(F_complex);
    else
            fprintf('Hmm, something went wrong!');
            sol.info
            yalmiperror(sol.problem)
    end
    X = zeros(n,n,N);
    for i = 1:N
        X(:,:,i) = F_0;
        t = t_0+(i-1)*time_step;
        for j = 1:M
            X(:,:,i) = X(:,:,i) + exp(1i*j*omega*t)*complex_sol(:,:,j);
            X(:,:,i) = X(:,:,i) + exp(-1i*j*omega*t)*conj(complex_sol(:,:,j));
        end
    end
    time = linspace(t_0,T,N);
end

