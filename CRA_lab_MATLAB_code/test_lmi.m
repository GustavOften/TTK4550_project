setlmis([]) 
N = 400;
q = 20;
X = [];
phi = linspace(0,pi,N);
A = zeros(3,3,N);
B = zeros(3,1,N);
Q = eye(3);
R = 1;
X = [];
sX = [];
omega = 1/2;
for i = 1:q
    [X(i),~, sX] = lmivar(1,[3 1]); % variable X, full symmetric
    [A(:,:,i), B(:,:,i)] = bf.get_linearization(phi(i), bf.function_for_dphi(phi(i)));
end
for i = 1:N
   lmiterm([i,1,1,0],Q);
   for k = 1:q
       lmiterm([i,1,1,X(k)],k*omega*sin(k*omega*phi(i)),1);
       lmiterm([i,1,1,X(k)],-A(:,:,i),cos(k*omega*phi(i)),'s');
       lmiterm([i,1,1,0],-Q);
       lmiterm([i,1,2,X(k)],-cos(k*omega*phi(i)),B(:,:,i));
       lmiterm([i,2,1,X(k)],-B(:,1,i)'*cos(k*omega*phi(i)),1);
       lmiterm([i,2,2,0],-R);
   end
end

LMIs = getlmis;

%c = mat2dec(LMIs,5*eye(3));
c = zeros(6,q);
for j = 1:N
    for k = 1:q
        c(1,k) = c(1,k) + cos(k*omega*phi(j));
        c(3,k) = c(3,k) + cos(k*omega*phi(j));
        c(6,k) = c(6,k) + cos(k*omega*phi(j));
    end
end
c = -reshape(c,q*6,1);
options = [1e-5,0,0,0,0] ;
[copt1,xopt1] = mincx(LMIs,c,options);

