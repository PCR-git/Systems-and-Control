
% MAE 270B
% Assignment 10
% Peter Racioppo

function [Contr, sysCL] = assign10_PeterRacioppo(Plant, Q, R, W, ny1)

A = Plant.a;
B = Plant.b;
C = Plant.c;
D = Plant.d;

ts = Plant.ts; % Sample time
if ts < 0 && t ~= -1
    disp('Error: invalid ts');
end

w = size(W,2);
B1 = B(:,1:w);
B2 = B(:,w+1:end);
C1 = C(1:ny1,:);
C2 = C(ny1+1:end,:);
D11 = D(1:ny1,1:w);
D12 = D(1:ny1,w+1:end);
D21 = D(ny1+1:end,1:w);

% ------------------
[K,~,~] = dlqr(A,B2,Q,R,0);

Rk = D21*W*D21';
K0 = Rk\D21*W*B1';
Qk = (B1 - K0'*D21)*W*(B1 - K0'*D21)';

[K1, P] = dlqr((A - K0'*C2)', C2' , Qk, Rk);
F = (K0 + K1)';

M = P*C2'/(Rk + C2*P*C2');
I = eye(size(M*C2));

% ------------------

Acontr = A-F*C2-B2*K*(I-M*C2);
Bcontr = F-B2*K*M;
Ccontr = K*(I-M*C2);
Dcontr = K*M;

Contr = ss(Acontr,Bcontr,Ccontr,Dcontr,ts);

% ------------------

Acl = [A-B2*Dcontr*C2, -B2*Ccontr;...
       Bcontr*C2     , Acontr];

Bcl = [B1-B2*Dcontr*D21, B2;...
       Bcontr*D21, zeros(size(Bcontr*D21,1),size(B2,2))];
   
Ccl = [C1-D12*Dcontr*C2, -D12*Ccontr;...
       C2, zeros(size(C2,1),size(-D12*Ccontr,2))];
   
Dcl = [D11-D12*Dcontr*D21, D12;...
       D21, zeros(size(D21,1),size(D12,2))];
   
sysCL = ss(Acl,Bcl,Ccl,Dcl,ts);

end

