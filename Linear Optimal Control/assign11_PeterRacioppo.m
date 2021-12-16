
% MAE 270B
% Homework 11
% Peter Racioppo

function [Contr, sysCL] = assign11_PeterRacioppo(Plant, Disturbance, Q, R, W, ny1)

A = Plant.a;
B = Plant.b;
C = Plant.c;
D = Plant.d;

ts = Plant.ts; % Sample time
if ts < 0 && t ~= -1
    disp('Error: invalid ts');
end

A_eta = Disturbance.a;
B_eta = Disturbance.b;
C_eta = Disturbance.c;

nw_eta = size(B_eta, 2);
nw = size(W, 2) - nw_eta;
% nu = size(R, 2);

ts=Plant.ts;

% ------------------

B1 = B(:,1:nw);
B2 = B(:,nw+1:end);
C1 = C(1:ny1,:);
C2 = C(ny1+1:end,:);
D11 = D(1:ny1,1:nw);
D12 = D(1:ny1,nw+1:end);
D21 = D(ny1+1:end,1:nw);

% ------------------
% New state space matrices

A_   = [A, B2*C_eta; zeros(size(A_eta,1),size(A,2)), A_eta];
B1_  = [B1, zeros(size(B1,1),size(B_eta,2));zeros(size(B_eta,1),size(B1,2)), B_eta];
B2_  = [B2;zeros(size(B_eta,1),size(B2,2))];
C2_  = [C2, zeros(size(C2,1),size(C_eta,2))];
D21_ = [D21, zeros(size(D21,1),nw_eta)];

% ------------------

[K,~,~] = dlqr(A,B2,Q,R,0);
K_til = [K C_eta];

Rk_ = D21_*W*D21_';
K0_ = Rk_\D21_*W*B1_';
Qk_ = (B1_ - K0_'*D21_)*W*(B1_ - K0_'*D21_)';

[K1_, P_] = dlqr((A_ - K0_'*C2_)', C2_' , Qk_, Rk_);
F_ = (K0_ + K1_)';

M_ = P_*C2_'/(Rk_ + C2_*P_*C2_');
I_ = eye(size(M_*C2_));

% ------------------

Acontr = A_- F_*C2_-B2_*K_til*(I_- M_*C2_);
Bcontr = F_- B2_*K_til*M_;
Ccontr = K_til*(I_- M_*C2_);
Dcontr = K_til*M_;

Contr = ss(Acontr,Bcontr,Ccontr,Dcontr,ts);

% ------------------
% Use old state space matrices

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
