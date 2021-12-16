
% MAE 270B
% Homework

% ------------------------------

% Assignment 1

function W = assign1_PeterRacioppo(A,B,C,D,K)

% (Both options work)
% Option 1
W = lyap((A-B*K)', (C-D*K)'*(C-D*K));

% Option 2
% sys = ss(A-B*K,B,C,D);
% W = gram(sys,'o');

end

% ------------------------------

% Assignment 2

clear;
disp('Assignment 2: ');
    
A = [0,     1,  0,     0;...
    -1, -0.02,  0,     0;...
     0,     0,  0,     1;...
     0,     0, -4, -0.02];
B1 = [0;1;0;1];
B2 = [0;1;0;-1];
C2 = [1,0,1,0];

Z = zeros(4,4);
Q1 = Z; Q1(1,1) = 1; Q1(2,2) = 1;
Q2 = Z; Q2(3,3) = 1; Q2(4,4) = 1;
Q3 = eye(4);
R = 1;

% [K1,P1,~] = lqr(A,B1,Q1,R);
% [K2,P2,~] = lqr(A,B1,Q2,R);
% [K3,P3,~] = lqr(A,B1,Q3,R);

[K1,P1] = lqr_PCR_1(A,B1,Q1,R)
[K2,P2] = lqr_PCR_1(A,B1,Q2,R)
[K3,P3] = lqr_PCR_1(A,B1,Q3,R)

sys1 = ss(A,[B1 B2],C2,0);
sys1_1 = ss(A-B1*K1,[B1 B2],C2,0);
sys1_2 = ss(A-B1*K2,[B1 B2],C2,0);
sys1_3 = ss(A-B1*K3,[B1 B2],C2,0);

figure;
opts = bodeoptions('cstprefs');
opts.PhaseWrapping = 'on';
opts.Grid = 'on';

bode(sys1,'b',opts);
hold on;
bode(sys1_3,'r',opts);

figure;
bode(sys1_1,'b',opts);
hold on;
bode(sys1_2,'r',opts);

% ------------------------------

% Functions

% Solve LQR using care function
function [K,P] = lqr_PCR_1(A,B,Q,R)
    
    [P,~,~] = care(A,B,Q);
    K = inv(R)*B'*P;

end

% Solves LQR manually
function [K,P] = lqr_PCR_2(A,B,Q,R)

    M = B*inv(R)*B';

    d = length(A);
    P = sym('P', [d d]);
    P = triu(P) + transpose(triu(P)) - diag(diag(triu(P)));
    Pvec1 = reshape(triu(P),[1,length(P)^2]);
    Pvec = transpose(nonzeros(Pvec1));
    
    Z = zeros(size(A));
    eqn1 = A'*P + P*A - P*M*P + Q == Z;
    soln = solve(eqn1,Pvec);
    
    if d == 2 || d== 3 || d== 4
        p11s = soln.P1_1; p12s = soln.P1_2; p22s = soln.P2_2;
    end
    
    if d == 3 || d== 4
        p13s = soln.P1_3; p23s = soln.P2_3; p33s = soln.P3_3;
    end
    
    if d == 4
        p14s = soln.P1_4; p24s = soln.P2_4;
        p34s = soln.P3_4; p44s = soln.P4_4;
    end

    for i = 1:length(p11s)

        if d == 2
            P = [p11s(i),p12s(i);...
                 p12s(i),p22s(i)];
        elseif d ==3
            P = [p11s(i),p12s(i),p13s(i);...
                 p12s(i),p22s(i),p23s(i);...
                 p13s(i),p23s(i),p33s(i)];
        elseif d ==4
            P = [p11s(i),p12s(i),p13s(i),p14s(i);...
                 p12s(i),p22s(i),p23s(i),p24s(i);...
                 p13s(i),p23s(i),p33s(i),p34s(i);...
                 p14s(i),p24s(i),p34s(i),p44s(i)];
        end
        
        Sgn = sign(real(eig(A-M*P)));

        if Sgn == -ones(1,length(Sgn))
            % real(eig(X))
            break
        end
    end
    
    K = double(inv(R)*B'*P);
    P = double(P);

end

% ------------------------------

% ------------------------------
% Assignment 3

clear;
disp('Assignment 3: ');

A = [0,1;-4,0];
B = [0;1];
C = [1,0];
Q = C'*C;
R = 1;

% Solves LQR manually
[K1,P1] = lqr_PCR_2(A,B,Q,R)

% Solves LQR using care function
[K2,P2] = lqr_PCR_1(A,B,Q,R)

% Solves LQR using lqr function
[K3,P3] = lqr(A,B,Q,R)

% P3 = lyap(A',Q-P*(B*B')*P)

% R0 = 0;
% Q0 = 0*Q;
% [K0,W0] = lqr_PCR_1(A,B,Q,R0)
% [K1,W1] = lqr_PCR_1(A,B,Q0,R)

% W0 + W1 = P
% Let W0 = W1 = P/2

W0 = P1/2
W1 = P1/2

% ------------------------------

% Functions

% Solve LQR using care function
function [K,P] = lqr_PCR_1(A,B,Q,R)
    
    [P,~,~] = care(A,B,Q);
    K = inv(R)*B'*P;

end

% Solves LQR manually
function [K,P] = lqr_PCR_2(A,B,Q,R)

    M = B*inv(R)*B';

    d = length(A);
    P = sym('P', [d d]);
    P = triu(P) + transpose(triu(P)) - diag(diag(triu(P)));
    Pvec1 = reshape(triu(P),[1,length(P)^2]);
    Pvec = transpose(nonzeros(Pvec1));
    
    Z = zeros(size(A));
    eqn1 = A'*P + P*A - P*M*P + Q == Z;
    soln = solve(eqn1,Pvec);
    
    if d == 2 || d== 3 || d== 4
        p11s = soln.P1_1; p12s = soln.P1_2; p22s = soln.P2_2;
    end
    
    if d == 3 || d== 4
        p13s = soln.P1_3; p23s = soln.P2_3; p33s = soln.P3_3;
    end
    
    if d == 4
        p14s = soln.P1_4; p24s = soln.P2_4;
        p34s = soln.P3_4; p44s = soln.P4_4;
    end

    for i = 1:length(p11s)

        if d == 2
            P = [p11s(i),p12s(i);...
                 p12s(i),p22s(i)];
        elseif d ==3
            P = [p11s(i),p12s(i),p13s(i);...
                 p12s(i),p22s(i),p23s(i);...
                 p13s(i),p23s(i),p33s(i)];
        elseif d ==4
            P = [p11s(i),p12s(i),p13s(i),p14s(i);...
                 p12s(i),p22s(i),p23s(i),p24s(i);...
                 p13s(i),p23s(i),p33s(i),p34s(i);...
                 p14s(i),p24s(i),p34s(i),p44s(i)];
        end
        
        Sgn = sign(real(eig(A-M*P)));

        if Sgn == -ones(1,length(Sgn))
            % real(eig(X))
            break
        end
    end
    
    K = double(inv(R)*B'*P);
    P = double(P);

end

% ------------------------------

% ------------------------------
% Assignment 4

clear;
disp('Assignment 4: ');

A = [0,1;9,0];
B = [0;1];
Q = 0*A;
R = 1;

[K,P] = lqr_PCR_1(A,B,Q,R)

H = [A,-B*inv(R)*B';0*A,-A'];

eigA_BK = eig(A-B*K);
eigH1 = eig(H);
eigH = eigH1(eigH1 < 0);

eigA1 = [eig(A);-eig(A')];
eigA = eigA1(eigA1 < 0);

% ------------------------------

% Functions

% Solve LQR using care function
function [K,P] = lqr_PCR_1(A,B,Q,R)
    
    [P,~,~] = care(A,B,Q);
    K = inv(R)*B'*P;

end

% ------------------------------

% Homework 5

% Problem definition of the state space model
% x1' = x3
% x2' = x4
% x3' = u
% x4' = (6*9.8)*x2 + 6*u
% X' = [0,0,1,0; 0,0,0,1; 0,0,0,0; 0,6*9.8,0,0]*X + [0;0;1;6]*u;

A = [0,0,1,0; 0,0,0,1; 0,0,0,0; 0,6*9.8,0,0];
B = [0;0;1;6];

q1 = 1;
q2 = 100;

% q1 = 4;
% q2 = 1000;

R = 1;
Q = zeros(4,4);
Q(1,1) = q1; Q(2,2) = q2;
[K,~,~] = lqr(A,B,Q,R);

C = eye(4);
D = zeros(4, 1);
sys1 = ss(A-B*K, B, C, D);

T = [0:0.01:10]';
v = zeros(size(T,1),1);
x0 = [1 0 0 0]';
x = lsim(sys1,v,T,x0);
u = -x*K';

idx = find(T==4);
xa4 = x(idx:end,1);
x1max = max(xa4) % Max value of x1 after 4 seconds
x2max = max(x(:,2)) % Max value of x2 at any time

figure;
subplot(211);
plot(T,x(:,1:2));
% xlim([4,10])
grid, xlabel('time (sec)');
legend('horizontal position of mass center (m)','angle from vertical (rad)');
subplot(212);
plot(T,u);
% xlim([4,10])
grid, xlabel('time (sec)');
legend('control command');

% ------------------------------

% Assignment 9

% Takes a discrete-time plant of the form
% x(t + 1) = Ax(t) + B1w(t) + B2u(t)
% y(t) = C2x(t) + D21w(t) + D22u(t)
% and weighting matrix W in the performance index for the Kalman predictor
% and filter. Builds the one-step Kalman predictor and the Kalman filter
% for estimating x, along with the corresponding error systems,
% all in state-space form.

function [sysKP, sysKP_err, sysKF, sysKF_err] = assign9_PeterRacioppo(Plant, W)

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
C1 = C(1:ny1,:);
C2 = C(ny1:end,:);
D21 = D(:,1:w);

R = D21*W*D21';
K0 = R\D21*W*B1';
Q = (B1 - K0'*D21)*W*(B1 - K0'*D21)';

[K1, P] = dlqr((A - K0'*C2)', C2' , Q, R);
F = (K0 + K1)';

M = P*C2'/(R + C2*P*C2');
I = eye(size(M*C2));
% Z = I*0;

if size(Plant.b, 2) ~= w
    % There exist u, B2, & D22
    B2 = B(:,w+1:end);
    D22 = D(:,w+1:end);
    
    AF = A - F*C2;
    BF = [F, (B2-F*D22)];
    CF = I - M*C2;
    DF = [M, -M*D22];
else
    AF = A - F*C2;
    BF = F;
    CF = I - M*C2;
    DF = M;
end

Be = F*D21 - B1;
De = M*D21;

% In standard state-space notation, the Kalman filter has the state-space
% realization (AF,BF,CF,DF).
% The predictor with output x_hat(t|t + 1) has the state-space
% realization (AF,BF,I,0).
% The predictor with output x_hat(t + 1|t) has the
% realization (AF,BF,AF,BF).

% The error system for the Kalman filter has the state-space
% realization (AF,B,CF,D).
% For the Kalman predictor with output x_hat(t|t + 1), the corresponding
% error system has output x(t|t + 1) and state-space
% realization (AF,B,I,0).
% For the Kalman predictor with output x_hat(t + 1|t), the corresponding
% error system has output x(t + 1|t) and state-space
% realization (AF,B,AF,B).

sysKP = ss(AF,BF,AF,BF,ts);
sysKP_err = ss(AF,Be,AF,Be,ts);
sysKF = ss(AF,BF,CF,DF,ts);
sysKF_err = ss(AF,Be,CF,De,ts);

end

% ------------------------------

% Assignment 10

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

% ------------------------------

% Homework 11

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

% ------------------------------