
% Peter Racioppo
% Linear Systems - HW8

clear all;
close all;
clc;

% -----------------------------------------------

% QS 1
% (a.) % Finding dimenion of minimal realization

% State-space system
A = [-2,0,1,0;1,-4,1,1;0,0,-1,0;1,0,1,-3];
B = [2;2;1;2];
C = [-1,0,1,1];

% Check that the system is asymp. stable
% eig(A) % A is Hurwitz

% Lyapunov Eqs. for Inf. Horizon Obs. and Con. Grammians
% A'*Go(inf) + Go(inf)*A = -C'*C
% A*Gc(inf) + Gc(inf)*A' = -B*B'

% Solving algebraic Lyapunov Eq.
syms go11 go12 go13 go14 go21 go22 go23 go24...
     go31 go32 go33 go34 go41 go42 go43 go44...
     gc11 gc12 gc13 gc14 gc21 gc22 gc23 gc24...
     gc31 gc32 gc33 gc34 gc41 gc42 gc43 gc44
Go_s = [go11,go12,go13,go14;go21,go22,go23,go24;go31,go32,go33,go34;go41,go42,go43,go44];
Gc_s = [gc11,gc12,gc13,gc14;gc21,gc22,gc23,gc24;gc31,gc32,gc33,gc34;gc41,gc42,gc43,gc44];
eqn1 = A'*Go_s + Go_s*A == -C'*C; % Obs. Lyap. Eq.
eqn2 = A*Gc_s + Gc_s*A' == -B*B'; % Con. Lyap. Eq.
S1 = solve(eqn1);
S2 = solve(eqn2);
Go = double([S1.go11,S1.go12,S1.go13,S1.go14;S1.go21,S1.go22,S1.go23,S1.go24;S1.go31,S1.go32,S1.go33,S1.go34;S1.go41,S1.go42,S1.go43,S1.go44]);
Gc = double([S2.gc11,S2.gc12,S2.gc13,S2.gc14;S2.gc21,S2.gc22,S2.gc23,S2.gc24;S2.gc31,S2.gc32,S2.gc33,S2.gc34;S2.gc41,S2.gc42,S2.gc43,S2.gc44]);

% Observability and Controllability Matrices
Obs = [C; C*A; C*A^2; C*A^3];
Con = [B, A*B, A^2*B, A^3*B];

r1 = rank(Go);  % Rank of Obs. Grammian = 2
r2 = rank(Gc);  % Rank of Con. Grammian = 2
r3 = rank(Obs); % Rank of Obs. Matrix = 2
r4 = rank(Con); % Rank of Con. Matrix = 2
R1 = rank(Go*Gc); % Rank of Grammian Product = 1
R2 = rank(Obs*Con); % Rank of Matrix Product = 1

% ---------------------
% (b.) Determining minimal realization

% Transforming to canonical observable coordinates
% Using Grammian
% [Uo,So,Vo] = svd(Go);
% A_t = Vo'*A*Vo;
% C_t = C*Vo;

% Using Obsv. Matrix
[~,~,Vo] = svd(Obs);
A_to = Vo'*A*Vo;
B_to = Vo'*B;
C_to = C*Vo;

% Truncating unobservable states
A_tos = A_to(1:2,1:2);
B_tos = B_to(1:2);
C_tos = C_to(1:2);
 
% [Uc,Sc,Vc] = svd(Con);
% A_t = Uc'*A*Uc;
% B_t = Uc'*B;

% Transforming canonical observable system to canonical controllable
% coordinates
Con2 = [B_tos, A_tos*B_tos];
[Uc2,Sc2,Vc2] = svd(Con2);
A_t = Uc2'*A_tos*Uc2;
B_t = Uc2'*B_tos;
C_t = C_tos*Uc2;

% Truncating uncontrollable states
A_r = A_t(1,1);
B_r = B_t(1);
C_r = C_t(1);

% Minimal Realization:
% w' = -w + u
% y = w

% ---------------------
% (c.) % Computing transfer function of minimal realization

% Transfer Function
% s*F(s) = -F(s) + u/s
% y/s = F(s)
% => u = s*(s+1)*F(s)
% => u = y*(s+1)
% => y/u = 1/(s+1)

% Also:
% y/u = C*(sI-A)^-1*B
% = 1/(s+1)

% ---------------------
% (d.) Plotting frequency response of full system & min. realization

% Create state space model
Ts = 0.01; % Time step
D = 0; % D matrix is zeros
sys1 = ss(A,B,C,D,Ts);

% % Compute transfer function
% z = tf('z',Ts);
% sys2 = C/(z*eye(4) - A)*B + D;
% 
% % Alternate method to compute transfer function
% [num,den] = ss2tf(A,B,C,D);
% sys3 = tf(num,den,Ts);

% State space model for minimal realization
D_r = 0;
sys_r = ss(A_r,B_r,C_r,D_r,Ts);

% Bode Plots
% bode(sys1,'b',sys2,'r--',sys3,'g--'); % Plots are the same
figure;
bode(sys1,'b',{0.01,100});
hold on;
bode(sys_r,'r--',{0.01,100});
grid on;

% (The Bode Plots are the same)

% ---------------------
% (e.) Diagonalizing A and confirming controllable/observable modes

[V,D] = eig(A);
A_d = V\A*V;
B_d = V\B;
C_d = C*V;

% The observable states correspond to lambda = -3,-1.
% The controllable states correspond to lambda = -2,-1.

% On the other hand,
[Uo,So,Vo] = svd(Obs);
A_to = Vo'*A*Vo;
% So the observable states correspond to lambda = -3,-1.

[Uc,Sc,Vc] = svd(Con);
A_tc = Uc'*A*Uc;
% So the controllable states correspond to lambda = -2,-1.

% Thus:
% Observable, Controllable System:     lambda = -1
% Observable, Uncontrollable System:   lambda = -3
% Unobservable, Controllable System:   lambda = -2
% Unobservable, Uncontrollable System: lambda = -4

% -----------------------------------------------
% QS 2: Kalman Decomposition

rank(Con); % The rank of the controllability matrix is 2, so the first
% 2 columns are linearly independent.
M1 = Con(:,1:2); % Basis for S_c
% rank(M1); % Check rank of M1 (rank(M1) = 2).
rank(Obs); % The rank of the observability matrix is 2.
[~,~,V] = svd(Obs); % A basis for S_o_bar is given by the last
% 2 columns of an SVD of the observability matrix.
M2 = V(:,3:4); % Basis for S_o_bar
% rank(M2) % Check. rank(M2) = 2.
M = [M1 M2]; % Basis for S_c v S_o_bar 
rank(M); % Check the rank of M.

% To find basis for S_c ^ S_o_bar ,
% Need to find a vector in the range of both M1 & M2
% Find v1, v2 s.t. M1*v1 = M2*v2
% => M1*v1 - M2*v2 = 0
% Let v = [v1 -v2]^T. Then [M1 M2]*v = M*v = 0
% That is, find v in the nullspace of M.
% Can take v1 as the first elements of v,
% since v2 is guaranteed to exist.

N = null(M); % Basis for nullspace of M
v1 = N(1:2); % Taking v1 from v
T2 = M1*v1; % Basis for S_c ^ S_o bar

T1 = M1(:,1); % Choose T1 as first column of basis for controllable subspace
M3 = [T1,T2]; rank(M3); % Check that T1 and T2 are linearly independent

T4 = M2(:,2); % Choose T4 as 2nd column of basis for unobservable subspace
M4 = [T4,T2]; rank(M4); % Check that T4 and T2 are linearly independent

% Finding a vector linearly independent from T1, T2, T4
M5 = [T1 T2 T4];
rref(M5);
T3 = [0;0;0;1]; % T3 is linearly independent

T = [T1 T2 T3 T4]; % Transformation matrix
rank(T); % Check that T is full rank

% Kalman Decomposition
A_t = T\A*T;
B_t = T\B;
C_t = C*T;

% syms s
% collect(C_t*inv(s*eye(4)-A_t)*B_t,s)

% Minimal Realization
A_t_r = A_t(1,1);
B_t_r = B_t(1);
C_t_r = C_t(1);

% Computing transfer function
syms s
MinTransf = C/(s-A_t_r)*B; % 1/(s + 1)
% The minimal realization has the same transfer function
% as the result from Problem 1.

% -----------------------------------------------
% QS 3: A Perturbed Version of (1)
% (a.) Plotting freq. response

Ap = [-2,0.1,1,0;1,-4,1,1;0,-0.1,-1,0;1,0,1.1,-3];
Bp = [2.1;2.1;1.1;1.9];
Cp = (1/1.23)*[-0.9,0.1,0.9,1.1];
Dp = 0;
Tsp = 0.01;

% Create state space model
sysp1 = ss(Ap,Bp,Cp,Dp,Tsp);

% Bode Plot
figure;
bode(sysp1,'b',{0.01,100});
grid on;
hold on;
bode(sys1,'r--',{0.01,100});

% ---------------------
% (b.) Confirming A is asymp. stable and minimal

Qp = eye(4); % Qp is positive definite
Pp = lyap(Ap',Qp); % Solving for P in Lyap. Eq.
% eig(Pp)
% The eigenvalues of Pp are greater than zero,
% so Pp is positive definite.
% Thus, the Lyapunov condition is satisfied,
% so the system is asymptotically stable.

% Alternately:
% eig(Ap)
% The eigenvalues of Ap are strictly negative,
% so Ap is Hurwitz,
% and the system is asymptotically stable.

% Computing observability and controllability Grammians
Go_p = lyap(Ap',Cp'*Cp);
Gc_p = lyap(Ap,Bp*Bp');

eigs_o = eig(Go_p);
eigs_c = eig(Gc_p);
% eigs_o(1) % 1.8975e-06
% eigs_c(1) % 7.5600e-07

% The Grammian eigenvalues are small, but nonzero,
% so the Grammians are non-singular,
% and the system is both observable and controllable,
% i.e. minimal.

% ---------------------
% (c.) Computing Hankel singular values

H = Go_p*Gc_p;
SvdH = sqrt(eig(H)); % Hankel singular values; 3 are small

% The absolute value of the difference of the frequency responses
% of the original system and the truncated system is no greater than:
% 2*(sum of truncated Hankel singular values)

ErrorBound = 2*sum(SvdH(2:4)); % 0.0031

% Computing transfer function
syms s
Ht = Cp/(s*eye(4) - Ap)*Bp;
DC_Gain = double(subs(Ht,s,0)); % DC Gain = 1.0002
% (The transfer function is approx. 1/(s+1)).

% Error bound on percent difference of transfer functions
PercentDiff = 100*ErrorBound/DC_Gain; % 0.3080

% Error of single-state system is about 0.31%.

% ---------------------
% (d.) Convert to balanced coordinates and truncate to 1 state

S_p = diag(SvdH); % Sigma matrix
GchGoGch_p = sqrtm(Gc_p)*Go_p*sqrtm(Gc_p);
[Up,Sp,Vp] = svd(GchGoGch_p);
% Up'*Up % Check Up is unitary

Sigma_2 = Up'*sqrtm(Gc_p)*Go_p*sqrtm(Gc_p)*Up; % Sigma squared
Sigma = real(sqrt(Sigma_2)); % Sigma matrix
Tp = sqrtm(Gc_p)*Up/(sqrtm(Sigma));  % T matrix

% % Check that the transformations diagonalize Go_p and Gc_p:
% Go_b = Tp'*Go_p*Tp
% Gc_b = inv(Tp)*Gc_p*inv(Tp')

% Pg 12 of Notes 15:
% A --> inv(T)*A*T
% B --> inv(T)*B
% C --> C*T

% The system in balanced coordinates:
Ap_b = Tp\Ap*Tp;
Bp_b = Tp\Bp;
Cp_b = Cp*Tp;

% Single-state system:
Ap_b1 = Ap_b(1,1);
Bp_b1 = Bp_b(1);
Cp_b1 = Cp_b(1);

Ts = 0.01; % Time step
sysp_full = ss(Ap_b,Bp_b,Cp_b,0,Ts); % Full state space model
% sysp_full2 = Cp_b/(s*eye(4) - Ap_b)*Bp_b
sysp_r = ss(Ap_b1,Bp_b1,Cp_b1,0,Ts); % Reduced-state space model

% Bode Plots
figure;
bode(sysp_full,'b',{0.01,100});
hold on;
bode(sysp_r,'r--',{0.01,100});
grid on;

% The Bode plots of the full and reduced-order systems are
% indistinguishable.

% -----------------------------------------------
% QS 4: Linearized EoM of a Satellite
% (a.) Controllability

% State-space system
syms w0 u1 u2
A4 = [0,1,0,0; 3*w0^2,0,0,2*w0; 0,0,0,1; 0,-2*w0,0,0];
B4 = [0,0;1,0;0,0;0,1];
C4 = [1,0,0,0;0,0,1,0];
D4 = 0;

Con4 = [B4, A4*B4, A4^2*B4, A4^3*B4]; % Contr. Matrix
% rref(Con4) % Con4 is full rank
rank(Con4); % Con4 has rank 4
% The system is controllable.

% ---------------------
% (b.) Controllability if u2 = 0

% Method 1:
B4 = [0,0;1,0;0,0;0,0]; % Reset B4
Con4 = [B4, A4*B4, A4^2*B4, A4^3*B4]; % Contr. Matrix
rank(Con4); % Con4_1 has rank 3
% The system is not controllable with only u1.

% % Method 2 (less computation):
% B4_s = [0;u1;0;u2];
% Con4_s = [B4_s, A4*B4_s, A4^2*B4_s, A4^3*B4_s];
% Con4_1 = subs(Con4_s,u2,0)
% rank(Con4_1) % Con4_1 has rank 3

% ---------------------
% (c.) Controllability if u1 = 0

% Method 1:
B4 = [0,0;0,0;0,0;1,0]; % Reset B4
Con4 = [B4, A4*B4, A4^2*B4, A4^3*B4]; % Contr. Matrix
rank(Con4); % Con4_1 has rank 4
% The system is controllable with only u2.

% % Method 2 (less computation):
% Con4_2 = subs(Con4_s,u1,0);
% rank(Con4_2) % Con4_2 has rank 4

% ---------------------
% (d.) Observability

B4 = [0,0;1,0;0,0;0,1]; % Reset B4
Obs4 = [C4; C4*A4; C4*A4^2; C4*A4^3]; % Obs. Matrix
rank(Obs4); % Obs 4 is full rank (4)
% The system is observable.

% ---------------------
% (e.) Observability if y2 = 0

C4 = [1,0,0,0;0,0,0,0]; % Reset C4
Obs4 = [C4; C4*A4; C4*A4^2; C4*A4^3]; % Obs. Matrix
rank(Obs4); % Obs 4 has rank 3
% The system is unobservable.

% ---------------------
% (f.)  Observability if y1 = 0

C4 = [0,0,0,0;0,0,1,0]; % Reset C4
Obs4 = [C4; C4*A4; C4*A4^2; C4*A4^3]; % Obs. Matrix
rank(Obs4); % Obs 4 has rank 4
% The system is observable.

% ---------------------
% (g.) Comparing transfer function poles and matrix ranks

% Transfer function def: C*inv(s*I-A)*B

% Resetting Matrices
B4 = [0,0;1,0;0,0;0,1];
C4 = [1,0,0,0;0,0,1,0];

TM = C4/Mid*B4; % Transfer matrix
% TM =
% [           1/(s^2 + w0^2),           (2*w0)/(s*(s^2 + w0^2))]
% [ -(2*w0)/(s*(s^2 + w0^2)), (s^2 - 3*w0^2)/(s^2*(s^2 + w0^2))]

% TM = [y1/u1, y1/u2; y2/u1, y2/u2];
% when one or the other input/output is zero.

% Number of poles vs rank of obs/con matrices:

% When u2 = 0, y2 = 0, the transfer function has 2 poles.
% In this case, the obs. and cont. matrices both have rank 3.
% But the intersection of the obs. and con. subspaces must be rank 2.
B4 = [0,0;1,0;0,0;0,0]; % Reset B4
Con4 = [B4, A4*B4, A4^2*B4, A4^3*B4]; % Con. matrix
C4 = [1,0,0,0;0,0,0,0]; % Reset C4
Obs4 = [C4; C4*A4; C4*A4^2; C4*A4^3]; % Obs. matrix
rank(Obs4*Con4); % Rank of their matrix product = 2

% When u1 = 0, y2 = 0, the transfer function has 3 poles.
% The obs. matrix has rank = 3, and the cont. matrix has rank = 4.

% When u2 = 0, y1 = 0, the transfer function has 3 poles.
% The obs. matrix has rank = 4, and the cont. matrix has rank = 3.

% When u1 = 0, y1 = 0, the transfer function has 4 poles.
% In this case, the obs. and cont. matrices both have rank 4.

% To summarize, the number of poles of each transfer function is equal
% to the rank of the intersection of the corresponding observability
% and controllability subspaces.

% -----------------------------------------------

