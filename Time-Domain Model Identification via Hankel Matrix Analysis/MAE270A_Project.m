% MAE 270A - Linear Systems
% Final Project
% Peter Racioppo
% 103953689

%% ------------------------------------------------
% Initializations

clear;
close all;

% Loading Data
load u1_impulse.mat
y11 = u1_impulse.Y(3).Data; y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data;  %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0);  %%% find index where pulse occurs

load u2_impulse.mat
y12 = u2_impulse.Y(3).Data; y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));

%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1); y12 = y12/max(u2);
y21 = y21/max(u1); y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);
ts = 1/40;  %%%% sample period
N = length(y11);  %%%% length of data sets
t = [0:N-1]*ts - 1;

% There are:
% m outputs, q inputs
% ns states, n samples

%% ------------------------------------------------

% Task #1: Model Identification
% (1)

n = 100; % Size of the Hankel Matrix
H = zeros(2*n,2*n); % Initializing Hankel Matrix
Ht = zeros(2*n,2*n); % Initializing Hankel_tilda Matrix

d_t = find(u1); % Find the index of the first nonzero element
hk = zeros(2,2,2*n); % Initializing hk
% Computing hk, the kth impulse response
for k = 1:1:2*n-1
    hk(:,:,k) = [y11(k+d_t),y12(k+d_t);y21(k+d_t),y22(k+d_t)];
end

% Constructing the Hankel matrices
for i = 1:1:n
    for j = 1:1:n
        H(2*i-1:2*i,2*j-1:2*j) = hk(:,:,i+j-1); % Hankel Matrix
        Ht(2*i-1:2*i,2*j-1:2*j) = hk(:,:,i+j); % Tilda Hankel Matrix
    end
end

[U,S,V] = svd(H);         % Computing an SVD of H
[Ut,St,Vt] = svd(Ht);     % Computing an SVD of Ht

% Check: Plot the singular values of the Hankel matrix
tS = 1:1:length(diag(S));
semilogy(tS,diag(S),'bo','MarkerFaceColor','b'); % Log-log plot
xlim([0 40]);
title('Hankel Matrix Singular Values','FontSize',25,'Interpreter','Latex');
xlabel('x','FontSize',18,'Interpreter','Latex');
ylabel('y','FontSize',18,'Interpreter','Latex');
grid on;

dim = [6,7,10,20,8]; % State dimensions of the systems

% 3D matrices to hold the A, B, and C matrices for each system
BigA = zeros(dim(4),dim(4),4);
BigB = zeros(dim(4),2,4);
BigC = zeros(2,dim(4),4);

i = 1;
for n = dim
    U1 = U(:,1:n);
    S1 = S(1:n,1:n);
    V1 = V(:,1:n);

    On = U1;     % Observability Matrix
    Cn = S1*V1'; % Controllability Matrix

    Oi = pinv(On); % Pseudoinverse of Obs. Mat.
    Ci = pinv(Cn); % Pseudoinverse of Con. Mat.

    A = Oi*Ht*Ci;  % A Matrix
    B = Cn(:,1:2); % B Matrix
    C = On(1:2,:); % C Matrix
    
    BigA(1:n,1:n,i) = A; % 3D matrix to hold A matrices for each system
    BigB(1:n,1:2,i) = B; % 3D matrix to hold B matrices for each system
    BigC(1:2,1:n,i) = C; % 3D matrix to hold C matrices for each system
    
    MaxAbsEig = max(abs(eig(A))); % 0.9526, 0.9144, 0.9141, 0.9977
    % Absolute values of eigenvalues are all less than 1

    i = i+1;
end

% Pulling As out of BigA
A1 = BigA(1:dim(1),1:dim(1),1);A2 = BigA(1:dim(2),1:dim(2),2);
A3 = BigA(1:dim(3),1:dim(3),3); A4 = BigA(1:dim(4),1:dim(4),4);

% Pulling Bs out of BigB
B1 = BigB(1:dim(1),1:2,1); B2 = BigB(1:dim(2),1:2,2);
B3 = BigB(1:dim(3),1:2,3); B4 = BigB(1:dim(4),1:2,4);

% Pulling Cs out of BigC
C1 = BigC(1:2,1:dim(1),1); C2 = BigC(1:2,1:dim(2),2);
C3 = BigC(1:2,1:dim(3),3); C4 = BigC(1:2,1:dim(4),4);

% D matrices, zero matrices in this case
D1 = zeros(2,2); D2 = zeros(2,2);
D3 = zeros(2,2); D4 = zeros(2,2);

% ----------------------
% (2)

% C*A^k-1*B

% hvec(1:2,1:2,1:401,i) holds a 3D matrix of 2x2 impulse response
% matrices for each timestep.
hvec(1:2,1:2,1:401,1) = f_ImpulseResponse(A1,B1,C1,t);
hvec(1:2,1:2,1:401,2) = f_ImpulseResponse(A2,B2,C2,t); 
hvec(1:2,1:2,1:401,3) = f_ImpulseResponse(A3,B3,C3,t); 
hvec(1:2,1:2,1:401,4) = f_ImpulseResponse(A4,B4,C4,t);

% Plotting each of the four elements of the impulse response matrices,
% for each of the four systems
f_PlotImpulse(y11,y12,y21,y22,t,d_t,hvec,1);
f_PlotImpulse(y11,y12,y21,y22,t,d_t,hvec,2);
f_PlotImpulse(y11,y12,y21,y22,t,d_t,hvec,3);
f_PlotImpulse(y11,y12,y21,y22,t,d_t,hvec,4);

% Plotting all four models together
figure;
subplot(2,2,1); hold on; grid on;
plot(t,y11,'r*','LineWidth',2);
f_PlotImpulseij(t,d_t,hvec,1,1);
xlim([-0.5 1.5]);
subplot(2,2,2); hold on; grid on;
plot(t,y21,'r*','LineWidth',2);
f_PlotImpulseij(t,d_t,hvec,1,2);
xlim([-0.5 1.5]);
subplot(2,2,3); hold on; grid on;
plot(t,y12,'r*','LineWidth',2);
f_PlotImpulseij(t,d_t,hvec,2,1);
xlim([-0.5 1.5]);
subplot(2,2,4); hold on; grid on;
plot(t,y22,'r*','LineWidth',2);
f_PlotImpulseij(t,d_t,hvec,2,2);
xlim([-0.5 1.5]);

% ----------------------
% (3) and (4)

% Estimating the frequency response for each channel
y11f = fft(y11(d_t:end-1))./fft(u1(d_t:end-1));
y12f = fft(y12(d_t:end-1))./fft(u2(d_t:end-1));
y21f = fft(y21(d_t:end-1))./fft(u1(d_t:end-1));
y22f = fft(y22(d_t:end-1))./fft(u2(d_t:end-1));

% Computing magnitudes and phases of the estimated freq. responses
me11 = abs(y11f); pe11 = angle(y11f);
me12 = abs(y12f); pe12 = angle(y12f);
me21 = abs(y21f); pe21 = angle(y21f);
me22 = abs(y22f); pe22 = angle(y22f);

N = length(y11f);
ts = 1/40; % Sampling time
w = (0:N-1)/(ts*N);

% Computing the magnitudes and phases of the frequency responses
% of each of the four models
[mag1,phase1,~] = f_FreqResponse(A1,B1,C1,D1,w,ts);
[mag2,phase2,~] = f_FreqResponse(A2,B2,C2,D2,w,ts);
[mag3,phase3,~] = f_FreqResponse(A3,B3,C3,D3,w,ts);
[mag4,phase4,~] = f_FreqResponse(A4,B4,C4,D4,w,ts);

% Picking out each input-output channel
[m1_11,p1_11] = f_ReshapeFreq(mag1,phase1,w,1,1);
[m1_12,p1_12] = f_ReshapeFreq(mag1,phase1,w,1,2);
[m1_21,p1_21] = f_ReshapeFreq(mag1,phase1,w,2,1);
[m1_22,p1_22] = f_ReshapeFreq(mag1,phase1,w,2,2);

[m2_11,p2_11] = f_ReshapeFreq(mag2,phase2,w,1,1);
[m2_12,p2_12] = f_ReshapeFreq(mag2,phase2,w,1,2);
[m2_21,p2_21] = f_ReshapeFreq(mag2,phase2,w,2,1);
[m2_22,p2_22] = f_ReshapeFreq(mag2,phase2,w,2,2);

[m3_11,p3_11] = f_ReshapeFreq(mag3,phase3,w,1,1);
[m3_12,p3_12] = f_ReshapeFreq(mag3,phase3,w,1,2);
[m3_21,p3_21] = f_ReshapeFreq(mag3,phase3,w,2,1);
[m3_22,p3_22] = f_ReshapeFreq(mag3,phase3,w,2,2);

[m4_11,p4_11] = f_ReshapeFreq(mag4,phase4,w,1,1);
[m4_12,p4_12] = f_ReshapeFreq(mag4,phase4,w,1,2);
[m4_21,p4_21] = f_ReshapeFreq(mag4,phase4,w,2,1);
[m4_22,p4_22] = f_ReshapeFreq(mag4,phase4,w,2,2);

% Plotting the magnitudes of the frequency responses
figure;
subplot(221);
f_PlotFreqResponse(m1_11,m2_11,m3_11,m4_11,me11,w);
subplot(222);
f_PlotFreqResponse(m1_12,m2_12,m3_12,m4_12,me12,w);
subplot(223);
f_PlotFreqResponse(m1_21,m2_21,m3_21,m4_21,me21,w);
subplot(224);
f_PlotFreqResponse(m1_22,m2_22,m3_22,m4_22,me22,w);

% Plotting the phases of the frequency responses
figure;
subplot(221);
f_PlotFreqResponse(p1_11,p2_11,p3_11,p4_11,pe11,w);
subplot(222);
f_PlotFreqResponse(p1_12,p2_12,p3_12,p4_12,pe12,w);
subplot(223);
f_PlotFreqResponse(p1_21,p2_21,p3_21,p4_21,pe21,w);
subplot(224);
f_PlotFreqResponse(p1_22,p2_22,p3_22,p4_22,pe22,w);

%% ------------------------------------------------
% Task #2: Transmission zeros of the MIMO model and zeros of each channel
% (1) Show that there are 5 finite transmission zeros

% [V,D] = eig(A,B) and [V,D] = eig(A,B,algorithm) return V as a matrix
% whose columns are generalized eigenvectors that satisfy A*V = B*V*D.

% A2, B2, C2, D2 are matrices for 7-state system
Aa2 = [A2,B2;-C2,-D2]; % Aa2: Matrix on LHS of generalized eigenvec/val eq.
% Bb2: Matrix on LHS of generalized eigenvec/val eq.
Bb2 = [eye(size(A2)),zeros(size(B2));zeros(size(C2)),zeros(size(D2))];

% This function calculates the transfer zeros of a system,
% ordered from largest to smallest magnitude.
TransZ = f_TransZeros(A2,B2,C2,D2);
% The 5 finite transfer zeros of the 7-state model are:
% -2.6310 + 0.0000i     0.9988 + 0.0000i
% -0.6078 + 0.6740i    -0.6078 - 0.6740i
% -0.3241 + 0.0000i

% ----------------------
% (2) Graph the eigenvalues and transmission zeros in the complex plane.

% Computing the eigenvalues of A, for the 7-state system
EigVec = eig(A2);

% Plotting the eigenvalues and transmission zeros of the 7-state
% model in the complex plane.
% Eigenvalues are denoted with an 'x' and transmission zeros with an 'o'
% All eigenvalues lie inside the unit circle centered at the origin
figure;
f_PlotTZ(TransZ,EigVec,0);

% ----------------------
% (3)

Lambdac = zeros(1,7); % Initializing Lambdac vector

% Converting discrete-time eigenvalues to continuous-time eigenvalues.
for i = 1:1:7
    Lambdac(i) = log(EigVec(i))/ts;
end
% Lambdac

% The eigenvalues are:
% -3.6842 +71.1394i  -3.6842 -71.1394i   -3.5800 +72.2737i 
% -3.5800 -72.2737i  -10.4604 + 0.0000i  -7.6568 + 0.0000i
% -5.6364 + 0.0000i

% Three of the continuous-time eigenvalues lie on the real axis and the
% remaining four are two complex-conjugate pairs.
% There are 2 damped oscillators in the model.
% The natural frequency of a pole corresponds to its modulus and the
% damping ratio to minus the cosine of its phase.

NaturalFreq = abs(Lambdac);
% The natural frequencies of the oscillators are 71.2348 Hz and 72.3623 Hz. 

% ----------------------
% (4)

B2t1 = [B2(:,1),zeros(1,length(B2))'];
B2t2 = [zeros(1,length(B2))',B2(:,2)];
C2t1 = [C2(1,:);zeros(1,length(B2))];
C2t2 = [zeros(1,length(B2));C2(2,:)];

% Transmission Zeros for each channel
TransZ1 = f_TransZeros(A2,B2t1,C2t1,D2);
TransZ2 = f_TransZeros(A2,B2t1,C2t2,D2);
TransZ3 = f_TransZeros(A2,B2t2,C2t1,D2);
TransZ4 = f_TransZeros(A2,B2t2,C2t2,D2);

% Plotting the channel transmission zeros
figure;
f_PlotTZ(TransZ1,EigVec,1);
f_PlotTZ(TransZ2,EigVec,2);
f_PlotTZ(TransZ3,EigVec,3);
f_PlotTZ(TransZ4,EigVec,4);

% Hankel singular values of each channel
HSV1 = f_DiscHSV(A2,B2t1,C2t1);
HSV2 = f_DiscHSV(A2,B2t1,C2t2);
HSV3 = f_DiscHSV(A2,B2t2,C2t1);
HSV4 = f_DiscHSV(A2,B2t2,C2t2);

% HSV1 = [0.0415,0.0264,0.0156,0.0000,0.0000,0.0000,0.0000];
% HSV2 = [0.0188,0.0081,0.0001,0.0001,0.0000,0.0000,0.0000];
% HSV3 = [0.0472,0.0375,0.0199,0.0000,0.0000,0.0000,0.0000];
% HSV4 = [0.0382,0.0376,0.0128,0.0063,0.0000,0.0000,0.0000];

% ----------------------
% (5)

% The 8-state model
A8 = BigA(1:dim(5),1:dim(5),5);
B8 = BigB(1:dim(5),1:2,5);
C8 = BigC(1:2,1:dim(5),5);
D8 = D2;

TransZ8 = f_TransZeros(A8,B8,C8,D8); % Transfer zeros
EigVec8 = eig(A8); % Eigenvalues
figure;
f_PlotTZ(TransZ8,EigVec8,0); % Plotting 8-model

B8t1 = [B8(:,1),zeros(1,length(B8))'];
B8t2 = [zeros(1,length(B8))',B8(:,2)];
C8t1 = [C8(1,:);zeros(1,length(B8))];
C8t2 = [zeros(1,length(B8));C8(2,:)];

% Transmission Zeros for each channel
TransZ1_8 = f_TransZeros(A8,B8t1,C8t1,D8);
TransZ2_8 = f_TransZeros(A8,B8t1,C8t2,D8);
TransZ3_8 = f_TransZeros(A8,B8t2,C8t1,D8);
TransZ4_8 = f_TransZeros(A8,B8t2,C8t2,D8);

% Plotting the channel transmission zeros
figure;
f_PlotTZ(TransZ1_8,EigVec8,1);
f_PlotTZ(TransZ2_8,EigVec8,2);
f_PlotTZ(TransZ3_8,EigVec8,3);
f_PlotTZ(TransZ4_8,EigVec8,4);

%% ------------------------------------------------
% Task #3: Block diagram from analysis of individual channels

% See Report

%% ------------------------------------------------
% Task #4: Impulse response identification from white noise inputs

% Loading data
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts = 1/40; % Time step
N = length(y1); % Length of output vector
t = [0:N-1]*ts - 1; % Time vector

% ----------------------
% (1) Verifying that each input sequence is approximately zero mean

u1_m = mean(u1); % -9.6619e-04
u2_m = mean(u2); % 0.0013

% ----------------------
% (2) Computing Ruu for various lags

step = 0.1; % Lag step size
LagVec = -5:step:5; % Vector of lag values
Len = length(u1); % Length of input vectors
lenp = (length(u1)-1)/2; % p in the summation
RuuVec = f_R(u1,u2,u1,u2,LagVec); % Ruu

% Plotting Ruu as the lag varies
figure;
hold on; grid on;
axis equal;
title('Output Auto-correlation vs Lag','FontSize',18,'Interpreter','Latex');
xlabel('Lag (s)','FontSize',14,'Interpreter','Latex');
ylabel('Ruu (s)','FontSize',14,'Interpreter','Latex');
for i = 1:1:4
    plot(LagVec,RuuVec(i,:));
end

% ----------------------
% (3) Showing that Ruu0 ~= 4*eye(2)

Ruu0 = reshape(RuuVec(:,(length(RuuVec)-1)/2+1),[2,2]);

% Ruu0 =  [3.9855    0.0437;...
%          0.0437    4.0194];

% Ruu0 ~= [4    0;...
%          0    4];

% ----------------------
% (4) Plotting Ryu

step = 0.025; % Step size
LagVec = -0.2:step:2; % Vector of lags
RyuVec = f_R(u1,u2,y1,y2,LagVec); % Ryu

var1 = Ruu0(1,1); % Variance of the 1st channel of u 
var2 = Ruu0(2,2); % Variance of the 2nd channel of u 

% Normalizing the columns of Ryu by var1 and var2, respectively
RyuVec = [RyuVec(1,:)/var1; RyuVec(2,:)/var2;...
          RyuVec(3,:)/var1; RyuVec(4,:)/var2];

% Plotting Ryu for the four channels, as the lag varies
figure;
hold on; grid on;
title('Cross-correlation vs Lag','FontSize',18,'Interpreter','Latex');
xlabel('Lag (s)','FontSize',14,'Interpreter','Latex');
ylabel('Ryu (s)','FontSize',14,'Interpreter','Latex');
for i = 1:1:4
%     subplot(2,2,i);
    plot(LagVec,RyuVec(i,:));    
%     grid on;
end
legend('11','12','21','22');


%% ------------------------------------------------
% Task #5: H_2 norm analysis of identified model

% (1) The RMS value of the scaled output data y

Y = [y1;y2]; % Concatenating y1, y2
U = [u1;u2]; % Concatenating u1, u2

Yrms = rms(Y'/2); % Root-mean-square of [y1 y2]^T
% y1 and y2 are scaled by 2
% Yrms = [0.1812, 0.1312];
Yrms_n = norm(Yrms); % ||Yrms|| = 0.2237

% ----------------------
% (2) % Computing ||P||_H2 using 7-state model

% Computing (||P||_H2)^2 using the sum in (9)
UpperLim = Len; % Upper limit for the sum
SumPh2_9 = 0; % Initializing the sum at 0
for k = 1:1:UpperLim
    summand = trace((B2')*(A2')^k*(C2')*C2*A2^k*B2); % Summand in (9)
    SumPh2_9 = SumPh2_9 + summand; % Adding to the sum
end
PH_t_9 = sqrt(SumPh2_9); % ||P||_H2 = 0.0946

% Computing (||P||_H2)^2 using the sum in (10)
UpperLim = Len; % Upper limit for the sum
SumPh2_10 = 0; % Initializing the sum at 0
for k = 1:1:UpperLim
    summand = trace(C2*A2^k*B2*(B2')*(A2')^k*(C2')); % Summand in (10)
    SumPh2_10 = SumPh2_10 + summand; % Adding to the sum
end
PH_t_10 = sqrt(SumPh2_10); % ||P||_H2 = 0.0946

% ----------------------
% (3) Approximating ||P||_H2 from (8) using
% experimental pulse response data

% The square of the norm of the root-mean-square of the output vector
Yrms_n2_1 = (norm(Yrms))^2; % ||Yrms||^2 = 0.0500

% Check: ||Yrms||^2 = trace(Ryy[0])
LagVec = 0; % Zero lag
RyyVec = f_R(y1/2,y2/2,y1/2,y2/2,LagVec); % Ryy[0] in vector form
% y1 and y2 are scaled by 2
Yrms_n2_2 = trace(reshape(RyyVec,[2 2])); % trace(Ryy[0]) = 0.0500
% Both methods return the same value

% (||P||_H2)^2 = ||Yrms||^2
PH_e = sqrt(Yrms_n2_1); % ||P||_H2 = 0.2237

%% ------------------------------------------------
% Task #6: H_infinity norm

% (1) - (3)

tol = 10e-4; % Tolerance in the bisection in the f_Hinf_disc function
MZero = 10e-15; % Machine zero, also for use in f_Hinf_disc

% Computing the H infinity norm of the identified discrete-time model,
% and the frequency at which it is achieved.
[P_Hinf,w_max] = f_Hinf_disc(A2,B2,C2,D2,ts,tol,MZero);
% P_Hinf = 0.4705
% w_max = 11.3356

 % (4) Plots of frequency responses and singular values

N = length(y11f);
w = (0:N-1)/(ts*N); % Vector of angular frequencies

% Computing magnitudes, phases, and maximum singular values
% of the frequency response matrix.
[mag,phase,maxS] = f_FreqResponse(A2,B2,C2,D2,w,ts);

% Picking out each input-output channel
[m11,p11] = f_ReshapeFreq(mag,phase,w,1,1);
[m12,p12] = f_ReshapeFreq(mag,phase,w,1,2);
[m21,p21] = f_ReshapeFreq(mag,phase,w,2,1);
[m22,p22] = f_ReshapeFreq(mag,phase,w,2,2);

% Set vector to hold empirical freq. response values
FreqR_e = zeros(2,2,length(y11f));
FreqR_e(1,1,:) = y11f; FreqR_e(1,2,:) = y12f;
FreqR_e(2,1,:) = y21f; FreqR_e(2,2,:) = y22f;

% Vector to hold max singluar values
maxS_e = zeros(1,length(y11f));
for i = 1:1:length(y11f)
   maxS_e(i) = svds(FreqR_e(1:2,1:2,i),1);
end

% Plotting the magnitudes of each component of the frequency response
% matrix, the singular values of this matrix, and the singular values
% of the empirical frequency response matrix.
figure;
grid on; hold on;
plot(w,m11,'b');
plot(w,m12,'r--');
plot(w,m21,'m.');
plot(w,m22,'g.-');
plot(w,maxS,'k--');
scatter(w,maxS_e,'k');
line([w_max w_max], [0 1.1*max(maxS)],'Color','black','LineStyle','--');
xlim([0 20]);
ylim([0 1.1*max(maxS)]);
lgd = legend('Freq Resp. 11','Freq Resp. 12','Freq Resp. 21',...
             'Freq Resp. 22','Max Model Singular Values',...
             'Max Empirical Singular Values','w_{max}');
lgd.FontSize = 14;
% lgd.Interpreter = 'Latex';
xlabel('Angular Freq. (rad)','FontSize',18,'Interpreter','Latex');
ylabel('Gain (db)','FontSize',18,'Interpreter','Latex');
title('Frequency Responses and Singular Values','FontSize',20,...
      'Interpreter','Latex');
xlim([0 20]);

%% ------------------------------------------------
% Functions

% Calculates impulse responses
function hvec = f_ImpulseResponse(A,B,C,t)
    hvec = zeros(2,2,length(t)); % Initialize vector
    for k = 1:1:length(t)
        hk = C*(A^(k-1))*B;  % Impulse response
        % Reshaping into 3D vector
        hvec(1,1,k) = hk(1,1); hvec(1,2,k) = hk(2,1);
        hvec(2,1,k) = hk(1,2); hvec(2,2,k) = hk(2,2);
    end
end

% Computes the ijth element of the kth Impulse Response
function hk_ij = f_Impulseij(hkvec,t,d_t,i,j)
    hk_ij = reshape(hkvec(i,j,:),[1,length(t)]);
    hk_ij = [zeros(1,d_t),hk_ij];
    hk_ij(402:end) = [];
end

% Plots the impulse responses of systems k = [1, 4]
function f_PlotImpulse(y11,y12,y21,y22,t,d_t,hvec_,k)
    hjvec = hvec_(1:2,1:2,1:401,k); % 4D mat of impulse response mats
    % Plotting y11
    figure;
    subplot(221);
    hold on; grid on;
    plot(t,y11,'r*','LineWidth',2);
    plot(t,f_Impulseij(hjvec,t,d_t,1,1));
    xlim([-0.5 1.5]);
    xlabel('Time (s)','FontSize',14,'Interpreter','Latex');
    ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
    % Plotting y21
    subplot(222);
    hold on; grid on;
    plot(t,y21,'r*','LineWidth',2);
    plot(t,f_Impulseij(hjvec,t,d_t,1,2));
    xlim([-0.5 1.5]);
    xlabel('Time (s)','FontSize',14,'Interpreter','Latex');
    ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
    % Plotting y12
    subplot(223);
    hold on; grid on;
    plot(t,y12,'r*','LineWidth',2);
    plot(t,f_Impulseij(hjvec,t,d_t,2,1));
    xlim([-0.5 1.5]);
    xlabel('Time (s)','FontSize',14,'Interpreter','Latex');
    ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
    % Plotting y22
    subplot(224);
    hold on;grid on;
    plot(t,y22,'r*','LineWidth',2);
    plot(t,f_Impulseij(hjvec,t,d_t,2,2));
    xlim([-0.5 1.5]);
    xlabel('Time (s)','FontSize',14,'Interpreter','Latex');
    ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
    set(gca,'FontSize',14);
end

% Plots the ijth element of the impulse response of systems k = [1, 4]
% each impulse response element of all 4 models can be compared on the
% same plot.
function f_PlotImpulseij(t,d_t,hvec_,i,j)
    hold on;
    for k = 1:1:4
        hjvec = hvec_(1:2,1:2,1:401,k); % 4D mat of impulse response mats
        plot(t,f_Impulseij(hjvec,t,d_t,i,j));
    end
    ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
    grid on;
end

% Calculates the Magnitude and Phase of the Frequency Response
function [mag,phase,maxS] = f_FreqResponse(A,B,C,D,w,ts)
    mag = zeros(2,2,length(w)); % Initialize vector
    phase = zeros(2,2,length(w)); % Initialize vector
    maxS = zeros(1,length(w));
    for i = 1:length(w)
       % Frequency response
       FreqResp = C*((exp(1i*2*pi*w(i)*ts)*eye(length(A))-A)\B)+D;
       mag(1:2,1:2,i) = abs(FreqResp); % Magnitude
       phase(1:2,1:2,i) = angle(FreqResp); % Phase
       maxS(1,i) = svds(FreqResp,1);
    end
end

% Reshapes the elements of the magnitude and phase matrices into vectors
function [mij,pij] = f_ReshapeFreq(mag,phase,w,i,j)
    mij = reshape(mag(i,j,:),[length(w),1]);
    pij = reshape(phase(i,j,:),[length(w),1]);
end

% Plots the frequency responses
function f_PlotFreqResponse(m1,m2,m3,m4,me,w)
    hold on; grid on;
    plot(w,m1); plot(w,m2);
    plot(w,m3); plot(w,m4);
    plot(w,me);
    xlim([0 20]);
    xlabel('Angular Freq. (rad)','FontSize',14,'Interpreter','Latex');
    ylabel('Gain (db)','FontSize',14,'Interpreter','Latex');
    legend('6-model','7-model','10-model','20-model','empirical');
end

% Calculates Transmission Zeros
function TransZ = f_TransZeros(Am,Bm,Cm,Dm)
    Aam = [Am,Bm;-Cm,-Dm]; % LHS matrix in gen. eigvec/val eq.
    % RHS matrix in gen. eigvec/val eq.
    Bbm = [eye(size(Am)),zeros(size(Bm));zeros(size(Cm)),zeros(size(Dm))];
    [~,D] = eig(Aam,Bbm); % Solving general eigenvector/value equation
    DiagD = diag(D); % Rearranging diagonal matrix D into a vector DiagD

    % Removing infinite elements:
    TransZ = zeros(1,9); % Initializing TransZ vector
    j = 1;
    for i = 1:1:9
        % If eigenvalue is not infinite
        if DiagD(i) ~= Inf+0.0000i && DiagD(i) ~= -Inf+0.0000i
            TransZ(j) = DiagD(i); % Keep eigenvalue
            j = j+1; % Go to next element  of TransZ
        end
    end
    TransZ = rmmissing(TransZ); % Remove missing elements
    TransZ = nonzeros(TransZ); % Remove zero elements
    % Arrange in descending order, by order of modulus
    TransZ = sort(TransZ,'descend','ComparisonMethod','auto');
end

% Plots Transfer Zeros and Eigenvalues in the Complex Plane
function f_PlotTZ(TZ,EVec,k)
    if k ~= 0
        subplot(2,2,k);
    end
    hold on;
    grid on;
    axis equal;
    % Plotting transfer zeros
    for i = 1:1:length(TZ)
        x = real(TZ(i)); % x-coord: real part
        y = imag(TZ(i)); % y-coord: imag part
        plot(x,y,'bo','MarkerSize',8,'Linewidth',2);
    end
    % Plotting eigenvalues
    for i = 1:1:7
        x = real(EVec(i)); % x-coord: real part
        y = imag(EVec(i)); % y-coord: imag part
        plot(x,y,'rx','MarkerSize',8,'Linewidth',2);
    end
    
    % Plotting a unit circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
    plot(xunit, yunit,'k--');
    
    xlabel('x','FontSize',14,'Interpreter','Latex');
    ylabel('y','FontSize',14,'Interpreter','Latex');
    title('Pole-Zero Plot','FontSize',18,'Interpreter','Latex');
end

% Computes the Hankel Singular Values of a discrete-time system
function HSV = f_DiscHSV(As,Bs,Cs)
    % Transpose matrices
    Ast = transpose(As); Bst = transpose(Bs);Cst = transpose(Cs);

    % Definitions of discrete-time Grammians
    % As1*Gc*As1t - Gc + Bs1*Bs1t == 0
    % As1t*Go*As1 - Go + Cs1t*Cs1 == 0
    
    Gc = dlyap(As,Bs*Bst);  % Controllability Grammian
    Go = dlyap(Ast,Cst*Cs); % Observability Grammian
    
    HSV = svd(Gc*Go); % Hankel singular values
end

% Computes autocorrelations and cross correlation:
% Ruu, Ryy, Ryu/Ruy
% Inputs v1,v2/w1,w2 are placeholders for 1st,2nd elements of u/y vectors.
function RVec = f_R(v1,v2,w1,w2,LagV)
    kVec = round(LagV*40,2); % Converts double to natural number
    Len = length(v1); % Length of input vectors
    lenp = (length(v1)-1)/2; % p in the summation
    RSumVec = zeros(4,length(LagV)); % Initialize vec to hold sum values
    for j = 1:1:length(kVec) % Loop over the lag values
        Lag = kVec(j); % Lag value for this iteration
        % Initialize variables for summations
        RSum11 = 0; RSum12 = 0;
        RSum21 = 0; RSum22 = 0;
        for i = 1:1:Len % Loop over elements of input vectors
            % If either input vector has an element out of range,
            % don't add value to sum.
            if i+Lag <= 0 || i+Lag >= Len
                void();
            else % Else, add u product to sum
                RSum11 = RSum11 + v1(i)*w1(i+Lag);
                RSum12 = RSum12 + v1(i)*w2(i+Lag);
                RSum21 = RSum21 + v2(i)*w1(i+Lag);
                RSum22 = RSum22 + v2(i)*w2(i+Lag);
            end
            % Place sums in RSumVec
            RSumVec(1,j) = RSum11;
            RSumVec(2,j) = RSum12;
            RSumVec(3,j) = RSum21;
            RSumVec(4,j) = RSum22;
        end
    end
    RVec = RSumVec/(2*lenp); % Auto/cross correlation
end

% Accepts as inputs the state-space matrices of a continuous-time system,
% upper and lower limits for the gamma search, and a tolerance,
% and then returns the H_inf norm of the system computed to within the
% specified tolerance.
function [gamma,eig_clp] = f_Hinf_cont(A,B,C,D,tol,MZero)

    n = length(A); % Dimension of A matrix
    sys = ss(A,B,C,D);
    Gc = gram(sys,'c'); % Controllability Grammian
    Go = gram(sys,'o'); % Observability Grammian
    lb = svds(D,1); % Lower bound for binary search
    % lb = max(svds(D,1),sqrt(trace(Go*Gc/n))); % Lower bound 2
    ub = svds(D,1) + 2*sqrt(trace(n*Go*Gc)); % Upper bound

    j = 0;
    while (ub-lb)/2 > tol % While interval larger than tolerance
        gamma = (ub+lb)/2; % Set gammas as middle point in interval
        Dg = gamma^2*eye(length(D))-D'*D; % D_gamma
        A11 = A+B*(Dg\D')*C;
        A12 = -B*(Dg\B');
        A21 = C'*C+C'*D*(Dg\D')*C;
        A22 = -A'-C'*D*(Dg\B');
        A_clp = [A11,A12;A21,A22]; % Closed-loop A matrix
        e_clp = eig(A_clp); % Eigenvalues of A_clp
        EigVec(1:length(e_clp),j+1) = e_clp;
        
        % Counting the number of strictly imaginary
        % eigenvalues of A_clp.
        count = 0;
        for i = 1:1:length(e_clp)
          % if real(e_clp(i)) == 0 && imag(e_clp(i)) ~= 0
            if abs(real(e_clp(i))) <= MZero...
            && abs(imag(e_clp(i))) >= MZero
                count = count+1;
            end
        end
        % If any eigenvalue of A_clp is purely imaginary,
        % the H_inf norm of P is greater than or equal to gamma,
        % so set the new lower bound in the binary search
        % to this value of gamma.
        if count > 0
            side(j+1) = 1; % Right side of gamma
            lb = gamma;
        % Otherwise, the H_inf norm of P is less than gamma,
        % so set the new upper bound in the binary search
        % to this value of gamma.
        else
            side(j+1) = 0; % Left side of gamma
            ub = gamma;
        end
        j = j+1;
        GVec(j) = gamma; % Vector to hold values of gamma
    end
    
    cc = find(side,1,'last');
    eig_clp = EigVec(1:length(e_clp),cc);
    
    % Plotting gamma at each of the iterations
    figure;
    plot(1:1:length(GVec),GVec);
    hold on;
    plot(1:1:length(GVec),GVec,'bo','MarkerFaceColor','b');
    grid on;
    xlabel('Iteration','FontSize',14,'Interpreter','Latex');
    ylabel('Gamma','FontSize',14,'Interpreter','Latex');
    title('Binary Search for Gamma','FontSize',18,'Interpreter','Latex');
end

% Accepts as inputs the state-space matrices of a discrete-time system,
% upper and lower limits for the gamma search, a tolerance, and the sample
% period ts associated with the system, and then returns the H_inf norm
% of the system computed to within the specified tolerance,
% and the approximate frequency at which the maximum gain is achieved.
function [gamma,w_max] = f_Hinf_disc(A,B,C,D,ts,tol,MZero)
    % Converting discrete system to equivalent continuous system
    [Ac,Bc,Cc,Dc,IA] = f_Disc2Cont(A,B,C,D); 

    % Calling the continuous-time function to compute gamma
    [gamma,eig_clp] = f_Hinf_cont(Ac,Bc,Cc,Dc,tol,MZero);
    
    % Computing the approximate frequency at which
    % the maximum gain is achieved.
    w = zeros(1,length(eig_clp));
    for i = 1:1:length(eig_clp)
        if abs(real(eig_clp(i))) <= 10e-10...
        && abs(imag(eig_clp(i))) >= 10e-10
            v = imag(eig_clp(i));
            % disp(v);
            w(i) = -1i*log((1+1i*v)/(1-1i*v))/ts;
        else
            void();
        end
    end
    % Frequency at which the the largest singular value
    % of the frequency response occurs
    w0 = max(real(w))/(2*pi);
    %b disp(w0)
    
    % The Nyquist frequency w_nyq was not covered by the
    % preceding procedure, so we must compare the frequency
    % responses at w0 and w_nyq.
    w_nyq = 1/(2*ts); % Nyquist Frequency
    P_0 = C*((1i*w0*IA-A)\B)+D; % Frequency response, at w0
    P_nyq = C*((1i*w_nyq*IA-A)\B)+D; % Frequency response, at w_nyq
    P_Hinf_0 = svds(P_0*exp(1i*w0*ts),1); % H_inf norm of P, for w0
    P_Hinf_nyq = svds(P_nyq*exp(1i*w_nyq*ts),1);  % H_inf norm of P, for w_nyq
    
    % disp(P_Hinf_0)
    % disp(P_Hinf_nyq)
    
    % Taking the frequency that maximizes the frequency response.
    if P_Hinf_0 > P_Hinf_nyq
        w_max = w0;
    else
        w_max = w_nyq;
    end
end

% Converts discrete system to equivalent continuous system
function [Ac,Bc,Cc,Dc,IA] = f_Disc2Cont(A,B,C,D)
    IA = eye(length(A)); % Identity matrix
    % Continuous-time state-space matrices
    Ac = -(IA+A)\(IA-A);
    Bc = sqrt(2)*((IA+A)\B);
    Cc = sqrt(2)*(C/(IA+A));
    Dc = D-C*((IA+A)\B);
end

% Void Function
function [] = void()
end

% ------------------------------------------------
