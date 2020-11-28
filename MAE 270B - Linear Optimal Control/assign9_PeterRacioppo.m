
% MAE 270B: Assignment 9
% Peter Racioppo

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
