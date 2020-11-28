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
