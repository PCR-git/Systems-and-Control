
% MAE 270B
% Peter Racioppo
% Assignment 2

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
    K = R\B'*P;

end

% ------------------------------


