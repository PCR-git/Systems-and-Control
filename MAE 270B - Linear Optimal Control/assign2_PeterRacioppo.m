
% MAE 270B
% Peter Racioppo

% Peter Racioppo
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