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
