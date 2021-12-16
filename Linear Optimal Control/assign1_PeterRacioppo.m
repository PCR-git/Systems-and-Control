
function W = assign1_PeterRacioppo(A,B,C,D,K)

% (Both options work)
% Option 1
W = lyap((A-B*K)', (C-D*K)'*(C-D*K));

% Option 2
% sys = ss(A-B*K,B,C,D);
% W = gram(sys,'o');

end