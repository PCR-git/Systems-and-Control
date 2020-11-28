
% Linear Systems
% HW 6

% Question 1

% syms t
% simplify(expm([0,t;-t,0]))

alpha = 1;
T =1 ;
X0 = eye(2,2);
XT = [cos(T),sin(T)*exp(-alpha);-sin(T),cos(T)*exp(-alpha)];
Eigs = eig(XT);
abs(Eigs);

% iter = 15;
iter = 5;
Pos1 = zeros(2,iter);
Pos2 = zeros(2,iter);
for i = 1:iter

Pos1i = XT^i*X0(:,1);
Pos1(:,i) = Pos1i;

Pos2i = XT^i*X0(:,2);
Pos2(:,i) = Pos2i;

end

figure;
plot(Pos1(1,:),Pos1(2,:),'b');
hold on;
plot(Pos2(1,:),Pos2(2,:),'r--');

axis equal;
grid on;
xlabel('x','FontSize',20,'FontAngle', 'italic');
ylabel('y','FontSize',20,'FontAngle', 'italic');
title('Phase Plane Portrait (15 Periods)');
set(gca, 'FontSize',20,'FontName', 'Times New Roman');

figure;
subplot(2,2,1);
plot(1:iter,Pos1(1,:),'b');
axis equal;
grid on;
xlabel('Periods','FontSize',20,'FontAngle', 'italic');
ylabel('x1','FontSize',20,'FontAngle', 'italic');
title('X Pos, Init. Cond. 1');
set(gca, 'FontSize',20,'FontName', 'Times New Roman');

subplot(2,2,2);
plot(1:iter,Pos1(2,:),'b');
axis equal;
grid on;
xlabel('Periods','FontSize',20,'FontAngle', 'italic');
ylabel('y1','FontSize',20,'FontAngle', 'italic');
title('Y Pos, Init. Cond. 1');
set(gca, 'FontSize',20,'FontName', 'Times New Roman');

subplot(2,2,3);
plot(1:iter,Pos2(1,:),'b');
axis equal;
grid on;
xlabel('Periods','FontSize',20,'FontAngle', 'italic');
ylabel('x2','FontSize',20,'FontAngle', 'italic');
title('X Pos, Init. Cond. 2');
set(gca, 'FontSize',20,'FontName', 'Times New Roman');

subplot(2,2,4);
plot(1:iter,Pos2(2,:),'b');
axis equal;
grid on;
xlabel('Periods','FontSize',20,'FontAngle', 'italic');
ylabel('y2','FontSize',20,'FontAngle', 'italic');
title('Y Pos, Init. Cond. 2');
set(gca, 'FontSize',20,'FontName', 'Times New Roman');

% ------------------------------------

% Question 2
% (a)

syms P11 P12 P21 P22 a b alpha x1

A = [0,1;-1,-0.5];

Q1 = [1,0;0,1];
Q2 = [2,1;1,1];
Q3 = [1,1;1,10];

X = [Q1,Q2,Q3];

P = [P11,P12;P21,P22];

LHS = P*A + transpose(A)*P;

for i = 1:3

Q = X(:,2*i-1:2*i);

eqn1 = LHS(1,1) == -Q(1,1);
eqn2 = LHS(1,2) == -Q(1,2);
eqn3 = LHS(2,1) == -Q(2,1);
eqn4 = LHS(2,2) == -Q(2,2);

Soln = solve([eqn1,eqn2,eqn3,eqn4],[P11 P12 P21 P22]);

P11_s = Soln.P11;
P12_s = Soln.P12;
P21_s = Soln.P21;
P22_s = Soln.P22;

P_s = [P11_s,P12_s;P21_s,P22_s];

end

P1 = [9/4, 1/2;...
      1/2,   2];

P2 = [5/2, 1;...
      1,   3];
 
P3 = [41/4, 1/2;...
      1/2,   11];
  
Pee = [P1,P2,P3];

% ------------------
% (b)

Eqn1 = (P12+P21)^2*x1^2 - 4*P22*(P11*x1^2-10) == 0;
solve(Eqn1);

Bounds = [(2*10^(1/2)*(-P22*(P12^2 + 2*P12*P21 + P21^2 - 4*P11*P22))^(1/2))/(P12^2 + 2*P12*P21 + P21^2 - 4*P11*P22);...
         -(2*10^(1/2)*(-P22*(P12^2 + 2*P12*P21 + P21^2 - 4*P11*P22))^(1/2))/(P12^2 + 2*P12*P21 + P21^2 - 4*P11*P22)];


for i = 1:3
    
P = Pee(:,2*i-1:2*i);

P11_subs = P(1,1);
P12_subs = P(1,2);
P21_subs = P(2,1);
P22_subs = P(2,2);

subs(Bounds,[P11,P12,P21,P22],[P11_subs,P12_subs,P21_subs,P22_subs]);

end

Bound1 = [-(2*10^(1/2)*34^(1/2))/17,(2*10^(1/2)*34^(1/2))/17];
Bound2 = [-(10^(1/2)*78^(1/2))/13,(10^(1/2)*78^(1/2))/13];
Bound3 = [-(10^(1/2)*4950^(1/2))/225,(10^(1/2)*4950^(1/2))/225];

Bound_Vec = [Bound1;Bound2;Bound3];

names = ['P1','P2','P3'];
for i = 1:3
    
P = Pee(:,2*i-1:2*i);

P11_subs = P(1,1);
P12_subs = P(1,2);
P21_subs = P(2,1);
P22_subs = P(2,2);

Boundi = Bound_Vec(i,:);
Boundi(1) = Boundi(1);
Boundi(2) = Boundi(2);

x1 = linspace(Boundi(1),Boundi(2),200);
x2_1 = subs((1/(2*P22))*(-(P12+P21)*x1 + sqrt((P12+P21)^2*x1.^2 - 4*P22*(P11*x1.^2-10))),[P11,P12,P21,P22],[P11_subs,P12_subs,P21_subs,P22_subs]);
x2_2 = subs((1/(2*P22))*(-(P12+P21)*x1 - sqrt((P12+P21)^2*x1.^2 - 4*P22*(P11*x1.^2-10))),[P11,P12,P21,P22],[P11_subs,P12_subs,P21_subs,P22_subs]);

figure;
plot(x1,x2_1,'r')
hold on;
plot(x1,x2_2,'r')

  for j = 1:5
  
  pos1 = [x1(40*j-20);x2_1(40*j-20)];
  pos2 = [x1(40*j-20);x2_2(40*j-20)];
  scale = 1/4;
  V1 = scale*A*pos1;
  V2 = scale*A*pos2;
  
  scatter(pos1(1),pos1(2),25,'b','filled');
  scatter(pos2(1),pos2(2),25,'b','filled');
  quiver(pos1(1),pos1(2),V1(1),V1(2),'b');
  quiver(pos2(1),pos2(2),V2(1),V2(2),'b');
  
  end

axis equal;
grid on;
xlabel('x1','FontSize',20,'FontAngle', 'italic');
ylabel('x2','FontSize',20,'FontAngle', 'italic');
title(names(2*i-1:2*i));
set(gca, 'FontSize',20,'FontName', 'Times New Roman');

end

