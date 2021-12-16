% Test Homework 11

load('OFB_Example_2');

Plant = Plant_Ex2;
Disturbance = ss(1,1,1,0);
ny1 = 1;

[Contr,sysCL] = assign11_PeterRacioppo(Plant,Disturbance,Q,R,W,ny1);

e1 = max(abs(Contr_Ex2_sg.a - Contr.a))
e2 = max(abs(Contr_Ex2_sg.b - Contr.b))
e3 = max(abs(Contr_Ex2_sg.c - Contr.c))
e4 = max(abs(Contr_Ex2_sg.d - Contr.d))

e5 = max(abs(sysCL_Ex2_sg.a - sysCL.a))
e6 = max(abs(sysCL_Ex2_sg.b - sysCL.b))
e7 = max(abs(sysCL_Ex2_sg.c - sysCL.c))
e8 = max(abs(sysCL_Ex2_sg.d - sysCL.d))

