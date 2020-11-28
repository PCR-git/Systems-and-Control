% test_assign10

clear all;
close all;
clc;

load('OFB_Example_1');
ny1 = 1;
[Contr, sysCL] = assign10_PeterRacioppo(Plant, Q, R, W, ny1);

max(max(sysCL_sg.a-sysCL.a))
max(max(sysCL_sg.b-sysCL.b))
max(max(sysCL_sg.c-sysCL.c))
max(max(sysCL_sg.d-sysCL.d))

max(max(sysCONTR_sg.a-Contr.a))
max(max(sysCONTR_sg.b-Contr.b))
max(max(sysCONTR_sg.c-Contr.c))
max(max(sysCONTR_sg.d-Contr.d))

