clear all;
close all;
clc

[Y,U,IsYout] = Simulate_data_Pessiglione2006();

[posterior_M1,out_M1] = invert_data_Pessiglione2006_M1(Y,U,IsYout);
[posterior_M2,out_M2] = invert_data_Pessiglione2006_M2(Y,U,IsYout);
[posterior_M3,out_M3] = invert_data_Pessiglione2006_M3(Y,U,IsYout);

F = [out_M1.F,out_M2.F,out_M3.F]

bar(F)
title('Approximated log_evidence')