function [fx, dF_dX] = f_SHC (Xt, ~, ~, inF)
% stable heteroclinic channels evolution function

deltat = inF.deltat;
x = Xt;

if ~ isfield (inF, 'K')
    inF.K = [1/2; 1/32];
end
if ~ isfield (inF, 'lambda')
    inF.lambda = 0.3;
end

if ~ isfield (inF, 'G0')
    inF.G0 = 50;
end
if ~ isfield (inF, 'beta')
    inF.beta = 0.5;
end

if ~ isfield (inF, 'ind1') || ~ isfield (inF, 'R1')
    inF.ind1 = 1 : 4;
    inF.R1{1} = ...
        [1   5   5   0.5
         0.5 1   5   5
         5   0.5 1   5
         5   5   0.5 1];
    inF.R1{2} = ...
        [1   5   5   0.5
         5   1   0.5 5
         0.5 5   1   5
         5   0.5 5   1];
    inF.R1{3} = ...
        [1   0.5 5   5
         5   1   0.5 5
         5   5   1   0.5
         0.5 5   1   5];
end
if ~ isfield (inF, 'ind2') || ~ isfield (inF, 'R2')
    inF.ind2 = 5 : 7;
    inF.R2 = ...
        [1   5   0.5
         0.5 1   5
         5   0.5 1 ];
end


% Separate states
x1 = x(inF.ind1);
x2 = x(inF.ind2);

% High level
[SX2, dsdx2] = VBA_sigmoid (x2, 'scale', inF.G0, 'slope', inF.beta);
ff{2} = inF.K(2) * (- inF.lambda * x2 - inF.R2 * SX2);

% Low level
[SX1, dsdx1] = VBA_sigmoid (x1, 'scale', inF.G0, 'slope', inF.beta);
R1eff = inF.R1{1} * SX2(1) + inF.R1{2} * SX2(2) + inF.R1{3} * SX2(3);
ff{1} = inF.K(1) * (- inF.lambda * x1 - R1eff * SX1);

f = [ff{1}; ff{2}];

df1 = zeros (size (f, 1), length (inF.ind1));
df2 = zeros (size (f, 1), length (inF.ind2));

df2(inF.ind2, :) = inF.K(2) * inF.lambda * (- eye (length (inF.ind2))) - inF.K(2) * inF.R2 * diag (dsdx2);
df1(inF.ind1, :) = inF.K(1) * inF.lambda * (- eye (length (inF.ind1))) - inF.K(1) * R1eff * diag (dsdx1);

df_dx = [df1, df2]';
df_dx(inF.ind2(1), inF.ind1) = (- inF.K(1) * dsdx2(1) * inF.R1{1} * SX1)';
df_dx(inF.ind2(2), inF.ind1) = (- inF.K(1) * dsdx2(2) * inF.R1{2} * SX1)';
df_dx(inF.ind2(3), inF.ind1) = (- inF.K(1) * dsdx2(3) * inF.R1{3} * SX1)';

fx = Xt + deltat * f;
dF_dX = eye (length (Xt)) + deltat * df_dx;

