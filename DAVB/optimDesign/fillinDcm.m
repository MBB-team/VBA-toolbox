function DCM = fillinDcm(DCM)
% fill in missing entries in DCM structure for review
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

DCM.M.m = size(DCM.c,2); % number of inputs
DCM.M.n = 5*size(DCM.Y.y,2); % number of hidden states
DCM.M.l = size(DCM.Y.y,2); % number of regions

DCM.y = zeros(size(DCM.Y.y)); % predicted data
DCM.R = DCM.Y.y; % residuals of the model
DCM.F = 'this model has not been estimated yet.';

DCM.Cp = zeros(DCM.M.l^2*(1+DCM.M.m+DCM.M.l)+DCM.M.l*(DCM.M.m+5));

DCM.Ep.A = DCM.a;
DCM.Ep.B = DCM.b;
DCM.Ep.C = DCM.c;
DCM.Ep.D = DCM.d;

DCM.Pp.A = 0.5*ones(size(DCM.a));
DCM.Pp.B = 0.5*ones(size(DCM.b));
DCM.Pp.C = 0.5*ones(size(DCM.c));
DCM.Pp.D = 0.5*ones(size(DCM.d));