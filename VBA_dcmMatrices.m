function [A,B,C,D,dim] = VBA_dcmMatrices(options,theta)
% [A,B,C,D,dim] = VBA_dcmMatrices(out,theta)
% DCM parameters are often given in a vector form (like as the output of a
% model inversion usin VBA_NLStateSpaceModel) that is harder to read than
% the original matrix form.
% This function helps reconstructing all the connectivity matrices A, B, C,
% and D describing a DCM given a DCM structure (indicator matrices) and a
% set of connectivity weights given in vector form.
%
% IN:
%   - options: structure containing a description of the DCM, eg.  
%              the output of the prepare_fullDCM(...) function
%   - theta  : parameters used to fill in the DCM matrices (eg posterior.muTheta or prior.muTheta)
% OUT:
%   - A,B,C,D: matrices describing the DCM with actual connectivity weights
%   - dim    : dimensions of the DCM


try 
    inF = options.inF{1};
catch
    inF = options.inF;
end

[n,nu] = size(inF.C) ;

try
nr = numel(inF.r);
catch
    nr=0;
end
dim.n=n;
dim.nr=nr;
dim.nu=nu;

% inF=out.options.inF;
% theta = posterior.muTheta;

A = inF.A;
indA = inF.indA;
if ~isempty(indA)
    A(A~=0) = theta(indA);
    A = A - exp(theta(inF.indself)).*eye(n);
else
    A = A - exp(theta(inF.indself)).*eye(n);
end

B = inF.B;
indB = inF.indB;
for i=1:nu
    if ~isempty(indB{i})
        B{i}(B{i}~=0) = theta(indB{i});
    end
end

C = inF.C;
indC = inF.indC;
if ~isempty(indC)
    C(C~=0) = theta(indC);
end

D = inF.D;
indD = inF.indD;
for i=1:n
    if ~isempty(indD{i})
        D{i}(D{i}~=0) = theta(indD{i});
    end
end

return

%% ########################################################################
try
hself = - diag(exp(theta(inF.indhself))) ;


if out.options.extended
hA = inF.hA;
indhA = inF.indhA;
if ~isempty(indhA)
    hA(hA~=0) = theta(indhA);
end

hB = inF.hB;
indhB = inF.indhB;
for i=1:nu
    if ~isempty(indhB{i})
        hB{i}(hB{i}~=0) = theta(indhB{i});
    end
end

hC = inF.hC;
indhC = inF.indhC;
if ~isempty(indhC)
    hC(hC~=0) = theta(indhC);
end

A = [A  zeros(n,nr);
     hA hself*eye(nr)] ;

C = [C ; hC] ;
end

for i=1:nu
    B{i} = [[B{i}; hB{i}] zeros(n+nr,nr)];
end
end

end
