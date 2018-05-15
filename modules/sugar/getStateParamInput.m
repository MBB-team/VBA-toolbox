function [state, parameter,input] = getStateParamInput(X,P,U,in)

% Get the hidden states, the free parameters and the inputs (using labels stored in options.in[FG]

if  isfield(in,'thetaLabel')
    parameter = priorPrettifyer(in.thetaLabel,P) ;
elseif isfield(in,'phiLabel')
    parameter = priorPrettifyer(in.phiLabel,P) ;
else
    parameter = P ;
end

if  isfield(in,'stateLabel')
    state = priorPrettifyer(in.stateLabel,X) ;
else
    state = X ;
end

if isfield(in, 'inputLabel')
    input = priorPrettifyer(in.inputLabel,U) ;
else
    input = U ;
end
    

