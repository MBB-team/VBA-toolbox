function [Sx,dsdx,dsdp] = sigm(x,in,Phi)

% legacy code
s = warning ('on');
warning ('*** The function `sigm` is now deprecated. Please see `VBA_sigmoid` for an alternative.') 
warning (s);

% fallback
if nargin < 2
    in = struct ();
end

if nargin < 3
    Phi = [];
end

if isfield(in, 'G0')
    in.scale = in.G0;
    in = rmfield(in,'G0');
end
 
if isfield(in, 'S0')
    in.offest = in.S0;
    in = rmfield(in,'S0');
end
    
if numel(Phi) > 0
    in.slope = exp(Phi(1));
    in.derivatives = {'slope'};
end

if numel(Phi) > 1
    in.center = Phi(2);
    in.derivatives = {'slope', 'center'};
end

[Sx,dsdx,dsdp] = VBA_sigmoid (x, in);

if numel(Phi) > 0
    dsdp(1,:) = dsdp(1,:) * exp(Phi(1));
end
