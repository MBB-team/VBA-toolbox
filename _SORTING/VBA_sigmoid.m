function [y, dsdx, dsdp] = VBA_sigmoid(x, varargin)

%% 
epsilon = 1e-8;

%% specify arguments
parser = inputParser;
parser.PartialMatching = false;

parser.addParameter ('inverse', false, @islogical);

parser.addParameter ('slope', 1, @isnumeric);
parser.addParameter ('center', 0, @isnumeric);

parser.addParameter ('offset', 0, @isnumeric);
parser.addParameter ('scale', 1, @(z) VBA_isInRange(z, [eps Inf]));

parser.addParameter ('lapseRate', 0, @(z) VBA_isInRange(z, [0 0.5]));

parser.parse (varargin{:});
params = parser.Results ;

if ~ ismember ('lapseRate', parser.UsingDefaults)
    if any (~ ismember ({'offset', 'scale'}, parser.UsingDefaults))
        error('*** VBA_sigmoid: you can not specify lapseRate and offset or scale at the same time');
    end
    params.offset = params.lapseRate;
    params.scale = 1 - 2 * params.lapseRate;
end
%% compute
    
if params.inverse
    % evaluate inverse sigmoid  
    x = params.scale .* (x - params.offset) .^-1  - 1;
    y = params.center - params.slope .^-1 .* log (x)  ;
    
    % skip derivatives, they are generally not used
    dsdx = [];
    dsdp = [];
else
    % evaluate sigmoid
    x = params.slope * (x - params.center);
    y = params.offset + params.scale ./ (1 + exp(- x));
    
    % ensure finite precision
%     minY = params.offset + epsilon;
%     y(y <= minY) = minY;
%     
%     maxY = params.offset + params.scale - epsilon;
%     y(y >= maxY) = maxY;
%     
    % compute derivatives
    % evaluate derivative wrt x

    if nargout < 2
        return
    end
    
    dsdx = params.slope * (y - params.offset) .* (1 - (y - params.offset) / params.scale);
    
    % evaluate derivatives wrt parameters
    
    if nargout < 3
        return
    end

    dsdp = [];
    
% pdiff = setdiff(parser.Parameters, parser.UsingDefaults);
% 
% % evaluate derivatives wrt parameters
% if size(Phi,1) >= 1
%     dsdp = zeros(size(Phi,1),length(x));
%     dsdp(1,:) = (x-th).*in.beta.*dsdx;
%     if size(Phi,1) == 2
%         dsdp(2,:) = -dsdx;
%     end
%     dsdp(isnan(dsdp)) = 0;
% else
%     dsdp = [];
% end


end
