function [y, dsdx, dsdp] = VBA_sigmoid(x, varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [y, dsdx, dsdp] = VBA_sigmoid(x, [name, value, ...])
% Apply a sigmoid transformation to x
%
% By default, the canonical sigmoid is used:
%
%   y = 1 / (1 + exp(- x))
%
% However, this function can be parametrized using name/value pairs of
% arguments to get, using all options:
%
%   y = offset + scale * 1 / (1 + exp(- slope * (x - center))
%
% IN:
%   - x: values to be transformed
%   - optional key/value pairs or structure (or both) that parametrize the 
%     sigmoid transformation:
%       - 'slope'     : inverse-temperature parameter
%       - 'center'    : absissa of the inflexion point
%       - 'scale'     : multiplicative gain of the transformation
%       - 'offset'    : additive gain 
%       - 'lapseRate' : shotcut for offset = lapseRate and scale = 1 - 2 * lapseRate
%                       this option is incompatible with offset or scale
%       - 'inverse'   : reversed transformation, ie. x = VBA_sigmoid(y, opt)
% 
%   Note that if the 'inverse' flag is set to true, the derivatives will
%   not be computed and dsdx = dsdp = [] will be returned instead.
%
% OUT:
%   - y: transformed values, with the same dimension as x.
%   - dsdx: derivative of the transformation wrt. x, taken at each point x.
%     Has the same dimension as x.
%   - dsdp: derivative of the transformation wrt. the parameters specified
%     as arguments. Derivatives are taken for parameters in alphabetical
%     order and aggregated along the first dimension. For example, calling
%     VBA_sigmoid([x1; x2], 'slope', 2, 'center', 3) will return dsdp as:
%     [ d_s/d_center(x1) d_s/d_center(x2) ;
%       d_s/d_slope(x1)  d_s/d_slope(x2)  ]
%     If x is multidimensional, size(dsdp) = [nb_params, size(x,1), size(x,2), ...]
%
% /////////////////////////////////////////////////////////////////////////


%% 
% epsilon = 1e-8;

%% parse arguments
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
    lx = params.scale .* (x - params.offset) .^-1  - 1;
    y = params.center - params.slope .^-1 .* log (lx)  ;
    
    % skip derivatives, they are generally not used
    dsdx = [];
    dsdp = [];
else
    % evaluate sigmoid
    lx = params.slope * (x - params.center);
    sx = 1 ./ (1 + exp(- lx));
    y = params.offset + params.scale * sx;
    
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

    dims = size(x);
    
    dsdp = cat (2, ...
        - vec (dsdx), ... % d_center
        1 - 2 * vec (sx), ... d_lapseRate
        ones(numel(x),1), ... d_offset
        vec (sx), ... d_scale
        ((vec (x) - params.center) / params.slope) .* vec (dsdx) ... d_slope
        );

    idxDefault = ismember ({'center','lapseRate','offset','scale','slope'}, parser.UsingDefaults);
    
    dsdp(:,idxDefault) = [];
    
    dsdp = reshape(dsdp', [size(dsdp,2) dims]);


end
