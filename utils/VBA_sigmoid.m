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
%       * 'slope'     : inverse-temperature parameter (default = 1)
%       * 'center'    : absissa of the inflexion point (default = 0)
%       * 'scale'     : multiplicative gain of the transformation (default = 1)
%       * 'offset'    : additive gain (default = 0)
%       * 'lapseRate' : shotcut for offset = lapseRate and scale = 1 - 2 * lapseRate
%                       this option is incompatible with offset or scale
%   - other optional key/value pairs or structure (or both)
%       * 'inverse' : reversed transformation, ie. x = VBA_sigmoid(y, opt) 
%                       (default = false)
%       * 'finite'  : boundaries that enforce precision to be finite and 
%                     derivative to stay numerically non zero (default = 1e-9)
%                     Set to 0 to deactivate
%       * 'derivatives' : cell array of parameter names wrt which dsdp must be
%                         computed (in given order). 
%                   
% 
%   Note that if the 'inverse' flag is set to true, the derivatives will
%   not be computed and dsdx = dsdp = [] will be returned instead.
%
% OUT:
%   - y: transformed values, with the same dimension as x.
%   - dsdx: derivative of the transformation wrt. x, taken at each point x.
%     Has the same dimension as x.
%   - dsdp: derivative of the transformation wrt. the parameters specified
%     as arguments. Derivatives are for parameters taken in alphabetical
%     order and aggregated along the first dimension. For example, calling
%     VBA_sigmoid([x1; x2], 'slope', 2, 'center', 3) will return dsdp as:
%     [ d_s/d_center(x1) d_s/d_center(x2) ;
%       d_s/d_slope(x1)  d_s/d_slope(x2)  ]
%     If x is multidimensional, size(dsdp) = [nb_params, size(x,1), size(x,2), ...]
%
% /////////////////////////////////////////////////////////////////////////


%% Globals
% =========================================================================
% truncature for finite sigmoid
epsilon = 1e-9;

%% Shrtcut
% =========================================================================
% quick version!
if nargin == 1 && nargout == 1
    y = epsilon + (1 - 2 * epsilon) ./ (1 + exp (- x));
    return
end

%% Parse arguments
% =========================================================================

% define inputParser
% -------------------------------------------------------------------------
persistent parser;

if isempty (parser)
    parser = getParser ();
end

% parse arguments
% -------------------------------------------------------------------------

parser.parse (varargin{:});
params = parser.Results ;

% apply shortcuts if needed
% -------------------------------------------------------------------------

if ~ ismember ('lapseRate', parser.UsingDefaults)
    % lapseRate is mutually ex
    if any (~ ismember ({'offset', 'scale'}, parser.UsingDefaults))
        error('*** VBA_sigmoid: you can not specify lapseRate and offset or scale at the same time');
    end
    params.offset = params.lapseRate;
    params.scale = 1 - 2 * params.lapseRate;
end

%% Compute transformation
% =========================================================================

%% Inversed case
if params.inverse
    % check that values are valid
    if ~ VBA_isInRange(x, [0 1])
        error('*** VBA_sigmoid: inverse sigmoid inputs must be between 0 and 1');
    end
    
    % inverse sigmoid transformation
    lx = params.scale * (x - params.offset) .^-1  - 1;
    y = params.center - params.slope ^-1 * log (lx) ;
    
    % skip derivatives, they are generally not used in inverse case
    dsdx = [];
    dsdp = [];
    
%% Normal case
else
    % evaluate sigmoid
    % ---------------------------------------------------------------------
    sx = 1 ./ (1 + exp(- params.slope * (x - params.center)));
    y = params.offset + params.scale * sx;
    
    % ensure finite precision ('finite' flag)
    % ---------------------------------------------------------------------
    if params.finite > 0
        minY = params.offset + epsilon;
        y = max(y,minY);
        maxY = params.offset + params.scale - epsilon;
        y = min(y,maxY);
    end
    
    % compute derivatives with respect to value
    % ---------------------------------------------------------------------

    % skip if not not required
    if nargout < 2
        return
    end
    
    % actual computation
    dsdx = params.slope * (y - params.offset) .* (1 - (y - params.offset) ./ params.scale);
    
    % compute derivatives with respect to parameters
    % ---------------------------------------------------------------------

    % skip if not not required
    if nargout < 3
        return
    end

    % concatenate derivatives for all parameters
    dims = size(x);
    
    dsdp = cat (2, ...
        - VBA_vec (dsdx), ... % d_center
        1 - 2 * VBA_vec (sx), ... d_lapseRate
        ones(numel(sx),1), ... d_offset
        VBA_vec (sx), ... d_scale
        ((VBA_vec (x) - VBA_vec(params.center)) / params.slope) .* VBA_vec (dsdx) ... d_slope
        );

    % keep only those passed s parameter
    derivables = {'center','lapseRate','offset','scale','slope'};
    if isempty(params.derivatives)
         params.derivatives = setdiff(derivables, parser.UsingDefaults);
    end
    dIdx = cellfun(@(l) find(strcmp (derivables, l)), params.derivatives);
    
    dsdp = dsdp(:,dIdx);
    
    % set derived parameter as first dimension
    if all(dims == 1)
        dsdp = dsdp';
    else
        dims(dims==1) = [];
        dsdp = reshape(dsdp', [size(dsdp,2) dims]);
    end


end

end

function parser = getParser()

parser = inputParser;
parser.PartialMatching = false;
parser.KeepUnmatched = true;

% define parameters
% -------------------------------------------------------------------------

% flags
parser.addParameter ('inverse', false, @islogical);
parser.addParameter ('finite', true, @(z) VBA_isInRange(z, [0 1e-2]));
parser.addParameter ('derivatives', {}, @iscellstr);

% x transfomations
parser.addParameter ('slope', 1, @isnumeric);
parser.addParameter ('center', 0, @isnumeric);

% sig transformation
parser.addParameter ('offset', 0, @isnumeric);
parser.addParameter ('scale', 1, @(z) VBA_isInRange(z, [eps Inf]));

% shortcut for lapse rate model
parser.addParameter ('lapseRate', 0, @(z) VBA_isInRange(z, [0 0.5]));
end
