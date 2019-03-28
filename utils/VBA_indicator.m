function v = VBA_indicator (x, S, inverse)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% x = VBA_indicator (x, S, inverse)
% indicator (aka incidence or characteristic) vector of x wrt set S.
%
% IN:
%   - x: K x N values to be converted. Indicator vectors are computed for 
%       each columns of x
%   - S: full set of elements in which values of x are looked for.
%       * If S is a scalar, values are looked for in [1 : S].
%       * If S is empty or undefined, values are looked for in [1 : max(x)].
%   - inverse: flag for applying the reverse mapping (default = false).
%
% OUT:
%   - v: |S| x N indicator vectors such that
%       v(i, t) = 1 if S(i) is in x(:,t)
%       v(i, t) = 0 otherwise
%
% /////////////////////////////////////////////////////////////////////////

% check input parameters
% =========================================================================
% default to direct mapping
if nargin < 3
    inverse = false;
end

% guess superset
if nargin < 2 || isempty (S)
    if inverse
        S = size (x,1);
    else
        S = max (VBA_vec (x));
    end
end
if  isscalar (S)
    S = 1 : S;
end
S = VBA_vec (S);

% check input consistency
if ~ inverse
    assert( all (VBA_vec (ismember (x, S))), '*** VBA_indicator: values cannot be outside S') 
else
    assert(VBA_isBinary (x(~isnan(x))), '*** VBA_indicator: indicator vector should be binary');  
end

% apply mapping
% =========================================================================



if ~ inverse
% value -> indicator
% -------------------------------------------------------------------------
    [K, N] = size (x);

    v = zeros (numel (S), N);
    for k = 1 : K
        v = max(v, +bsxfun(@eq, x(k, :), S));
    end

    % reject ill-specified values
    v(:, any (isnan (x))) = NaN;
else
% indicator -> value
% -------------------------------------------------------------------------
    x = mat2cell(x, size(x,1), ones(1, size(x,2)));
    v = cellfun(@(xt) findInS (xt, S), x, 'UniformOutput', false);

    % if number of values are consistent, transform to matrix
    try
        v = cat(2, v{:});
    end 
end  
end

function idx = findInS (x, S)
    if any (isnan (x))
    	idx = nan;
    else
        idx = S(logical (x));
    end
end
