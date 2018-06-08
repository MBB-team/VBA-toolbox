function uu = VBA_getU(u,options,dim,flag)
% reshape micro-resolution driving input u
% function uu = VBA_getU(u,options,dim,flag)
% IN:
%   - u: macro-/micro- resolution input
%   - options: the options structure (see VBA_check.m)
%   - dim: the dim structure (see VBA_check.m)
%   - flag: gives the transformation direction (from micro- to macro- or
%   back...)
% OUT:
%   - uu: reshaped micro-resolution driving input

if options.microU && ~isequal(options.decim,1)
    switch flag
        case '2macro'
            u0 = u;
            nut = dim.u*options.decim;
            uu = zeros(nut,dim.n_t);
            for t=1:dim.n_t
                uu(:,t) = VBA_vec(u0(:,(t-1)*options.decim+1:t*options.decim));
            end
        case 'back2micro'
            u0 = u;
            uu = zeros(dim.u,options.decim*dim.n_t);
            for t=1:dim.n_t
                ut = reshape(u0(:,t),dim.u,options.decim);
                uu(:,(t-1)*options.decim+1:t*options.decim) = ut;
            end
    end
else % do not change input
    uu = u;
end


