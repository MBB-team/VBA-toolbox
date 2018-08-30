function [pX,gX,pY,gY,X,Y] = get_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np,lx,ly)
% legacy code
s = warning ('on');
warning ('*** The function `get_MCMC_predictiveDensity` is now deprecated. Please use `VBA_MCMC_predictiveDensity` instead (same syntax).') 
warning (s);

% fallback
switch nargin
    case 10
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np,lx,ly);
    case 9
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np,lx);
    case 8
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np);
    case 7
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N);
    case 6
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim);
    otherwise
        error('VBA_MCMC_predictiveDensity: wrong number of arguments');
end



