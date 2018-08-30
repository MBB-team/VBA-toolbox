function [pX,gX,pY,gY,X,Y,U] = get_MCMC_predictiveDensity_fb(f_fname,g_fname,u,n_t,options,dim,fb,N,np,lx,ly)

% legacy code
s = warning ('on');
warning ('*** The function `get_MCMC_predictiveDensity_fb` is now deprecated. Please use `VBA_MCMC_predictiveDensity_fb` instead (same syntax).') 
warning (s);

% fallback
switch nargin
    case 11
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity_fb(f_fname,g_fname,u,n_t,options,dim,fb,N,np,lx,ly);
    case 10
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity_fb(f_fname,g_fname,u,n_t,options,dim,fb,N,np,lx);
    case 9
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity_fb(f_fname,g_fname,u,n_t,options,dim,fb,N,np);
    case 8
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity_fb(f_fname,g_fname,u,n_t,options,dim,fb,N);
    case 7
        [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity_fb(f_fname,g_fname,u,n_t,options,dim,fb);
    otherwise
        error('VBA_MCMC_predictiveDensity_fb: wrong number of arguments');
end



