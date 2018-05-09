function [u] = get_U_basis(T,dt,n,btype)
% build input basis function set on temporal grid
% function [u] = get_U_basis(T,dt,n,btype)
% IN:
%   - T: end time (in sec) {1e2}
%   - dt: temporal resolution (in sec) {1}
%   - n: number of basis functions {16}
%   - btype: 'Fourier' or 'RBF' {'Fourier'}
% OUT:
%   - u: nX(T/dt) matrix of basis functions evaluated on the temporal
%   grid [0:dt:floor(n_t*dt)].

try; T; catch; T = 1e2; end
try; dt; catch; dt = 1; end
try; n; catch; n = 16; end
try; btype; catch; btype = 'Fourier'; end

if n == 0
    u = [];
else
    ut = dt:dt:T;
    switch btype
        case 'Fourier'
            u_fname = @u_Fourier;
            in.W = 0:(n-1);
            in.T = T;
            in.deltat = dt;
        case 'Fourier_complete'
            u_fname = @u_FourierComplete;
            in.W = 0:(n-1);
            in.T = T;
            in.deltat = dt;
            n = 2*(n-1)+1;
        case 'RBF'
            u_fname = @u_RBF;
            in.centres = [(T/n)/2:(T/n):T];
            in.sig = (T/n);
        otherwise
            u = [];
            return
    end
    [u] = u_fname([],eye(n),ut,in);
end