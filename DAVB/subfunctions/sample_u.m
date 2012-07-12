function [U,fU,gridf] = sample_u(nu,iti,f0,dt,T,dtu,stdu)
% sample input with random jitter and permutation of conditions
% [U,futh] = sample_u(nu,iti,f0,dt,T,dtu,stdu)
% IN:
%   - nu: number of conditions
%   - iti: mean inter-stimuli interval
%   - f0: ratio of null-events (must be 0<f0<1)
%   - dt: micro-time resolution for input train
%   - T: total time of the input train
%   - dtu: event duration
%   - stdu: standard-deviation of temporal jitter
% OUT:
%   - U: nuX(T/dt) train of stimuli
%   - fU: Fourier transform of U
%   - gridf: the frequency grid over which the Fourier transform was
%   derived



try dtu; catch; dtu = 1; end
try stdu; catch; stdu = 1; end

n_t = ceil(T./dt);
uth = zeros(1,n_t);
ts0 = ceil(iti./dt):floor((iti+dtu)./dt):n_t;
ts = floor(ts0 + (stdu/dt)*rand(size(ts0)));
ts(ts<=0) = [];
ts(ts>=n_t) = [];
uth(ts) = 1;

[C] = get_nkdraws(2,nu,0)-1; % get all combinations of categorical events
ntypes = size(C,2)-1; % null events!
C = C(:,2:ntypes+1);
C = C(:,randperm(ntypes));

indt = find(uth==1);
nevents = length(indt);
n_nulls = floor(nevents.*f0);
C = [repmat(C,1,ceil((1-f0).*nevents./ntypes)),zeros(nu,n_nulls)];
C = C(:,randperm(size(C,2)));
C = C(:,1:nevents);

U = zeros(nu,n_t);
U(:,indt) = C;

n = ceil(dtu./dt);
if n >1
    for i=1:nu
        indt = find(U(i,:)==1);
        indt = repmat(indt,n,1) + repmat([0:n-1]',1,length(indt));
        U(i,indt) = 1;
    end
end


if nargout < 2
    return
end

Fs = 1./dt;                    % Sampling frequency
NFFT = 2^nextpow2(n_t); % Next power of 2 from length of y
Y = fft(U',NFFT)/n_t;
gridf = Fs/2*linspace(0,1,NFFT/2+1);
fU = 2*abs(Y(1:NFFT/2+1,:));
% Plot single-sided amplitude spectrum.



