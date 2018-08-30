function [PI,out] = MixtureOfBinomials(y,K,options)

% This function inverts a binomials mixture model (BMM).
% [PI,out] = MixtureOfBinomials(y,K,options)
%
% IN:
%   - y: the Dxn binary data matrix, where n is the number of samples, and
%   D the dimension of the data.
%   - K: the number of classes to be used. As a general advice, K should be
%   significantly lower than both D and n.
%   - options: a optional structure (default = []), that contains:
%       .algo: flag for either 'VB' or 'Gibbs' variant (default = 'VB').
%       .maxIt: the maximum number of iterations of the algo
%       .TolFun: the minimum increase in free energy (for VB variant)
%       .b: a Kx1 vector, containing the prior sample counts for the
%       K classes (default = ones(K,1))
%       .c: a 2x1 vector, containing the prior sample counts for a possible
%       bias in favour of 1 versus 0 elements in the data matrix
% OUT:
%   - PI: a nxK matrix containing the estimated probability, for each data
%   point, of belonging to each of the classes.
%   - out: a structure containing the following entries:
%       .logEv: the approximation to the model log-evidence. This is
%       supposed to be a lower bound for VB, and an upper bound for the
%       Gibbs variant of the algorithm
%       .Elambda: a DxK matrix containing the estimated average pattern of
%       each class
%       .Ealpha: a Kx1 vector containing the estimated frequency of each
%       class within the sampled data



[D,n] = size(y);

try
    algo = options.algo;
catch
    algo = 'VB';
end

try
    maxIt = options.maxIt;
catch
    switch algo
        case 'VB'
            maxIt = 128;
        case 'Gibbs'
            maxIt = 1e4;
    end
end
try
    minIt = options.minIt;
catch
    minIt = 8;
end

try
    TolFun = options.TolFun;
catch
    TolFun = 1e-8;
end

try
    b = options.b;
catch
    b = ones(K,1);
end

try
    c = options.c;
catch
    c = ones(2,1);
end


out = [];

% Initialize
switch algo

    case 'Gibbs'

        addpath([fileparts(mfilename('fullpath')),filesep,'..',filesep,'sampling']);
        
        pi = zeros(n,K);
        gamma1 = ones(D,K);
        gamma2 = ones(D,K);
        f = ones(K,1);
        lambda = 0.5*ones(D,K);
        xi = zeros(n,K);
        PI = zeros(n,K);
        ALPHA = zeros(K,maxIt);
        LAMBDA = zeros(D,K,maxIt);
        L = zeros(maxIt,1);

        % draw sample from prior on class frequencies
        [alpha,seed] = dirichlet_sample (K,b,1e3);
        alpha = regularize(alpha);

        h = waitbar(0,'Gibbs sampling: please wait...');
        

    case 'VB'

                
        gamma1 = c(1).*(1 + 1e-3*rand(D,K));
        gamma2 = c(2).*(1 + 1e-3*rand(D,K));
        
        f = b;
        
        handles.f = figure;
        handles.a(1) = subplot(2,1,1,'parent',handles.f);
        handles.a(2) = subplot(2,1,2,'parent',handles.f);
        set(handles.a,'nextplot','add');
        title(handles.a(1),'Free energy')
        title(handles.a(2),'Derivative of the Free energy')

        % initialize free energy with priors
        lnPi = -y'*ones(D,K) - (1-y)'*ones(D,K) ...
            + repmat(psi(f)' - psi(sum(f)),n,1);
        F(1) = sum(sum(0.5.*lnPi)) - n*K*0.5*log(0.5);
        plot(handles.a(1),1,F(1),'o');
        drawnow

end


stop = 0;
it = 1;
while ~stop

    switch algo

        case 'Gibbs'

            % update class labels
            for i=1:n
                % get sufficient statistics
                lnPi = log(lambda')*y(:,i) + ...
                    log(1-lambda')*(1-y(:,i)) + ...
                    log(alpha);
                lnPi = lnPi' - max(lnPi);
                pi(i,:) = exp(lnPi)./sum(exp(lnPi));
                % draw a sample
                [tmp,seed] = multinomial_sample(1,K,pi(i,:)',seed);
                xi(i,:) = tmp';
            end

            % update class average
            for k=1:K
                for d=1:D
                    % Get sufficient statistics
                    gamma1(d,k) = c(1) + y(d,:)*xi(:,k);
                    gamma2(d,k) = c(2) + (1-y(d,:))*xi(:,k);
                    % draw a sample
                    [dir,seed] = dirichlet_sample(2,[gamma1(d,k),gamma2(d,k)],seed);
                    lambda(d,k) = dir(1);
                end
            end
            lambda = regularize(lambda);

            % update class frequencies
            % get sufficient statistics
            f = sum(xi,1)' + b;
            % draw sample
            [alpha,seed] = dirichlet_sample(K,f,seed);
            alpha = regularize(alpha);

            % calculate likelihood
            ll = logLikelihood(y,xi,lambda);
            
            
            % store samples
            PI = (it*PI + xi)./(it+1);
            LAMBDA(:,:,it) = lambda;
            ALPHA(:,it) = alpha;
            L(it) = exp(ll);
            
            
            if it==maxIt
                stop = 1;
                out.logEv = log(mean(L));
                out.Elambda = mean(LAMBDA,3);
                out.Ealpha =  mean(ALPHA,3);
            elseif it/100==floor(it./100)
                waitbar(it/maxIt,h)
            end


        case 'VB'

            % update class labels
            ElogLambda = psi(gamma1) - psi(gamma1+gamma2);
            ElogOneMinusLambda = psi(gamma2) - psi(gamma1+gamma2);
            ElogAlpha = psi(f) - psi(sum(f));
            lnPi = y'*ElogLambda + (1-y)'*ElogOneMinusLambda ...
                + repmat(ElogAlpha',n,1);
            maxlp = max(lnPi,[],2);
            pi = exp(lnPi-repmat(maxlp,1,K));
            sumPi = sum(pi,2);
            pi = pi./repmat(sumPi,1,K);
            pi = regularize(pi);
            

            % update class average
            for k=1:K
                for d=1:D
                    gamma1(d,k) = c(1) + y(d,:)*pi(:,k);
                    gamma2(d,k) = c(2) + (1-y(d,:))*pi(:,k);
                end
            end        
            

            % update class frequencies
            f = sum(pi,1)' + b;


            % calculatate free energy
            ElogLambda = psi(gamma1) - psi(gamma1+gamma2);
            ElogOneMinusLambda = psi(gamma2) - psi(gamma1+gamma2);
            ElogAlpha = psi(f) - psi(sum(f));
            lnPi = y'*ElogLambda + (1-y)'*ElogOneMinusLambda ...
                + repmat(ElogAlpha',n,1);
            F(it+1) = sum(sum(pi.*lnPi)) ...
                + gammaln(sum(b)) - sum(gammaln(b)) + sum((b-1)'*ElogAlpha) ...
                + sum(sum( gammaln(sum(c)) - gammaln(c(1)) - gammaln(c(2)) ...
                + (c(1)-1).*ElogLambda + (c(2)-1).*ElogOneMinusLambda )) ...
                - sum(sum(pi.*log(pi))) ...
                - gammaln(sum(f)) + sum(gammaln(f)) - sum((b-1)'*ElogAlpha) ...
                - sum(sum( gammaln(gamma1+gamma2) - gammaln(gamma1) - gammaln(gamma2) ...
                + (gamma1-1).*ElogLambda + (gamma2-1).*ElogOneMinusLambda ));

            % convergence criterion
            dF = F(end) - F(end-1);
            if ( dF <= TolFun || it == maxIt ) && it > minIt
                stop = 1;
                PI = pi;
                out.logEv = F(end);
                out.Elambda = gamma1./(gamma1+gamma2);
                out.Ealpha = f./sum(f);
                close(handles.f)
            else
                try
                    plot(handles.a(1),it+1,F(end),'o');
                    plot(handles.a(2),it,dF,'o');
                    drawnow
                end
            end


        otherwise
            error('Unknown algorithm')
            return
    end

    it = it + 1;


end

try close(h); end
return


%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%

function logL = logLikelihood(y,xi,lambda)
tmp = y'*log(lambda) + (1-y')*log(1-lambda);
logL = sum(xi(:).*tmp(:));




function [f] = psi(z)

siz = size(z);
z=z(:);
zz=z;
f = 0.*z; % reserve space in advance
%reflection point
p=find(real(z)<0.5);
if ~isempty(p)
    z(p)=1-z(p);
end
%Lanczos approximation for the complex plane
g=607/128; % best results when 4<=g<=5

c = [  0.99999999999999709182;
    57.156235665862923517;
    -59.597960355475491248;
    14.136097974741747174;
    -0.49191381609762019978;
    .33994649984811888699e-4;
    .46523628927048575665e-4;
    -.98374475304879564677e-4;
    .15808870322491248884e-3;
    -.21026444172410488319e-3;
    .21743961811521264320e-3;
    -.16431810653676389022e-3;
    .84418223983852743293e-4;
    -.26190838401581408670e-4;
    .36899182659531622704e-5];
n=0;
d=0;
for k=size(c,1):-1:2
    dz=1./(z+k-2);
    dd=c(k).*dz;
    d=d+dd;
    n=n-dd.*dz;
end
d=d+c(1);
gg=z+g-0.5;
%log is accurate to about 13 digits...
f = log(gg) + (n./d - g./gg) ;
if ~isempty(p)
    f(p) = f(p)-pi*cot(pi*zz(p));
end
p=find(round(zz)==zz & real(zz)<=0 & imag(zz)==0);
if ~isempty(p)
    f(p) = Inf;
end
f=reshape(f,siz);
return

function lambda = regularize(lambda)
i0 = find(lambda==0);
lambda(i0) = 1e-3;
i1 = find(lambda==1);
lambda(i1) = 1-1e-3;
