function [DKL,DJS] = VB_KL(m1,v1,m2,v2,distrib)
% computes the KL and JS divergences between two distributions
% FORMAT [DKL] = VB_KL(m1,v1,m2,v2,disp)
% IN:
%   - m1,m2: the posterior means
%   - v1,v2: the posterior covariance matrices
%   - distrib: the type of distribution ({'Normal'} or 'Gamma')
% OUT:
%   - DKL: the KL-divergence D(p1||p2)
%   - DJS: the Jensen-Shannon divergence (symmetrized KL)

if isequal(m1,m2) && isequal(full(v1),full(v2))
    DKL = 0;
    DJS = 0;
    return
end

try, distrib; catch distrib='Normal'; end
DJS = [];

switch distrib

    case 'Normal'
        
        n = size(m1,1);
        iv2 = VB_inv(v2,[]);
        vv = v1*iv2;
        DKL = 0.5*(-VBA_logDet(vv) + trace(vv) + (m1-m2)'*iv2*(m1-m2) - n);
        
        if nargout > 1
            DJS = JensenShannon({m1,m2},{v1,v2},0.5*ones(2,1));
        end
        
    case 'Gamma'
        
        % derive standard parameters of the Gamma distribution
        b1 = m1/v1;
        b2 = m2/v2;
        a1 = b1*m1;
        a2 = b2*m2;
        
        DKL = a1*log(b1) - a2*log(b2) + gammaln(a2) - gammaln(a1) ...
            +(a1-a2)*(psi(a1)-log(b1)) - a1*(1-b2/b1);
        
        
    otherwise
        
        DKL = [];
        
end


function [DJS] = JensenShannon(mus,Qs,ps)
n = length(mus);
Vy = zeros(size(Qs{1}));
muy = zeros(size(mus{1}));
sH = 0;
for i=1:n
    muy = muy + ps(i).*mus{i};
    [e] = eig(full(Qs{i}));
    logDet = sum(log2(e));
    sH = sH + 0.5*ps(i).*logDet;
end
for i=1:n
    tmp = mus{i} - muy;
    tmp = tmp*tmp' + Qs{i};
    Vy = Vy + ps(i).*tmp;
end
[e] = eig(full(Vy));
Hy = 0.5*sum(log2(e));
DJS = Hy - sH;




function [f] = psi(z)
% psi(x) = d[log(gamma(x))]/dx
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






