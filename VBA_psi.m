function [f] = VBA_psi(z)
% psi(x) = d[log(gamma(x))]/dx

try
    f = psi(z+1)-1./z; % for numerical purposes
catch
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
end
return



