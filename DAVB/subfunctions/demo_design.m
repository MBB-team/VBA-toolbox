% demo for design optimization

dt = 1e-1;
n = 5*60/dt;
TR = 3;
N = 32;
ncon = 8;


ind = 1:floor(TR/dt):(n-floor(TR/dt));
ni = length(ind);
gridt = [1:1:ni]*floor(TR/dt);

NFFT = 2^nextpow2(n); 
gridf = 1/(2*dt)*linspace(0,1,NFFT/2+1);
% tf = find(gridf>1/8);
tf = 2;%tf(1);


mo = zeros(length(gridt),1);
so = mo;
mn = mo;
sn = mo;
mmf = mo;
smf = mo;
mcf = mo;
scf = mo;

[hrf] = spm_hrf(dt);
[X0] = get_U_basis(n*dt,dt,ncon,'Fourier_complete')';
ncon = size(X0,2);
X0 = X0(ind,:);

A = -eye(ni) + diag(ones(ni-1,1),1);
A(end,:) = [];
iA = pinv(A);
S = iA*iA' + eye(ni);
iS = pinv(S);

% figure,imagesc(iS),colorbar,pause

c = [1;0;zeros(ncon-1,1)];


for i=1:length(gridt)
    
    opt = zeros(N,1);
    nn = opt;
    mf = opt;
    cf = opt;
    
    for j=1:N
        
        X1 = zeros(n,1);
        i1 = 1:gridt(i):n;
        i1 = floor(i1 + 1e-1*gridt(i)*randn(1,length(i1)));
        i1(i1<=0) = [];
        i1(i1>n) = [];
        X1(i1) = 1;
        
        nn(j) = sum(X1>0);
        
        X1 = conv(double(X1),hrf,'same');
        
        Y = fft(X1,NFFT)/n;
        py = abs(Y(1:NFFT/2+1));
        py = py./sum(py);
        cf(j) = gridf(tf:end)*py(tf:end);
        [mp,im] = max(py(tf:end));
        mf(j) = gridf(im+tf-1);
        
        X = X1(ind,:);
        X = [X,X0];
        Cov = inv(X'*iS*X + 1e-8*eye(ncon+1));
        opt(j) = 1/trace(c'*Cov*c);
        
    end
    
    mo(i) = mean(opt);
    so(i) = std(opt)./sqrt(N);
    
    mn(i) = mean(nn);
    sn(i) = std(nn)./sqrt(N);
    
    mmf(i) = mean(mf);
    smf(i) = std(mf)./sqrt(N);
    
    mcf(i) = mean(cf);
    scf(i) = std(cf)./sqrt(N);
    
    
end


hf = figure;
subplot(2,2,1),errorbar(gridt,mo,so),title('design efficiency')
subplot(2,2,2),errorbar(gridt,mn,sn),title('# events')
subplot(2,2,3),plot(mn,mo,'.'),xlabel('# events'),ylabel('design efficiency')
subplot(2,2,4),errorbar(gridt,mmf,smf),title('design max frequency')

