function [mw,vw] = extractKernels(posterior,out)
inG = out.options.inG;
mw = zeros(inG.dim.n_t,inG.dim.nu);
vw = zeros(inG.dim.n_t,inG.dim.nu);
for i=1:inG.dim.nu
    ind{i} = (i-1)*inG.dim.n_t+1:i*inG.dim.n_t;
    mw(:,i) = posterior.muPhi(ind{i});
    vw(:,i) = diag(posterior.SigmaPhi(ind{i},ind{i}));
end
