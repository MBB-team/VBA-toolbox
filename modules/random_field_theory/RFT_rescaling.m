function lnu = RFT_rescaling(nu)
if nu <3
    disp('Error: RFT smoothness rescaling is not defined for d.o.f. <3!')
    lnu = [];
end
dt = 1e-2;
gt = -50:dt:50;
ft = (gt.^2 + nu-1).^2.*tpdf(gt,nu).^3./myp(gt,nu).^2;
ft = ft./((nu-1).*(nu-2));
ft = ft(~isinf(ft));
gt = gt(~isinf(ft));
lnu = sum(ft.*dt);


function p = myp(t,nu)
p = VBA_spm_Npdf(spm_invNcdf(spm_Tcdf(t,nu),0,1));