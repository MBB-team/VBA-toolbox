function [HC,EHC,gridx] = RFT_Euler(X,verbose)
% derive induced Euler-characteristic
mi = min(X);
ma = max(X);
dm = (ma-mi).*1e-2;
gridx = mi:dm:ma;
[out0] = RFT_expectedTopo(gridx,L,fwhm,1);
EHC = out0.HC;
HC = NaN(length(gridx),1);
if verbose
    fprintf(1,'RFT: deriving Euler-Hadwiger characteristic ...')
    fprintf(1,'%6.2f %%',0)
end
for i=1:length(gridx)
    [out0] = RFT_clusters(X,gridx(i),0);
    HC(i) = length(out0.clusters);
    if verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',round(100*i/length(gridx)))
    end
end
if verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end
disp(' ')