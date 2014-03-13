function [SigmaX] = laggedCovX(SigmaX,muX,t,...
    dim,iR,indIn,iRp,sigmaHat,dG_dX,dF_dX,iQy,disp)

% computes lagged covariance of hidden states x(t) and x(t+1) 

Sinter = SigmaX.inter{t};

R = VBA_inv(iR{t},indIn{t});
if t == dim.n_t-1
    Rp = VBA_inv(iRp{t+1},indIn{t+1});
    Kt = Rp*dG_dX{t+1}*VBA_inv((VBA_inv(iQy{t+1},[])/sigmaHat ...
        + dG_dX{t+1}'*Rp*dG_dX{t+1}),indIn{t+1});
    SigmaX.inter{t} = (eye(dim.n) - Kt*dG_dX{t+1}')*dF_dX{t}'*R;
else
    Jt = R*dF_dX{t}*iRp{t+1};
    R2 = VBA_inv(iR{t+1},indIn{t+1});
    SigmaX.inter{t} = R2*Jt' + Jt*(SigmaX.inter{t+1}-...
        dF_dX{t}'*R2)*Jt';
end

SigmaX.inter{t} = SigmaX.inter{t}' + muX(:,t)*muX(:,t+1)';
% NB: SigmaX.inter{t} = E[(x(t)-mu(t))*(x(t+1)-mu(t+1))']

try
    if disp
        hff = figure;
        ha = subplot(2,2,1,'parent',hff);
        plot(ha,Sinter(:),SigmaX.inter{t}(:),'b.')
        grid(ha,'on')
        axis(ha,'tight')
        xlabel(ha,'previous')
        ylabel(ha,'actual')
        ha = subplot(2,2,2,'parent',hff);
        imagesc(Sinter-SigmaX.inter{t},'parent',ha),colorbar
        xlabel(ha,'previous-actual')
        ha = subplot(2,2,3,'parent',hff);
        imagesc(Sinter,'parent',ha),colorbar
        xlabel(ha,'previous')
        ha = subplot(2,2,4,'parent',hff);
        imagesc(SigmaX.inter{t},'parent',ha),colorbar
        xlabel(ha,'actual')
        pause
        close(hff)
    end
end