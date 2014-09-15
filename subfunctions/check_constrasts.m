function [betas]=check_constrasts(y,SPM)

regId = vertcat(SPM.Sess.col);
regId(:,end-5:end) = [];
[nSession, nRegress] = size(regId);

designMat = [SPM.xX.X(:,vec(regId)) SPM.xX.X(:,(end-nSession+1):end)];

for i=1:size(y,1)
    b = glmfit(designMat,y(i,:)','normal','Constant','off'); 
    b((end-nSession+1):end)=[];
    betas{i} = reshape(b,nSession,nRegress) ;
end
    
