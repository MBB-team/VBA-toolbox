function analysis = demo_sensitivity_analyser(results)


out=[];
for i=1:length(results)
    % prediction errors are...
   
    try
    error_W(i,:) = abs(results(i).theta' - results(i).muTheta_W') ;
    error_WO(i,:) = abs(results(i).theta' - results(i).muTheta_WO') ;
    catch
        out = [out i];
    end
     if any(results(i).muTheta_WO(1:6)<0)
        out = [out i];
    end

end
error_W(out,:)=[];
error_WO(out,:)=[];

cerror_W = [error_W(:,[1 3 5]) ; error_W(:,[2 4 6])];
cerror_WO = [error_WO(:,[1 3 5]) ; error_WO(:,[2 4 6])];

results(out)=[];
thetas = horzcat(results(:).theta)';

% error covariance
analysis.improvement.med = median(error_WO - error_W);%./mean(abs(thetas)) ;
analysis.improvement.mea = mean(error_WO - error_W);%./mean(abs(thetas)) ;
analysis.improvement.cmed = median(cerror_WO - cerror_W);%./mean(abs(thetas)) ;
analysis.improvement.cmea = mean(cerror_WO - cerror_W);%./mean(abs(thetas)) ;
analysis.error_WO = error_WO;
analysis.error_W =  error_W;
analysis.cerror_WO = cerror_WO;
analysis.cerror_W =  cerror_W;

analysis.w.corr = corrcov(cov(error_W)) ;
analysis.wo.corr = corrcov(cov(error_WO)) ;


sigma_w = zeros(6);
for t=1:length(results)
    sigma_w = sigma_w + (results(t).SigmaTheta_W(1:6,1:6));
    sigma_wo = sigma_w + (results(t).SigmaTheta_WO(1:6,1:6));
end
analysis.w.sigma = sigma_w/length(results);
analysis.wo.sigma = sigma_wo/length(results);

analysis.improvement.cov2 = analysis.wo.sigma-analysis.w.sigma;
analysis.improvement.cov2(1:7:end) = NaN;

% connection strenght effect
% for t=1:size(thetas,2)
%     r_w=corrcoef(abs(thetas(:,t)),error_W(:,t));
%     r_wo=corrcoef(abs(thetas(:,t)),error_WO(:,t));
%     analysis.w.ampl(t) = r_w(2);
%     analysis.wo.ampl(t) = r_wo(2);
% end

% effect strenght effect
% contrasts = vertcat(results(:).effect) ;

% for t=1:size(thetas,2)
%     r_w=corrcoef(contrasts(:,3),error_W(:,t));
%     r_wo=corrcoef(contrasts(:,3),error_WO(:,t));
%     analysis.w.contrast(t) = r_w(2);
%     analysis.wo.contrast(t) = r_wo(2);
% end



end