function results=VBA_lesionAnalysis_robust(posterior,out)

N=200;

posterior_temp = posterior;

results=VBA_lesionAnalysis(posterior,out);
    close all

for i=1:N-1
    posterior_temp.muTheta = posterior.muTheta + randn(numel(posterior.muTheta),1) .* sqrt(diag(posterior.SigmaTheta));
    results_temp=VBA_lesionAnalysis(posterior_temp,out);
    results = average_results(results,results_temp,N);
    close all
end

VBA_lesionAnalysisDisplay(results)
end


function r=average_results(r,temp,N)

r.normal=average_kernel(r.normal,temp.normal,N);
for i=1:numel(r.lesion)
    r.lesion(i)=average_kernel(r.lesion(i),temp.lesion(i),N);
end

end


function k=average_kernel(k,temp,N)

if isempty(k)
    k = temp;
else
    k.y = k.y + temp.y/N;
    for i=1:numel(k.kernel)
        k.kernel(i).params = k.kernel(i).params + temp.kernel(i).params/N;
        k.kernel(i).timeseries = k.kernel(i).timeseries + temp.kernel(i).timeseries/N;
        for j=1:numel(k.kernel(1).landmarks)
        k.kernel(i).landmarks(j).tMax = k.kernel(i).landmarks(j).tMax + temp.kernel(i).landmarks(j).tMax/N;
        k.kernel(i).landmarks(j).aMax = k.kernel(i).landmarks(j).aMax + temp.kernel(i).landmarks(j).aMax/N;
        end
        k.kernel(i).sigma.landmarks = k.kernel(i).sigma.landmarks + temp.kernel(i).sigma.landmarks/N;
        k.kernel(i).sigma.params = k.kernel(i).sigma.params + temp.kernel(i).sigma.params/N;
    end
    
end
end

