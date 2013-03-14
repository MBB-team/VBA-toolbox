function [results] = launcher(N,As,Bs,hAs,noise,reps)


parfor_progress(N);
parfor i=1:N
    try
        [w,wo]=demo_negfeedback(As,Bs,hAs,noise,reps); 
    catch
        w=NaN;
        wo=NaN;
    end
     results(i).noise=noise;
     results(i).reps=reps;
     results(i).A=As;
     results(i).B=Bs;
     results(i).hA=hAs;
     results(i).w=w;
     results(i).wo=wo;
     parfor_progress;
end

parfor_progress(0);

end
