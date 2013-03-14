function [ROC,AUC]=demo_negfeedback_analyser(results)

noise_lvls = unique([results(:).noise]);
reps_lvls = [5 10 15] ;%unique([results(:).reps]);

models = {'1  0','0 -1','1 -1'} ; %unique(cellfun(@num2str,{results(:).A},'UniformOutput',false));
schemes = {'1  0','0  1'}; %unique(cellfun(@num2str,{results(:).hA},'UniformOutput',false));

get_model_id = @(r) find(cellfun(@(x) all(x==num2str(r.A)),models));
get_scheme_id = @(r) find(cellfun(@(x) all(x==num2str(r.hA)),schemes));

%% compute individual scores
scores = [];
out=[];
for t=1:length(results)
    
    scores(t).noise=results(t).noise;
    scores(t).reps=results(t).reps;

    try
        for lb={'w','wo'}
            n=lb{1};
            scores(t).(['model_' n]) = get_score(results(t).(n).Fnested,'model',get_model_id(results(t)));
            scores(t).(['scheme_' n]) =get_score(results(t).(n).Fnested,'scheme',get_scheme_id(results(t)));
        end
    catch
        out(end+1)=t;
    end
end
scores(out)=[];

%% construct ROC curves for each group
alphas=0:.5:500;
cpt=1;
for noise = noise_lvls
    for reps = reps_lvls
        temp_score = only(scores,'noise',noise);
        temp_score = only(temp_score,'reps',reps);
        
        for lb={'w','wo'}
            n=lb{1};
            sens_model.(n) = zeros(1,length(alphas));
            spec_model.(n) = zeros(1,length(alphas));
            sens_scheme.(n) = zeros(1,length(alphas));
            spec_scheme.(n) = zeros(1,length(alphas));  
        
            for ai=1:length(alphas)
                model_scores = [temp_score(:).(['model_' n])];
                model_TP = sum(model_scores >= alphas(ai));
                model_FP = sum(-model_scores >= alphas(ai)) ;
                sens_model.(n)(ai) =  model_TP/sum(model_scores>0);
                spec_model.(n)(ai) = model_FP/sum(model_scores<0);
                
                scheme_scores = [temp_score(:).(['scheme_' n])];
                scheme_TP = sum(scheme_scores >= alphas(ai));
                scheme_FP = sum(-scheme_scores >= alphas(ai)) ;
                sens_scheme.(n)(ai) =  scheme_TP/sum(scheme_scores>0);
                spec_scheme.(n)(ai) = scheme_FP/sum(scheme_scores<0);
            end
        end
        ROC(cpt).noise=noise;
        ROC(cpt).reps=reps;
        ROC(cpt).sens_model_w=sens_model.w;
        ROC(cpt).spec_model_w=spec_model.w;
        ROC(cpt).sens_model_wo=sens_model.wo;
        ROC(cpt).spec_model_wo=spec_model.wo;
        ROC(cpt).sens_scheme_w=sens_scheme.w;
        ROC(cpt).spec_scheme_w=spec_scheme.w;
        ROC(cpt).sens_scheme_wo=sens_scheme.wo;
        ROC(cpt).spec_scheme_wo=spec_scheme.wo;
        cpt=cpt+1;
    end
end

%% compute area under curve
cpt=1;
for noise = noise_lvls
    for reps = reps_lvls
        AUC(cpt).noise=noise;
        AUC(cpt).reps=reps;
        AUC(cpt).model_w = -trapz(ROC(cpt).spec_model_w,ROC(cpt).sens_model_w);
        AUC(cpt).model_wo = -trapz(ROC(cpt).spec_model_wo,ROC(cpt).sens_model_wo);
        AUC(cpt).scheme_w = -trapz(ROC(cpt).spec_scheme_w,ROC(cpt).sens_scheme_w);
        AUC(cpt).scheme_wo = -trapz(ROC(cpt).spec_scheme_wo,ROC(cpt).sens_scheme_wo);
        cpt=cpt+1;
    end
end

end

function so=only(s,f,id)
   is_id = @(x) strcmp(x,num2str(id)) ;
   idx = cellfun(is_id,cellfun(@num2str,{s.(f)},'UniformOutput',false));
   so=s(idx);
end

function [score]=get_score(Fs,type,id)
    % return score for the model/scheme (depending on type) id.

    if     strcmp(type,'scheme')
        Ffam = sum(Fs,2);
    elseif strcmp(type,'model')
        Ffam = sum(Fs,1);
    else
        error('***');
    end
    
    win  = find(Ffam==max(Ffam),1);
    not_win = setdiff(1:length(Ffam),win);

    score = Ffam(win) - max(Ffam(not_win));
    if win~=id
        score = -score;
    end
    
  
  
end