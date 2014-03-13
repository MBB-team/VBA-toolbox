% OTO : loop over subjects

ana = 'load';
full_results = 1;
n = 19;
% flag = 1;
for i=1:n
    for j=1:2
        for flag = 1:2
            
            %         fname{1} = ['OTO_AVL_',num2str(i),'_cue',num2str(j),'.mat'];
            %         fname{1} = ['OTO_AVL_',num2str(i),'_cue',num2str(j),'_v2_static.mat'];
            %         fname{1} = ['OTO_AVL_',num2str(i),'_cue',num2str(j),'_v2.mat'];
            %         fname{1} = ['OTO_AVL_',num2str(i),'_cue',num2str(j),'_v3.mat'];
            %         fname{1} = ['OTO_AVL_',num2str(i),'_cue',num2str(j),'_v4.mat'];
            fname{1} = [...
                'OTO_AVL_',num2str(i),...
                '_cue',num2str(j),...
                '_flag',num2str(flag),'_v6.mat'];
            switch ana
                
                case 'load'
                    % load analysis file
                    for f=1:length(fname)
                        load(fname{f})
                        FE(f,i,j) = out.F;
                    end
                    [tmp,bm(i,j)] = max(FE(:,i,j));
                    load(fname{bm(i,j)})
                    % store estimated parameters
                    theta(:,i,j) = [posterior.muTheta(2);posterior.muPhi];
                    vp(:,i,j) = full([posterior.SigmaTheta(2,2);...
                        diag(posterior.SigmaPhi)]);
                    F(i,j,flag) = out.F;
                    % re-display VB analysis output
                    if full_results & flag==2
                        sname = ...
                            ['subject #',num2str(i),...
                            ' ; cue # ',num2str(j),...
                            ' ; flag # ',num2str(flag)];
%                         hf(j) = VBA_ReDisplay(posterior,out);
%                         set(hf(j),'tag','df','name',sname)
                        % display posterior risk and RT histograms
                        [Q,ha] = risk_OTO(posterior,out);
                        set(get(ha(1),'parent'),...
                            'name',sname)
                    end
                    
                case 'do'
                    fname = fname{1};
                    [y,u,iQy,nmissed(i,j)] = read_data_AVL('data',i,j);
                    disp('----------')
                    if ~exist(fname,'file')
                        % do and save analysis
                        disp(['analysing ',fname])
                        save(fname,'fname')
                        [posterior,out] = ana_behav_AVL(y,u,flag,iQy);
                        save(fname,'posterior','out')
                    else
                        % skip and go to next data file
                        disp(['skipping ',fname])
                    end
                    disp('----------')
                    clear fname
            end
        end
    end
    if isequal(ana,'load') && full_results
        disp('---paused!---')
        pause
        try; close(hf); end
    end
end