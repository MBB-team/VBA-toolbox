% In this script, I compute a few analysis on the raw data from Stefano
% - learning curve for both conditions mixed
% - learning curve for each condition individually


clear all
close all

%-----NOTES on how file is organized
% session1 = training session
% session2-3 = test sessions
%-- DATA
% column 1 = index of session
% column 2 = index of trial
% column 3 = condition (1=Gain (1€), 2=Neutral (0€), 3=Loss (-1€))
% column 4 = checktime
% column 5 = 1=Go, -1=NoGo
% column 6 = choice (1=correct, -1=incorrect)
% (correct = most rewarded option in cond=1 | less punished option in cond=3)
% colonne 7 = feedback (cond 1 : 1 = 1€ , -1=0€ | cond 3, 1=0€, -1=-1€)
% colonne 8 = reaction time, (millisec if Go, 0 if NoGo)
% colonne 9 = number of hesitations


%-------------- Loading data of each subject


controls1=[1001:1015 1103 1106 1108 1109]; % controles agés (Stefano, tache mêre)
SESSIONS = [2,3]; % the two test sessions, session 1 was the training session.
INVERSION = cell(1,length(controls1));

i_s =0;

for i_subject = controls1
    i_s = i_s+1;
    
    root = 'C:\Users\vincent.adam\Desktop\Vincent ADAM\Matlab\Demo_toolbox\Data\Data_tache_mere\VINCENT\'; % location of data folder
    subfolder = 'CONTROLS_1\';
    
    %-- LOADING THE TWO SESSIONS --%
    
    for i_session = SESSIONS; % index of session
        
        filename = ['ILtestSub',num2str(i_subject),'Session',num2str(i_session)]; % as saved by stefano
        file_root = [root, subfolder, filename]; % position of file on hard-drive
        
        load(file_root); % loading file
        
        cond = data(:,3); % vector of conditions
        i_cond1 = find(cond==1); % indices of trial where condition is 1
        i_cond2 = find(cond==2); % indices of trial where condition is 2
        
        yc1 = data(i_cond1,6)'>0; % binary choice (0: correct, 1:incorrect)
        yc2 = data(i_cond2,6)'>0; % binary choice (0: correct, 1:incorrect)
        
        Yc1((i_s-1)*2+(i_session-1),:)=yc1;
        Yc2((i_s-1)*2+(i_session-1),:)=yc2;
        
        
    end
    
    
end
    
    
    %-------------- Computing learning curves (both conditions mixed)
    
    LC = mean([Yc1;Yc2],1);
    sLC = std([Yc1;Yc2],1);
    
    %-------------- Computing learnig curves (each condition separately)
    
    LC1 = mean(Yc1,1);
    LC2 = mean(Yc2,1);
    sLC1 = std(Yc1,1);
    sLC2 = std(Yc2,1);
    
    %---------- Plotting results
    
    figure(1)
    hold on
    plot(LC,'k');   % plot(LC+sLC/2,'--k');    plot(LC-sLC/2,'--k')
    plot(LC1,'r');  %  plot(LC1+sLC1/2,'--r');    plot(LC1-sLC1/2,'--r')
    plot(LC2,'b');  %  plot(LC2+sLC2/2,'--b');    plot(LC2-sLC2/2,'--b')
    legend('Mixed conditions','Gain','Loss')
    title('Learning curves, all subjects together')
    ylabel('proportion of choice of the best alternative')
    
  