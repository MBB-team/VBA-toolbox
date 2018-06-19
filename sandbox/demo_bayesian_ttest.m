%% demo_bayesian_ttest
% this script demonstrates the use of the 'bayesian_ttest' function for 2
% random normal samples of equal variances.
%

% initialize
clear all;
clc;
rng('default'); % reproducible random-number generation

% parameters
sigma = 1;
nsample = [10 30 100];
d = [0.1 0.5 1];
n_iteration = 1e1; % MUST be increased to a larger number (eg 1e3) for real test
n_null = numel(nsample)*numel(n_iteration);
n_alternative = n_null*numel(d);

% variable definitions
report=struct;
% - null hypothesis
report.null.nsample = nan(n_null,1);
report.null.h_student = nan(n_null,1);
report.null.h_bayes = nan(n_null,1);
report.null.ep_student = nan(n_null,1);
report.null.ep_bayes = nan(n_null,1);

% - alternative hypothesis
report.alternative.nsample = nan(n_alternative,1);
report.alternative.d = nan(n_alternative,1);
report.alternative.h_student = nan(n_alternative,1);
report.alternative.h_bayes = nan(n_alternative,1);
report.alternative.ep_student = nan(n_alternative,1);
report.alternative.ep_bayes = nan(n_alternative,1);


%% iterations
for iteration = 1:n_iteration
    
    % display
    fprintf('# iteration %d/%d \n',iteration,n_iteration);

    for i_nsample = 1:numel(nsample)
        
        % indexingscatte
        index = i_nsample ...
                + (iteration-1)*numel(nsample) ;
       
        % H0 - no effect
        % -----------------------------------------

        %%% random sampling
        report.null.nsample(index) = nsample(i_nsample);
        x = randn(report.null.nsample(index),2).*sigma ;

        %%% frequentist student-test
        [h,p] = ttest2(x(:,1),x(:,2));
        report.null.h_student(index) = h;
        report.null.ep_student(index) = 1-p;

        %%% bayesian regularised test
        [h,ep] = bayesian_ttest(x(:,1),x(:,2));
        report.null.h_bayes(index) = h;
        report.null.ep_bayes(index) = ep;

        % H1 - effect size = d
        % -----------------------------------------
        for i_effect_size = 1:numel(d)
            
            % indexing
            index = i_effect_size ...
                  + (i_nsample-1)*numel(d)...
                  + (iteration-1)*numel(d)*numel(nsample) ;

            %%% random sampling
            report.alternative.nsample(index) = nsample(i_nsample);
            report.alternative.d(index) = d(i_effect_size);
            x = randn(report.alternative.nsample(index),2).*sigma ...
                + repmat([0 report.alternative.d(index)],report.alternative.nsample(index),1);

            %%% frequentist student-test
            [h,p] = ttest2(x(:,1),x(:,2));
            report.alternative.h_student(index) = h;
            report.alternative.ep_student(index) = 1-p;

            %%% bayesian regularised test
            [h,ep] = bayesian_ttest(x(:,1),x(:,2));
            report.alternative.h_bayes(index) = h;
            report.alternative.ep_bayes(index) = ep;
        end
    end

end

% display reports
disp(report.null);
disp(report.alternative);

%% graphics
set(groot,'defaultLineLineWidth',1.5);
scaling = 0.75;
f=figure;
%%% false positive rate
subplot(2,2,1);hold on;
x_bin = findgroups(report.null.nsample);
x_mean = splitapply(@VBA_nanmean,report.null.nsample,x_bin);
y_mean = splitapply(@VBA_nanmean,report.null.h_student,x_bin);
plot(x_mean,y_mean,'b');
y_mean = splitapply(@VBA_nanmean,report.null.h_bayes,x_bin);
plot(x_mean,y_mean,'r');
legend('student','bayes');
xlabel(' sample size (n) ');
ylabel(' false positive rate (%) ');
title(' H0: null hypothesis ');
ylim([0 1]);

%%% true positive rate
subplot(2,2,2);hold on;
x_bin = findgroups(report.alternative.nsample);
x_mean = splitapply(@VBA_nanmean,report.alternative.nsample,x_bin);
for i_effect_size = 1:numel(d)
    subset = (report.alternative.d == d(i_effect_size));
    y_mean = splitapply(@VBA_nanmean,report.alternative.h_student(subset),x_bin(subset));
    p(i_effect_size) = plot(x_mean,y_mean,'b','LineWidth',scaling*i_effect_size);
    y_mean = splitapply(@VBA_nanmean,report.alternative.h_bayes(subset),x_bin(subset));
    plot(x_mean,y_mean,'r','LineWidth',scaling*i_effect_size);
end
legend([p(1) p(2) p(3)],{'d=0.1','d=0.5','d=1'});
xlabel(' sample size (n) ');
ylabel(' true positive rate (%) ');
title(' H1: alternative hypothesis ');
ylim([0 1]);

%%% confidence
subplot(2,2,3);hold on;
x_bin = findgroups(report.null.nsample);
x_mean = splitapply(@VBA_nanmean,report.null.nsample,x_bin);
y_mean = splitapply(@VBA_nanmean,report.null.ep_student,x_bin);
plot(x_mean,y_mean,'b');
y_mean = splitapply(@VBA_nanmean,report.null.ep_bayes,x_bin);
plot(x_mean,y_mean,'r');
xlabel(' sample size (n) ');
ylabel(' confidence = 1 - p (%) ');
ylim([0 1]);

subplot(2,2,4);hold on;
x_bin = findgroups(report.alternative.nsample);
x_mean = splitapply(@VBA_nanmean,report.alternative.nsample,x_bin);
for i_effect_size = 1:numel(d)
    subset = (report.alternative.d == d(i_effect_size));
    y_mean = splitapply(@VBA_nanmean,report.alternative.ep_student(subset),x_bin(subset));
    plot(x_mean,y_mean,'b','LineWidth',scaling*i_effect_size);
    y_mean = splitapply(@VBA_nanmean,report.alternative.ep_bayes(subset),x_bin(subset));
    plot(x_mean,y_mean,'r','LineWidth',scaling*i_effect_size);
end
xlabel(' sample size (n) ');
ylabel(' confidence = 1 - p (%) ');
ylim([0 1]);

%% save
% save('demo_bayesian_ttest','report');
