%This MATLAB code is associated with the following manuscript: Barrick,
%S.K., S.R. Clippinger, L. Greenberg, M.J. Greenberg. 2019. Computational
%tool to study perturbations in muscle regulation and its application to 
%heart disease.

%This is a script for global fitting of the data and estimation of
%uncertainties by bootstrapping.

%%

%Global fitting and uncertainty estimation

%This script performs global fitting and bootstrapping simulations of the
%pooled fluorescence titration data collected at three calcium concentrations:
%2 mM EGTA (nocal), pCa 3 (cal), and pCa 6.25 (midcal).

%NOTE: Each technical replicate should be normalized separately, using the
%script named "Script_normalization_replicate." These normalized curves
%should be combined into a pooled data set consisting of the individually
%normalized data for each technical replicate. This pooled data set is the
%input for this script.

%This assigns the measured values to variables for the rest of the
%program.  The user can change the name so their data is assigned to the proper
%variables.
x1=s1_nocal; %This is the myosin concentration for the data collected at low calcium (2 mM EGTA)
y1=fl_nocal; %This is the fractional change in fluorescence for the data collected at low calcium (2 mM EGTA)
x2=s1_cal; %This is the myosin concentration for the data collected at high calcium (pCa 3)
y2=fl_cal; %This is the fractional change in fluorescence for the data collected at high calcium (pCa 3)
x3=s1_midcal; %This is the myosin concentration for the data collected at intermediate calcium (pCa 6.25)
y3=fl_midcal; %This is the fractional change in fluorescence for the data collected at intermediate calcium (pCa 6.25)

%This assigns the variables appropriate values.
KBnocal=0.290; %User should change this to their value for KB from the stopped flow experiments
KS=18; %This value is fixed as described in McKillop and Geeves (1993).
KBcal=20; %This value is fixed as described in McKillop and Geeves (1993).

%These are the parameters for the confidence interval calculations.
N_bs=1000; %Number of bootstrapping simulations. 1000 is usually sufficient.
CI_bound=0.95; %This is the percent value for the confidence intervals.
tail=2; %This defines whether the p-value is obtained for a 1- or 2-tailed
%test. Use 2 by default.

%This removes any non-numerical values in the imported data set (for
%example, if you have blank spaces due to different numbers of data points
%for each calcium concentration).
q=find(isnan(y1));
x1(q)=[];
y1(q)=[];
q=find(isnan(y2));
x2(q)=[];
y2(q)=[];
q=find(isnan(y3));
x3(q)=[];
y3(q)=[];

%This clears the variables used during fitting to make it easier to
%consecutively perform multiple iterations of the fitting.
clear y_error* best_fit_params *_sq_min bs* bootstrap_params*

%This sets the options for the fitting.
hybridopts = optimoptions('patternsearch','TolFun',1e-8,'TolX',1e-7,'TolMesh',1e-7,'display','none','MaxIter',50000,'MaxFunEvals',60000,'TimeLimit',120);
options = optimoptions('simulannealbnd','TolFun',1e-8,'TimeLimit',10,'display','none','InitialTemperature',0.5,'MaxIter',100000); 
options = optimoptions('simulannealbnd',options,'HybridFcn',{@patternsearch,hybridopts});

%Note that the 'optimoptions' command is not compatible with some earlier
%versions of MATLAB. If you get an error related to this command, use the
%following code to set the options instead:
%hybridopts = psoptimset('TolFun',1e-8,'TolX',1e-7,'TolMesh',1e-7,'display','none','MaxIter',50000,'MaxFunEvals',60000,'TimeLimit',120);
%options = saoptimset('TolFun',1e-8,'TimeLimit',10,'display','none','InitialTemperature',0.5,'MaxIter',100000); 
%options = saoptimset(options,'HybridFcn',{@patternsearch,hybridopts});

%The following code performs the global fit via simulated annealing and
%least-squares minimization of y_error. It returns the set of best-fit parameters
%(best_fit_params) as well as the value of the sum of squares at the
%minimum (sum_sq_min). The order of parameters in best_fit_params is the
%same as in y_error. Fixed parameters are the values of KBnocal, KBcal, and
%KS. Fitted parameters are nH (the Hill coefficient), KW, KT at low
%(KTnocal), high (KTcal), and intermediate calcium (KTmidcal), and the
%normalization amplitudes (A,B,C) for the pooled data sets.

%This defines y_error for the nocal, cal, and midcal conditions as the sum
%of squared differences between Equation 1 (see user guide) and pooled
%fractional change in fluorescence data for that condition (y1, y2, y3).
y_error_nocal = @(X1,KW,KTnocal,nH,A) sum(((A.*KW.*x1.*(1 + KW.*x1.*(1 + KS)).^(nH - 1).*(KTnocal.*(1 + KS).^nH + 1))./((KTnocal.*(1 + KW.*x1.*(1 + KS)).^nH + (1 + KW.*x1).^nH + 1./KBnocal).*(1 + KS).^(nH - 1)) - y1).^2);
y_error_cal = @(X2,KW,KTcal,nH,B) sum(((B.*KW.*x2.*(1 + KW.*x2.*(1 + KS)).^(nH - 1).*(KTcal.*(1 + KS).^nH + 1))./((KTcal.*(1 + KW.*x2.*(1 + KS)).^nH + (1 + KW.*x2).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y2).^2);
y_error_midcal = @(X3,KW,KTmidcal,nH,C) sum(((C.*KW.*x3.*(1 + KW.*x3.*(1 + KS)).^(nH - 1).*(KTmidcal.*(1 + KS).^nH + 1))./((KTmidcal.*(1 + KW.*x3.*(1 + KS)).^nH + (1 + KW.*x3).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y3).^2);

%The error from all three curves is summed to create y_error, which is later
%minimized in the global fit.
y_error = @(X1,X2,X3,KW,KTnocal,KTcal,KTmidcal,nH,A,B,C) y_error_nocal(X1,KW,KTnocal,nH,A)+y_error_cal(X2,KW,KTcal,nH,B)+y_error_midcal(X3,KW,KTmidcal,nH,C);

%This performs the global fit by minimizing y_error. Here, we are interested
%in the best-fit parameter values.
[best_fit_params,sum_sq_min] = simulannealbnd(@(z) y_error(x1,x2,x3,z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8)),[0.03,0.12,0.24,0.18,7,1,1,1],[0,0,0,0,0,0,0,0],[10,10,10,10,20,10,10,10],options);

%This prints to the command window a table containing the best-fit values
%of the parameters floated during the fitting (best_fit_values), as well as
%the value of the sum of squares at the minimum (sum_sq_min).
paramIDs = {'KW','KTnocal','KTcal','KTmidcal','nH','A','B','C'};
best_fit_values = array2table(best_fit_params,'VariableNames',paramIDs)
sum_sq_min

%The following code performs the bootstrapping simulations for error estimation.
%It outputs a matrix, bootstrap_params, that contains the fit values
%for each of the simulated data sets. The columns give the values of the
%parameters and each row contains values obtained from a separate simulation.

%This randomly selects indices of the original, pooled data set to
%generate the simulated (i.e., bootstrapped) data sets.
bs_ind_nocal = ceil(length(x1)*rand(length(x1),N_bs));
bs_ind_cal = ceil(length(x2)*rand(length(x2),N_bs));
bs_ind_midcal = ceil(length(x3)*rand(length(x3),N_bs));

%This generates fits to the N_bs simulated data sets. Each row of
%bootstrap_fit_params contains the best-fit parameters for one round of
%resampling. The order of parameters in the columns of bootstrap_params is
%the same as in y_error.
parfor i=1:N_bs
    x_ind_nocal = x1(bs_ind_nocal(:,i));
    y_ind_nocal = y1(bs_ind_nocal(:,i));
    x_ind_cal = x2(bs_ind_cal(:,i));
    y_ind_cal = y2(bs_ind_cal(:,i));
    x_ind_midcal = x3(bs_ind_midcal(:,i));
    y_ind_midcal = y3(bs_ind_midcal(:,i));

y_error_nocal = @(X1,KW,KTnocal,nH,A) sum(((A.*KW.*X1.*(1 + KW.*X1.*(1 + KS)).^(nH - 1).*(KTnocal.*(1 + KS).^nH + 1))./((KTnocal.*(1 + KW.*X1.*(1 + KS)).^nH + (1 + KW.*X1).^nH + 1./KBnocal).*(1 + KS).^(nH - 1)) - y_ind_nocal).^2);
y_error_cal = @(X2,KW,KTcal,nH,B) sum(((B.*KW.*X2.*(1 + KW.*X2.*(1 + KS)).^(nH - 1).*(KTcal.*(1 + KS).^nH + 1))./((KTcal.*(1 + KW.*X2.*(1 + KS)).^nH + (1 + KW.*X2).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y_ind_cal).^2);
y_error_midcal = @(X3,KW,KTmidcal,nH,C) sum(((C.*KW.*X3.*(1 + KW.*X3.*(1 + KS)).^(nH - 1).*(KTmidcal.*(1 + KS).^nH + 1))./((KTmidcal.*(1 + KW.*X3.*(1 + KS)).^nH + (1 + KW.*X3).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y_ind_midcal).^2);
y_error = @(X1,X2,X3,KW,KTnocal,KTcal,KTmidcal,nH,A,B,C) y_error_nocal(X1,KW,KTnocal,nH,A)+y_error_cal(X2,KW,KTcal,nH,B)+y_error_midcal(X3,KW,KTmidcal,nH,C)
[bootstrap_params(i,:),bootstrap_sum_sq_min(i,:)] = simulannealbnd(@(z) y_error(x_ind_nocal,x_ind_cal,x_ind_midcal,z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8)),[0.1,0.12,0.24,0.18,7,1,1,1],[0,0,0,0,0,0,0,0],[10,10,10,10,20,10,10,10],options);
end
%delete(gcp) %This optional step closes the parallel pool when the loop is done running.

%This calculates the confidence intervals (CIs) based on the bootstrapped data.
%It prints to the command window a table containing the confidence intervals
%for the best-fit values of the parameters floated during the fitting (param_CIs).
%The first row of param_CIs contains the measured values for each parameter.
%The second row contains the distance from the lower bound of the CI to the
%measured value, and the third row contains the distance from the upper bound
%of the CI to the measured value.
clear CI paramCIs param_CIs

rowNames = {'value','-','+'};
paramIDs = {'KW','KTnocal','KTcal','KTmidcal','nH','A','B','C'};
CI = prctile(bootstrap_params,[100*(1-CI_bound)/tail,100*(1-(1-CI_bound)/tail)]);
paramCIs = [best_fit_params;best_fit_params-CI(1,:);CI(2,:)-best_fit_params];
param_CIs = array2table(paramCIs,'RowNames',rowNames,'VariableNames',paramIDs)