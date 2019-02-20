%This MATLAB code is associated with the following manuscript: Barrick,
%S.K., S.R. Clippinger, L. Greenberg, M.J. Greenberg. 2019. Computational
%tool to study perturbations in muscle regulation and its application to 
%heart disease.

%This script performs the fitting of the data to the McKillop and Geeves
%model without global fitting. The code extends the original method to
%include error estimation of the fit parameters.

%%

%Here, we have data collected at two calcium concentrations: 2 mM EGTA (nocal)
%and pCa 3 (cal).

%This assigns the measured values to variables for the rest of the
%program. The user can change the names so their data are assigned to the proper
%variables.
x1=s1_nocal; %This is the myosin concentration for the data collected at low calcium (2 mM EGTA)
y1=fl_nocal; %This is the fractional change in fluorescence for the data collected at low calcium (2 mM EGTA)
x2=s1_cal; %This is the myosin concentration for the data collected at high calcium (pCa 3)
y2=fl_cal; %This is the fractional change in fluorescence for the data collected at high calcium (pCa 3)

%This assigns the variables appropriate values.
KBnocal=0.290; %User should change this to their value for KB from the stopped-flow experiments.
KS=18; %This value is fixed as described in McKillop and Geeves (1993).
KBcal=20; %This value is fixed as described in McKillop and Geeves (1993).

%These are the parameters for the confidence interval calculations.
N_bs=1000; %Number of bootstrapping simulations.  1000 is usually sufficient.
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

%The following code performs the fitting of the two curves via simulated
%annealing and least-squares minimization. Note that this is not global fitting,
%as each curve is fit individually.  It returns the set of best-fit parameters
%(best_fit_params) as well as the value of the sum of squares at the
%minimum (sum_sq_min). The order of parameters in best_fit_params is the
%same as in y_error. Fixed parameters are the values of KBnocal, KBcal, and
%KS. Fitted parameters are nH (the Hill coefficient), KW, KT at low
%(KTnocal) and high (KTcal) calcium (KTmidcal), and the
%normalization amplitudes (A,B).

%This defines y_error for the nocal and cal conditions as the sum
%of squared differences between Equation 1 (see user guide) and pooled
%fractional change in fluorescence data for that condition (y1, y2).
y_error_nocal = @(X1,KW,KTnocal,nH,A) sum(((A.*KW.*x1.*(1 + KW.*x1.*(1 + KS)).^(nH - 1).*(KTnocal.*(1 + KS).^nH + 1))./((KTnocal.*(1 + KW.*x1.*(1 + KS)).^nH + (1 + KW.*x1).^nH + 1./KBnocal).*(1 + KS).^(nH - 1)) - y1).^2);
y_error_cal = @(X2,KW,KTcal,nH,B) sum(((B.*KW.*x2.*(1 + KW.*x2.*(1 + KS)).^(nH - 1).*(KTcal.*(1 + KS).^nH + 1))./((KTcal.*(1 + KW.*x2.*(1 + KS)).^nH + (1 + KW.*x2).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y2).^2);

%This line performs the fitting by minimizing y_error_nocal.
[best_fit_params_nocal,sum_sq_min_nocal] = simulannealbnd(@(z) y_error_nocal(x1,z(1),z(2),z(3),z(4)),[0.03,0.12,7,1],[0,0,0,0],[10,10,20,10],options);

%This line performs the fitting by minimizing y_error_cal.
[best_fit_params_cal,sum_sq_min_cal] = simulannealbnd(@(z) y_error_cal(x2,z(1),z(2),z(3),z(4)),[0.03,0.24,7,1],[0,0,0,0],[10,10,20,10],options);

%This prints to the command window a table containing the best-fit values
%of the parameters floated during each fit (best_fit_values_nocal and
%best_fit_values_cal), as well as the value of the sum of squares at the
%minimum for each fit (sum_sq_min_nocal and sum_sq_min_cal).
paramIDs_nocal = {'KW','KTnocal','nH','A'};
paramIDs_cal = {'KW','KTcal','nH','B'};
best_fit_values_nocal = array2table(best_fit_params_nocal,'VariableNames',paramIDs_nocal)
sum_sq_min_nocal
best_fit_values_cal = array2table(best_fit_params_cal,'VariableNames',paramIDs_cal)
sum_sq_min_cal

%The following code performs the bootstrapping simulations for error estimation.
%It outputs a matrix bootstrap_params that contains the fit values
%for each of the simulated data sets. The columns give the values of the
%parameters and each row contains values obtained from a separate simulation.

%This randomly selects indices of the original, pooled data set to
%generate the simulated (i.e., bootstrapped) data sets.
bs_ind_nocal = ceil(length(x1)*rand(length(x1),N_bs));
bs_ind_cal = ceil(length(x2)*rand(length(x2),N_bs));

%This generates fits to the N_bs simulated data sets. Each row of
%bootstrap_fit_params contains the best-fit parameters for one round of
%resampling. The order of parameters in the columns of bootstrap_params is
%the same as in y_error.
parfor i=1:N_bs
    x_ind_nocal = x1(bs_ind_nocal(:,i));
    y_ind_nocal = y1(bs_ind_nocal(:,i));
    x_ind_cal = x2(bs_ind_cal(:,i));
    y_ind_cal = y2(bs_ind_cal(:,i));

y_error_nocal = @(X1,KW,KTnocal,nH,A) sum(((A.*KW.*X1.*(1 + KW.*X1.*(1 + KS)).^(nH - 1).*(KTnocal.*(1 + KS).^nH + 1))./((KTnocal.*(1 + KW.*X1.*(1 + KS)).^nH + (1 + KW.*X1).^nH + 1./KBnocal).*(1 + KS).^(nH - 1)) - y_ind_nocal).^2);
[bootstrap_params_nocal(i,:),bootstrap_sum_sq_min(i+1,:)] = simulannealbnd(@(z) y_error_nocal(x_ind_nocal,z(1),z(2),z(3),z(4)),[0.03,0.12,7,1],[0,0,0,0],[10,10,20,10],options)

y_error_cal = @(X2,KW,KTcal,nH,B) sum(((B.*KW.*X2.*(1 + KW.*X2.*(1 + KS)).^(nH - 1).*(KTcal.*(1 + KS).^nH + 1))./((KTcal.*(1 + KW.*X2.*(1 + KS)).^nH + (1 + KW.*X2).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y_ind_cal).^2);
[bootstrap_params_cal(i,:),bootstrap_sum_sq_min(i+1,:)] = simulannealbnd(@(z) y_error_cal(x_ind_cal,z(1),z(2),z(3),z(4)),[0.03,0.24,7,1],[0,0,0,0],[10,10,20,10],options)

end
%delete(gcp) %This optional step closes the parallel pool when the loop is done running.

%This calculates the confidence intervals based on the bootstrapped data.
%It returns param_CI_nocal and param_CI_cal. The first row contains the
%measured values for each parameter. The second row contains the distance
%from the lower bound of the CI to the measured value, and the third row
%contains the distance from the upper bound of the CI to the measured value.
clear CI_nocal paramCIs_nocal CI_cal paramCIs_cal param_CI_nocal param_CI_cal

rowNames = {'value','-','+'};
paramIDs_nocal = {'KW','KTnocal','nH','A'};
paramIDs_cal = {'KW','KTcal','nH','B'};
CI_nocal = prctile(bootstrap_params_nocal,[100*(1-CI_bound)/tail,100*(1-(1-CI_bound)/tail)]);
paramCIs_nocal = [best_fit_params_nocal;best_fit_params_nocal-CI_nocal(1,:);CI_nocal(2,:)-best_fit_params_nocal];
CI_cal = prctile(bootstrap_params_cal,[100*(1-CI_bound)/tail,100*(1-(1-CI_bound)/tail)]);
paramCIs_cal = [best_fit_params_cal;best_fit_params_cal-CI_cal(1,:);CI_cal(2,:)-best_fit_params_cal];
param_CI_nocal = array2table(paramCIs_nocal,'RowNames',rowNames,'VariableNames',paramIDs_nocal)
param_CI_cal = array2table(paramCIs_cal,'RowNames',rowNames,'VariableNames',paramIDs_cal)