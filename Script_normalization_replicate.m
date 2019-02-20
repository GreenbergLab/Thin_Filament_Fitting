%This MATLAB code is associated with the following manuscript: Barrick,
%S.K., S.R. Clippinger, L. Greenberg, M.J. Greenberg. 2019. Computational
%tool to study perturbations in muscle regulation and its application to 
%heart disease.

%This is a script for normalization of a single technical replicate.

%%

%Normalization of technical replicates

%This script performs fitting of fluorescence titration data to obtain the
%normalization constants A, B, and C for data collected at 2 mM EGTA (nocal),
%pCa 3 (cal), and pCa 6.25 (midcal), respectively. For each titration, the
%myosin concentration and fractional change in fluorescence at that
%concentration are recorded.

%This assigns the measured values to variables for the rest of the program.
%The user can change the names so their data are assigned to the proper
%variables.
x1=s1_nocal; %This is the myosin concentration for the data collected at low calcium (2 mM EGTA)
y1=fl_nocal; %This is the fractional change in fluorescence for the data collected at low calcium (2 mM EGTA)
x2=s1_cal; %This is the myosin concentration for the data collected at high calcium (pCa 3)
y2=fl_cal; %This is the fractional change in fluorescence for the data collected at high calcium (pCa 3)
x3=s1_midcal; %This is the myosin concentration for the data collected at intermediate calcium (pCa 6.25)
y3=fl_midcal; %This is the fractional change in fluorescence for the data collected at intermediate calcium (pCa 6.25)

%This assigns the variables appropriate values.
KBnocal=0.290; %User should change this to their value for KB from the stopped-flow experiments
KS=18; %This value is fixed as described in McKillop and Geeves (1993).
KBcal=20; %This value is fixed as described in McKillop and Geeves (1993).

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
clear y_error* best_fit_params sum_sq_min *_norm

%The following code performs the fit via simulated annealing and minimization of
%y_error. Here, we are interested in obtaining the normalization amplitudes
%(A, B, C) to enable us to normalize each technical replicate before combining
%the data into a pooled data set for global fitting and error estimation.

%This defines y_error for the nocal, cal, and midcal conditions as the sum
%of squared differences between Equation 1 (see user guide) and the fractional
%change in fluorescence data for that condition (y1, y2, y3).
y_error_nocal = @(X1,KW,KTnocal,nH,A) sum(((A.*KW.*x1.*(1 + KW.*x1.*(1 + KS)).^(nH - 1).*(KTnocal.*(1 + KS).^nH + 1))./((KTnocal.*(1 + KW.*x1.*(1 + KS)).^nH + (1 + KW.*x1).^nH + 1./KBnocal).*(1 + KS).^(nH - 1)) - y1).^2);
y_error_cal = @(X2,KW,KTcal,nH,B) sum(((B.*KW.*x2.*(1 + KW.*x2.*(1 + KS)).^(nH - 1).*(KTcal.*(1 + KS).^nH + 1))./((KTcal.*(1 + KW.*x2.*(1 + KS)).^nH + (1 + KW.*x2).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y2).^2);
y_error_midcal = @(X3,KW,KTmidcal,nH,C) sum(((C.*KW.*x3.*(1 + KW.*x3.*(1 + KS)).^(nH - 1).*(KTmidcal.*(1 + KS).^nH + 1))./((KTmidcal.*(1 + KW.*x3.*(1 + KS)).^nH + (1 + KW.*x3).^nH + 1./KBcal).*(1 + KS).^(nH - 1)) - y3).^2);

%The error from all three curves is summed to create y_error, which is later
%minimized in the global fit.
y_error = @(X1,X2,X3,KW,KTnocal,KTcal,KTmidcal,nH,A,B,C) y_error_nocal(X1,KW,KTnocal,nH,A)+y_error_cal(X2,KW,KTcal,nH,B)+y_error_midcal(X3,KW,KTmidcal,nH,C);

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

%This performs the global fit by minimizing y_error.
[best_fit_params,sum_sq_min] = simulannealbnd(@(z) y_error(x1,x2,x3,z(1),z(2),z(3),z(4),z(5),z(6),z(7),z(8)),[0.03,0.12,0.24,0.18,7,1,1,1],[0,0,0,0,0,0,0,0],[10,10,10,10,20,10,10,10],options);

%This returns a table containing the normalization amplitudes for nocal (A),
%cal (B), and midcal (C), as well as the value of the sum of squares at the
%minimum (sum_sq_min).
normIDs = {'A','B','C'};
norm_constants = array2table(best_fit_params(6:8),'VariableNames',normIDs)
sum_sq_min

%The normalized fluorescence values and corresponding myosin concentrations
%are then output to the workspace. The normalized values for each technical
%replicate should be pooled to form a pooled data set, used for the global fitting.
nocal_norm = [s1_nocal y1./best_fit_params(6)]; %Normalized data for nocal
cal_norm = [s1_cal y2./best_fit_params(7)]; %Normalized data for cal
midcal_norm = [s1_midcal y3./best_fit_params(8)]; %Normalized data for midcal