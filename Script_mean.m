%This MATLAB code is associated with the following manuscript: Barrick,
%S.K., S.R. Clippinger, L. Greenberg, M.J. Greenberg. 2019. Computational
%tool to study perturbations in muscle regulation and its application to 
%heart disease.

%This is a script to calculate the mean value for a population of objects
%and the associated uncertainties.

%%

%Here, we input a range of measured values into x1.
x1=data;

%These are the parameters for the confidence interval calculations.
N_bs=1000; %Number of bootstrapping simulations. 1000 is usually sufficient.
CI_bound=0.95; %This is the percent value for the confidence intervals.
tail=2; %This defines whether the p-value is obtained for a 1- or 2-tailed
%test. Use 2 by default.

%This removes any non-numerical values in the imported data set.
q=find(isnan(x1));
x1(q)=[];

%This calculates the mean of the real data.
best_fit_param=mean(x1);

%The following code performs the bootstrapping simulations for error estimation.
%It outputs a matrix, bootstrap_param, that contains the mean value
%for each of the simulated data sets. Each row contains the value obtained
%from a separate simulation.

%This randomly selects indices of the original data set used to
%generate the simulated (i.e., bootstrapped) data sets.
bs_ind=ceil(length(x1)*rand(length(x1),N_bs));

%This generates fits to the N_bs simulated data sets. Each row of
%bootstrap_param contains the mean for one round of resampling.
bootstrap_param = zeros(length(data)-1, 1);

parfor i=1:N_bs
    x_ind=x1(bs_ind(:,i));
 
bootstrap_param(i,:)=mean(x_ind)
end

%This calculates the confidence intervals based on the bootstrapped data.
%It returns param_CI. The first row contains the measured mean value.  The
%second row contains the distance from the lower bound of the CI to the
%measured value, and the third row contains the distance from the upper
%bound of the CI to the measured value.
clear CI paramCI param_CI

rowNames = {'mean','-','+'};
CI = prctile(bootstrap_param,[100*(1-CI_bound)/tail,100*(1-(1-CI_bound)/tail)]);
paramCI = [best_fit_param;best_fit_param-CI(1);CI(2)-best_fit_param];
param_CI = array2table(paramCI,'RowNames',rowNames)