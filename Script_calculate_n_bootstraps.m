%This MATLAB code is associated with the following manuscript: Barrick,
%S.K., S.R. Clippinger, L. Greenberg, M.J. Greenberg. 2019. Computational
%tool to study perturbations in muscle regulation and its application to 
%heart disease.

%This script tests whether a sufficient number of bootstraps have been
%performed.

%%

%The input is a matrix of values obtained from bootstrapping (bootstrap_params). A large number of bootstraps should be used. 

input=bootstrap_params;

ind=1; %This is the index of the parameter to be compared. For example, 
%if best_fit_params is ordered as (KW,KTnocal,KTcal,KTmidcal,nH,A,B,C), 
%then ind=1 selects KW, ind=2 selects KTnocal, etc.

%These are the parameters for the confidence interval calculations.
CI_bound=0.95; %This is the percent value for the confidence intervals.
tail=2; %This defines whether the p-value is obtained for a 1- or 2-tailed
%test. Use 2 by default.

%This clears the variables used to make it easier to run the script consecutively
%for multiple indices.
clear data_subset data_mean CI_data_subset

data=input(:,ind); %This selects the desired index of the output matrix

%This loop calculates the mean value obtained if only i bootstraps were
%run, along with the confidence intervals for this number of bootstraps.
data_mean = zeros(length(data)-1, 0);
CI_data_subset = zeros(2, length(data)-1);
for i=2:length(data)
    clear data_subset
    data_subset = data(1:i);
    data_mean(i-1) = mean(data_subset);
    CI_data_subset(:,i-1) = prctile(data_subset,[100*(1-CI_bound)/tail,100*(1-(1-CI_bound)/tail)]);
end

x_axis=2:length(data); %This generates the x-axis.

%This generates a plot. The x-axis is the number of bootstraps. The
%black line is the mean parameter value from that number of bootstraps.
%The red and blue lines show the upper and lower bounds, respectively, of the
%95% confidence intervals. If a sufficient number of bootstraps have been
%performed, the values of the mean and bounds of the confidence intervals
%should change only minimally at large values of x (number of bootstraps).

figure(1)
hold off
plot(x_axis,data_mean,'k.')
hold on
plot(CI_data_subset(1,:),'b')
plot(CI_data_subset(2,:),'r')