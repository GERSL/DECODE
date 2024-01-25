% A simple tutorial file to interface with RF
% Options copied from http://cran.r-project.org/web/packages/randomForest/randomForest.pdf

%run plethora of tests
clc
close all
clear options
clear extra_options

%compile everything
if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
%    compile_windows
else
%    compile_linux
end

total_train_time=0;
total_test_time=0;

%load the twonorm dataset
load data/twonorm

%modify so that training data is NxD and labels are Nx1, where N=#of
%examples, D=# of features

X = inputs';
Y = outputs;

[N D] =size(X);

X_new = X;
Y_new = []; %lose labels for Y before sending it to the program

extra_options.do_trace = 1; %(Default = 0)
model = classRF_train(X_new,Y_new,100,0,extra_options);

%copied this from MDSplot.R from the randomForest's R source
%basically use cmdscale command for multidimensional scaling and
%then plot points in 2 dimensions.
rf_mds = cmdscale(1-model.proximity);
plot(rf_mds(Y==1,1),rf_mds(Y==1,2),'*k');hold on
plot(rf_mds(Y==-1,1),rf_mds(Y==-1,2),'*r')