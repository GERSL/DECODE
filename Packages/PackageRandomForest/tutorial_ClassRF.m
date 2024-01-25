% A simple tutorial file to interface with RF
% Options copied from http://cran.r-project.org/web/packages/randomForest/randomForest.pdf
%run plethora of tests
clear extra_options
clc
close all

%compile everything
if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
%   compile_windows
else
%   compile_linux
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
%randomly split into 250 examples for training and 50 for testing
randvector = randperm(N);

X_trn = X(randvector(1:250),:);
Y_trn = Y(randvector(1:250));
X_tst = X(randvector(251:end),:);
Y_tst = Y(randvector(251:end));

% example 1:  simply use with the defaults
    model = classRF_train(X_trn,Y_trn);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 1: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
 
% example 2:  set to 100 trees
    model = classRF_train(X_trn,Y_trn, 100);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 2: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 3:  set to 100 trees, mtry = 2
    model = classRF_train(X_trn,Y_trn, 100,2);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 3: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 4:  set to defaults trees and mtry by specifying values as 0
    model = classRF_train(X_trn,Y_trn, 0, 0);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 4: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 5: set sampling without replacement (default is with replacement)
    extra_options.replace = 0 ;
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 5: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 6: Using classwt (priors of classes)
    clear extra_options;
    extra_options.classwt = [1 1]; %for the [-1 +1] classses in twonorm
    % if you sort the labels in training and arrange in ascending order then
    % for twonorm you have -1 and +1 classes, with here assigning 1 to
    % both classes
    % As you have specified the classwt above, what happens that the priors are considered
    % also is considered the freq of the labels in the data. If you are
    % confused look into src/rfutils.cpp in normClassWt() function

    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 6: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 7: modify to make class(es) more IMPORTANT than the others
    %  extra_options.cutoff (Classification only) = A vector of length equal to
    %                       number of classes. The 'winning' class for an observation is the one with the maximum ratio of proportion
    %                       of votes to cutoff. Default is 1/k where k is the number of classes (i.e., majority
    %                       vote wins).    clear extra_options;
    extra_options.cutoff = [1/4 3/4]; %for the [-1 +1] classses in twonorm
    % if you sort the labels in training and arrange in ascending order then
    % for twonorm you have -1 and +1 classes, with here assigning 1/4 and
    % 3/4 respectively
    % thus the second class needs a lot less votes to win compared to the first class
    
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 7: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    fprintf('   y_trn is almost 50/50 but y_hat now has %f/%f split\n',length(find(Y_hat~=-1))/length(Y_tst),length(find(Y_hat~=1))/length(Y_tst));
    

%  extra_options.strata = (not yet stable in code) variable that is used for stratified
%                       sampling. I don't yet know how this works.

% example 8: sampsize example
    %  extra_options.sampsize =  Size(s) of sample to draw. For classification, 
    %                   if sampsize is a vector of the length the number of strata, then sampling is stratified by strata, 
    %                   and the elements of sampsize indicate the numbers to be drawn from the strata.
    clear extra_options
    extra_options.sampsize = size(X_trn,1)*2/3;
    
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 8: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    
% example 9: nodesize
    %  extra_options.nodesize = Minimum size of terminal nodes. Setting this number larger causes smaller trees
    %                   to be grown (and thus take less time). Note that the default values are different
    %                   for classification (1) and regression (5).
    clear extra_options
    extra_options.nodesize = 2;
    
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 9: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
        

% example 10: calculating importance
    clear extra_options
    extra_options.importance = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 10: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    
    %model will have 3 variables for importance importanceSD and localImp
    %importance = a matrix with nclass + 2 (for classification) or two (for regression) columns.
    %           For classification, the first nclass columns are the class-specific measures
    %           computed as mean decrease in accuracy. The nclass + 1st column is the
    %           mean decrease in accuracy over all classes. The last column is the mean decrease
    %           in Gini index. For Regression, the first column is the mean decrease in
    %           accuracy and the second the mean decrease in MSE. If importance=FALSE,
    %           the last measure is still returned as a vector.
    figure('Name','Importance Plots')
    subplot(2,1,1);
    bar(model.importance(:,end-1));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Accuracy');
    
    subplot(2,1,2);
    bar(model.importance(:,end));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Gini index');
    
    
    %importanceSD = The ?standard errors? of the permutation-based importance measure. For classification,
    %           a D by nclass + 1 matrix corresponding to the first nclass + 1
    %           columns of the importance matrix. For regression, a length p vector.
    model.importanceSD

% example 11: calculating local importance
    %  extra_options.localImp = Should casewise importance measure be computed? (Setting this to TRUE will
    %                   override importance.)
    %localImp  = a D by N matrix containing the casewise importance measures, the [i,j] element
    %           of which is the importance of i-th variable on the j-th case. NULL if
    %          localImp=FALSE.
    clear extra_options
    extra_options.localImp = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 11: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

    model.localImp
    
% example 12: calculating proximity
    %  extra_options.proximity = Should proximity measure among the rows be calculated?
    clear extra_options
    extra_options.proximity = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 12: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

    model.proximity
    

% example 13: use only OOB for proximity
    %  extra_options.oob_prox = Should proximity be calculated only on 'out-of-bag' data?
    clear extra_optionsv
    extra_options.proximity = 1; %(0 = (Default) Don't, 1=calculate)
    extra_options.oob_prox = 0; %(Default = 1 if proximity is enabled,  Don't 0)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 13: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));


% example 14: to see what is going on behind the scenes    
%  extra_options.do_trace = If set to TRUE, give a more verbose output as randomForest is run. If set to
%                   some integer, then running output is printed for every
%                   do_trace trees.
    clear extra_options
    extra_options.do_trace = 1; %(Default = 0)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 14: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 14: to see what is going on behind the scenes    
%  extra_options.keep_inbag Should an n by ntree matrix be returned that keeps track of which samples are
%                   'in-bag' in which trees (but not how many times, if sampling with replacement)
%
    clear extra_options
    extra_options.keep_inbag = 1; %(Default = 0)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 15: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    
    model.inbag

% example 16: getting the OOB rate. model will have errtr whose first
% column is the OOB rate. and the second column is for the 1-st class and
% so on
    model = classRF_train(X_trn,Y_trn);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 16: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    
    figure('Name','OOB error rate');
    plot(model.errtr(:,1)); title('OOB error rate');  xlabel('iteration (# trees)'); ylabel('OOB error rate');
    

% example 17: getting prediction per tree, votes etc for test set (returns
% prediction_per_tree a Ntst x Ntree matrix)
    model = classRF_train(X_trn,Y_trn);
    
    test_options.predict_all = 1;
    [Y_hat, votes, prediction_per_tree] = classRF_predict(X_tst,model,test_options);
    fprintf('\nexample 17: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    
% example 18: proximity of test cases
% proximity values in proximity_ts
    model = classRF_train(X_trn,Y_trn);
    
    test_options.predict_all = 1;
    test_options.proximity = 1;
    [Y_hat, votes, prediction_per_tree, proximity_ts] = classRF_predict(X_tst,model,test_options);
    fprintf('\nexample 18: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

% example 19: to verbosely print the forest creation status, i.e. print after creation of each tree
    clear extra_options
    extra_options.print_verbose_tree_progression = 1; %(Default = 0)
   
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('\nexample 19: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    
% example 20: get the nodes for test examples
%   terminal node indicators (ntest by ntree matrix)
    clear extra_options
    extra_options.nodes = 1;
    [Y_hat, tmp, tmp, tmp, nodes] = classRF_predict(X_tst,model,extra_options);
    fprintf('\nexample 20: node information now in node variable\n');

% example 21: using categories
    fprintf('\nexample 21: example of using categories\n');
    load data/twonorm
    X = inputs';
    Y = outputs;

    default_scale = 8; %X values will be scaled between 0-(default_scale-1)
    for i=1:size(X,2)
        min_0_X = double(X(:,i) - min(X(:,i))); % minimum value will be now 0.
        X_i_between_0_1 = min_0_X/max(min_0_X); %now the range will be 0-1
        X_i_between_0_scale_2_to_n = round(X_i_between_0_1 * (default_scale-1));
        X(:,i) = X_i_between_0_scale_2_to_n;
    end
    % X is fully categorized and ranges from values 0:7

    % construct the category information for X
    % .categories is a 1 x D vector saying what features are categorical and what are not categorical
    % true = categorical, false = numeric

    extra_options.categorical_feature = ones(1,size(X,2)); % all features are categorical

    % to choose features 1:5 as non-categorcal, use below code
    % extra_options.categorical_feature(1:5) = false; 

    % without categories
    model=classRF_train(X,Y, 500, 3);
    y_hat = classRF_predict(X, model);
    err_without_categories = length(find(y_hat~=Y))/length(Y)
    avg_size_of_trees = mean(model.ndbigtree)

    % with categories
    model=classRF_train(X,Y, 500, 3, extra_options);
    y_hat = classRF_predict(X, model);
    err_with_categories = length(find(y_hat~=Y))/length(Y)
    avg_size_of_trees = mean(model.ndbigtree)

    fprintf('\terror with/without categories should be similar, categorical trees might be smaller than numerical trees\n')

% example 22:  Stratified Sampling, change the sampling rate for individual classes by changing sampsize
    clear extra_options
    extra_options.sampsize = [10 100]; % -1 is sampled 10 times, 1 is samples 100 times (sampsize is set with values given for classes [sorted in ascending order])
    fprintf('\nexample 22: ');
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    fprintf('Because we changed the sampling times for the classes, \n   +1 is sampled %d times,\n   -1 is sampled %d times \n(if the sampling was the same then sampling would have been equal)\n', round(mean(model.oob_times(find(Y_trn==1)))), round(mean(model.oob_times(find(Y_trn==-1)))) )
    fprintf('seems like sampsize has an inverse role when strata is NOT involved, the larger the sampsize the lower the probability of something being sampled\n')
    
% example 23:  Stratified Sampling, instead of changing the sampsize variable, we use a different variable as the strata
    clear extra_options
    extra_options.strata = [ones(125,1); ones(125,1)*2]; %first 125 examples are in strata 1 the second 125 examples are in strata 2
    extra_options.sampsize = [200 10]; % strata 1 is sampled 10 times, strata 2 is sampled at 100 times
    fprintf('\nexample 23: ');
    model = classRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = classRF_predict(X_tst,model);
    fprintf('error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    fprintf('Because we changed the sampling times for the strata, \n   first 125 examples are sampled %d times,\n   second 125 examples are sampled %d times \n(if the sampling was the same then sampling would have been equal)\n', round(mean(model.oob_times(1:125))), round(mean(model.oob_times(126:end))) )
    fprintf('seems like sampsize has an inverse role when strata is involved, the larger the sampsize the lesser probability of something being sampled\n')
        