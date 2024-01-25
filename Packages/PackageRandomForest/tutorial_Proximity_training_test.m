%tutorial_Proximity_training_test

%how to get proximity values. 
% 2 cases
% 1. Test Data only. i.e. proximity between test examples and training data only. You do get
%    proximity values for training examples but there is no proximity values
%    inbetween the test and the training
% 2. Test and Training data both. i.e. proximity between test and training.
%    Proximity values on test, training and inbetween them.


if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
   compile_windows
else
   compile_linux
end
%load the twonorm dataset 
load data/twonorm
 
%modify so that training data is NxD and labels are Nx1, where N=#of
%examples, D=# of features

X = inputs';
Y = outputs;

[N D] =size(X);
%randomly split into 250 examples for training and 50 for testing
randvector = randperm(N);

X_trn = X(randvector(1:10),:); %10 examples in training
Y_trn = Y(randvector(1:10));
X_tst = X(randvector(295:end),:); %5 examples in test
Y_tst = Y(randvector(295:end));

%%% use training and test as the same
%  X_tst = X_trn;
%  Y_tst = Y_trn;

% example 1: proximity of test cases
% First training and then finding the proximity on test separately
% proximity values in proximity_ts

% 1. Test Data only. i.e. proximity between test examples only. You do get
%    proximity values for training examples but there is no proximity values
%    between the test and the training


    model1 = classRF_train(X_trn,Y_trn,2000, 0);
    
    test_options.predict_all = 1;
    test_options.proximity = 1;
    [Y_hat, votes, prediction_per_tree, proximity_ts] = classRF_predict(X_tst,model1,test_options);
    fprintf('\nexample 1: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));
    

% example 2: proximity of test cases
% Finding the proximity on training and test simultaneously


% 2. Test and Training data both. i.e. proximity between test and training.
%    Proximity values on test, training and inbetween them.
% just pass X_tst and Y_tst in addition to the training data with
% classRF_train

    clear extra_options
    extra_options.proximity = 1; %(0 = (Default) Don't, 1=calculate)
    model2 = classRF_train(X_trn,Y_trn, 2000, 0, extra_options,X_tst,Y_tst);
    
    Y_hat = classRF_predict(X_tst,model2);
    fprintf('\nexample 2: error rate %f\n',   length(find(Y_hat~=Y_tst))/length(Y_tst));

    
%this matrix has the proximity on the test, first Ntst x Ntst values are between the test the rest are for training
%i.e the size of this matrix is Ntst x (Ntst + Ntrn). Ntst=size of test, Ntrn = size of training
model2.proximity_tst

%this matrix has the proximity values between only the training examples.
%size of this matrix is Ntst x Ntst
proximity_ts
