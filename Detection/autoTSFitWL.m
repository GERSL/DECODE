function [fit_cft,rmse,v_dif]=autoTSFitWL(x,y,df,wl)
% Revisions: 
% v1.0 Using lasso for timeseries modeling (01/27/2013)
% Auto Trends and Seasonal Fit between breaks
% INPUTS:
% x - Julian day [1; 2; 3];
% y - predicted reflectances [0.1; 0.2; 0.3];
% df - degree of freedom (num_c)
% wl - tide water level
% OUTPUTS:
% fit_cft - fitted coefficients;

n=length(x); % number of clear pixels
% num_yrs = 365.25; % number of days per year
w=2*pi/365.25; % num_yrs; % anual cycle
% fit coefs
fit_cft = zeros(9,1);

% global lamda_lasso

%% LASSO Fit
% build X
X = zeros(n,df-1);
X(:,1) = x;

if df >= 5
    X(:,2)=cos(w*x);
    X(:,3)=sin(w*x);
end

if df >= 7
    X(:,4)=cos(2*w*x);
    X(:,5)=sin(2*w*x);
end

if df >= 9
    X(:,6)=cos(3*w*x);
    X(:,7)=sin(3*w*x);
end
X(:,end)= wl;
% lasso fit with lambda = 20
fit = glmnet_fast(X,y,glmnetSetL(20));  

% curr_cft=[fit.a0;fit.beta];
fit_cft(1:df) = [fit.a0;fit.beta]; % curr_cft;


yhat=autoTSPredWL(x,fit_cft,wl);

% rmse=median(abs(y-yhat));
v_dif = y-yhat;

% rmse=norm(v_dif)/sqrt(n-df);
rmse=norm(v_dif)/sqrt(n-df);
% f(x) = a0 + b0*x + a1*cos(x*w) + b1*sin(x*w) (df = 4)
%% Double check the problem caused by water level
if rmse > 1000 || sum(v_dif>0)/sum(v_dif<0)>2 || sum(v_dif<0)/sum(v_dif>0)>2
    X(:,end)= 0; 
    % lasso fit with lambda = 20
    fit = glmnet_fast(X,y,glmnetSetL(20));  
    
    % curr_cft=[fit.a0;fit.beta];
    fit_cft(1:df) = [fit.a0;fit.beta]; % curr_cft;
    
    
    yhat=autoTSPredWL(x,fit_cft,wl);
    
    % rmse=median(abs(y-yhat));
    v_dif = y-yhat;
    
    % rmse=norm(v_dif)/sqrt(n-df);
    rmse=norm(v_dif)/sqrt(n-df);
end

end
