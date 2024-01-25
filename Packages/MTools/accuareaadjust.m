function [n,p,wi,accu,area] = accuareaadjust(n,maparea)
% This is the function for adjust accuracy and area with standard error

% create confusion matrix
% n = confusionmat(reference,classification);

% number of classes
nc = length(maparea);

% weigts of each mapped class
wi = maparea/sum(maparea);

% sum of each rows
ni = sum(n,2);

% sume of each columns
nj = sum(n,1);

% adjusted n based on mapped area 
p = n;
for i = 1:nc
    for j = 1:nc
        p(i,j) = wi(i)*n(i,j)/ni(i);
    end
end

% calcuate adjusted overall accuracy, user's, and producer's accuracies
overall = 0;
U = zeros(1,nc);
P = zeros(1,nc);
pi = sum(p,2);
pj = sum(p,1);
for i = 1:nc
    overall = overall + p(i,i);
    U(i) = p(i,i)/pi(i);
    P(i) = p(i,i)/pj(i);
end
accu.overall = overall;

% estimated variance of overall accuracy
VO = 0;
for i = 1:nc
    VO = VO + wi(i)^2*U(i)*(1-U(i))/(ni(i)-1);
end
accu.overall_stderr = sqrt(VO);

% estimaed variance of user's accuracy of map class i
accu.user = U;

VU = zeros(1,nc);
for i = 1:nc
    VU(i) = U(i)*(1-U(i))/(ni(i)-1);
end
accu.user_stderr = sqrt(VU);

% estimated total number of pixels of reference class j
Nj = zeros(1,nc);

for j = 1:nc
    for i = 1:nc
        Nj(j) = Nj(j) + maparea(i)*n(i,j)/ni(i);
    end
end

% estimated variance of producer's accuracy of reference class j
accu.producer = P;

VP = zeros(1,nc);
for j = 1:nc
    VP(j) = (1/Nj(j)^2)*maparea(j)^2*(1-P(j))^2*U(j)*(1-U(j))/(nj(j)-1);
    for i = 1:nc
        if i ~= j
            VP(j) = VP(j) + (1/Nj(j)^2)*P(j)^2*maparea(i)^2*(n(i,j)/ni(i))*(1-n(i,j)/ni(i))/(ni(i)-1);
        end
    end
end
accu.producer_stderr = sqrt(VP);

% area adjusment
area.adjust_area = sum(maparea)*sum(p,1);

% estimated variance of area
VA = zeros(1,nc);

for k = 1:nc
    for i = 1:nc
        VA(k) = VA(k) + wi(i)^2*(n(i,k)/ni(i))*(1-n(i,k)/ni(i))/(ni(i)-1);
    end
end

area.adjust_area_stderr = sum(maparea)*sqrt(VA);

% end of the function
end

    
    
