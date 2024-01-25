%this predicts for all the trees in the RF
function [Y_new, Y_new_per_tree] = classRF_predict_matlab(X,classRF_model)
    %Y_new is the winner vote from entire RF
    %Y_new_per_tree is the winner vote for individual tree in RF
    

    %treemap is a variable that holds the treemap. if you need more infor about it i might have to look deeper
    %but i dont remember it from the top of my head.
    classRF_model.treemap2 = zeros(classRF_model.nrnodes*2,classRF_model.ntree);
    for treenum=1:classRF_model.ntree
        classRF_model.treemap2(:,treenum)=[classRF_model.treemap(:,(treenum-1)*2+1); classRF_model.treemap(:,(treenum-1)*2+2)];
    end
    
    %countts number of examples
    countts = zeros(size(X,1),classRF_model.nclass);
    %jts prediction for examples x trees
    jts = zeros(size(X,1),classRF_model.ntree);
    
    %now iterate through each of the tree and get the prediction for the tree
    for j=1:classRF_model.ntree
        [tmp_jts,nodex] = predictTree(X,classRF_model,j);
        jts(:,j) = tmp_jts;
        for i=1:size(X,1)
            countts(i,tmp_jts(i)) = countts(i,tmp_jts(i)) +1;
        end
    end
    
    
    %below is strictly needed if the classweights are not equal to 1, then crit comes into the picture
    %ties are also broken here
    tied_index = zeros(size(X,1),1);
    jet = zeros(size(X,1),1);
    for n=1:size(X,1)
        cmax = 0;
        ntie = 1;
        for j=1:classRF_model.nclass
            crit = (countts(n,j)/classRF_model.ntree) / classRF_model.cutoff(j);
            if(crit > cmax)
                jet(n) = j;
                cmax = crit;
            end
            if(crit == cmax)
                ntie = ntie+1;
                if (rand > 1.0/ntie)
                    jet(n)=j;
                    tied_index(n)=1;
                end
            end
        end
    end
    Y_hat = jet;
    Y_hat_per_tree = jts;
    
    Y_new = double(Y_hat);
    Y_new_per_tree = double(Y_hat_per_tree);
    new_labels = classRF_model.new_labels;
    orig_labels = classRF_model.orig_labels;
    
    for i=1:length(orig_labels)
        Y_new(find(Y_hat==new_labels(i)))=Inf;
        Y_new(isinf(Y_new))=orig_labels(i);
        Y_new_per_tree(find(Y_hat_per_tree==new_labels(i)))=Inf;
        Y_new_per_tree(isinf(Y_new_per_tree))=orig_labels(i);
    end
end




%this predicts per tree
% pass it the examples to predict, the RF model and the tree number in the RF ensemble
function [jts,nodex] = predictTree(X,classRF_model,treenum)
    %predict this for tree 'treenum'
    
    NODE_TERMINAL=-1;
    NODE_TOSPLIT =-2;
    NODE_INTERIOR=-3;
    
%    iterate through each example
    for i=1:size(X,1)
        k=1;
        %start with k=1 and then go on till we reach a terminal node
        %nodestatus is i think numnodes x numberof trees
        while (classRF_model.nodestatus(k,treenum) ~= NODE_TERMINAL)
    	    %m is the variable that was used to split
            m = classRF_model.bestvar(k,treenum);
            
            %now that we know m we can find if X(current_example,m) <= the split value
            %then either go right or left in the tree.
                        
            % if the X's value is less then the split go left else go right
            if X(i,m)  <= classRF_model.xbestsplit(k,treenum)
                kold = k;
                %k = treemap((k-1)*2+1);
                k=classRF_model.treemap2((k-1)*2+1,treenum);
            else
                kold = k;
                %k = treemap((k-1)*2+2);
                k=classRF_model.treemap2((k-1)*2+2,treenum);
            end
        end
        jts(i) = classRF_model.nodeclass(k,treenum);
        nodex(i) = k+1;
    end
end