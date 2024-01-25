function [AC,APR,AAPC] = calculateAAPC(datain,Sen,mode)
    
    if ~exist('mode','var')
        mode = 'first';
    end
    if strcmp(mode,'median')
        % Plot sen's slope
        vv = median(datain(:,2));        
        middata = datain(round(length(datain)/2),1);
        predictions = vv + Sen*(datain(:,1)-middata);
        predictions = [datain(:,1) predictions];
    elseif strcmp(mode,'first')
        % Use the initial observation and sen's slope to predict the
        % further observations
        initialValue = datain(1,2);
        predictions = initialValue + Sen*(datain(:,1)-datain(1,1));
        predictions = [datain(:,1) predictions];
    end
    % 
    if Sen~=0
        % Annual change 
        AC = (predictions(end,2)-predictions(1,2))/length(datain);
        % annual percentage growth rate
        APR = (predictions(end,2)-predictions(1,2))/(predictions(1,2)+eps)/length(datain)*100;
        % Calculate average annual (compound) growth rate
        if predictions(end,2)*predictions(1,2)<0
            AAPC = NaN;
        else
            AAPC = power(predictions(end,2)/(predictions(1,2)+eps),1/length(predictions))-1;
            AAPC = AAPC * 100;
        end
    else
        AC = 0;
        APR = 0;
        AAPC = 0;
    end

% 
%     senplot = [datain(:,1) predictions(:,2)];
%     plot_Sen = plot(datain(:,1),datain(:,2),'o','Color', 'r','DisplayName', 'slope L8');
%     hold on;
%     plot_Sen = plot(datain(:,1),predictions(:,2),'-','Color', 'g','LineWidth',2,'DisplayName', 'Sen slope L8');
%     hold on;
end