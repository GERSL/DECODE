function predictions = senPredict(datain,Sen,mode,isPlot)
    if ~exist('mode','var')
        mode = 'first';
    end
    if ~exist('isPlot','var')
        isPlot = 0;
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
%     if isPlot
% %     senplot = [datain(:,1) predictions(:,2)];
% %     plot_Sen = plot(datain(:,1),datain(:,2),'o','Color', 'r','DisplayName', 'slope L8');
% %     hold on;
% %     plot_Sen = plot(datain(:,1),predictions(:,2),'-','Color', 'g','LineWidth',2,'DisplayName', 'Sen slope L8');
% %     hold on;
%     end
end