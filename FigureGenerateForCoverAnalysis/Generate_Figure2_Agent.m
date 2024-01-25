% This function generate Figure 2 with six subfigures; Including (1) agent;
% (2) annual cover change;(3)-(6) different covers
% Source code:

clc
clear
close all
restoredefaultpath

fontName = 'Arial';
fontSize = 14;
coverColors = {'#009e73','#f0e442','#cc79a7','#e69f00','#0072b2','#000000'};
changeColors = {'#009e73','#d55e00','#0072b2','#f0e442','#000000'};
classNames = {'Tidalmarsh','Mangrove','Dieback','Tidalflats','Openwater','Others'};

folderpath_mfile = fileparts(fileparts(mfilename('fullpath')));
addpath(folderpath_mfile);
addpath(fullfile(folderpath_mfile,'Packages/MTools')); % "Good practice" to calculate the uncertainty


pathFigure = fullfile(fileparts(mfilename('fullpath')),'Figures');
if ~exist(pathFigure,'dir')
    mkdir(pathFigure)
end

pathDECODE = fileparts(fileparts(globalsets.PathDECODE));

folderCoverAnalyses = 'CoverAnalyses';
pathCoverAnalyses = fullfile(pathDECODE,folderCoverAnalyses); 

pathChangeStatistic = fullfile(pathCoverAnalyses,'CoverChange');
pathMapStatistics = fullfile(pathCoverAnalyses,'MapStatistics');
pathConditionChange = fullfile(pathCoverAnalyses,'Disturbance');

%% Good Practice adjusted weights
confusionMatrixTB = readtable(fullfile(pathMapStatistics,'goodPractice_CoverMap.csv'),'ReadVariableNames',true,'VariableNamingRule','preserve');
realArea = confusionMatrixTB{ismember(confusionMatrixTB.Predicted,'realArea [%]'),2:end};
realAreaStd = confusionMatrixTB{ismember(confusionMatrixTB.Predicted,'realAreaStd [%]'),2:end};


figCover = figure('Position',[1 1 1600 1600]);
t = tiledlayout(3,2,'TileSpacing','loose','Padding','compact');
%% ----------------------------------------------------------Plot the agent----------------------------------------------------------
% Load the data
agentStatistics = readtable(fullfile(pathChangeStatistic,'coverConversionValidationSample_Driver.xlsx'),'VariableNamingRule','preserve');

agentTypes = {'Press disturbance','Extreme weather events','Recovery post events','Direct human activities','Direct human activities'}; % Six major agents
agentColors = [0,45,70;80,40,0;0,60,50;95,90,25;95,90,25]/100;

coverChangeAgentsOrigin = readtable(fullfile(pathChangeStatistic,'coverChangeAgentBootstrap.csv'),'ReadRowNames',true,'VariableNamingRule','preserve');
gapRow = array2table(nan(1,5),'variablenames',coverChangeAgentsOrigin.Properties.VariableNames);
gapRow.agent = 'Press disturbance';
gapRow.uniqueID = 11;
coverChangeAgents = [coverChangeAgentsOrigin(1:4,:); gapRow;gapRow;coverChangeAgentsOrigin(5:8,:);gapRow;gapRow;coverChangeAgentsOrigin(9:end,:)];

means = coverChangeAgents.mean;
stds = coverChangeAgents.std;
uniqueIDs = coverChangeAgents.uniqueID;


% Plot the bar
ax1 = nexttile(1);
agentBar = bar(means','FaceColor','flat','HandleVisibility','off');
displayColor = zeros(length(means),3);
displayNames = cell(1,length(uniqueIDs));
for i_u = 1:length(uniqueIDs)
    displayColor(i_u,:) = agentColors(mod(uniqueIDs(i_u),10),:);
    displayNames(i_u) = cellstr(agentTypes{mod(uniqueIDs(i_u),10)});
end
agentBar.CData = displayColor;
agentBar.EdgeAlpha = 0.8;
hold on;
%% To display the legend with a virtual bar
legendBar = bar([(1:length(agentTypes))/1000000;(1:length(agentTypes))/10000000]);
displayColor = zeros(length(agentTypes),3);
displayNames = cell(1,length(agentTypes));
for i_u = 1:length(agentTypes)
    displayColor(i_u,:) = agentColors(i_u,:);
    displayNames(i_u) = cellstr(agentTypes{i_u});
    legendBar(i_u).FaceColor = displayColor(i_u,:);
end

legend(legendBar(1:4),displayNames{1:4},'Location','southwest','box','off');
hold on;

errorbar(means,stds,'LineStyle','none','Color','k','LineWidth',2,'HandleVisibility','off');
lowerRange = round(min(means),-1);
upperRange = round(max(means),-1)+10;
ylim([lowerRange upperRange]);

% Rectangular
% Y = [upperRange,upperRange,lowerRange,lowerRange];
% 
% X = [0.5,4.49,4.49,0.5];
% patch('XData',X,'YData',Y,'EdgeColor',coverColors{1},'FaceColor','None','LineWidth',2,'LineStyle','--','HandleVisibility','off');
% 
% X = [4.51,8.49,8.49,4.51];
% patch('XData',X,'YData',Y,'EdgeColor',coverColors{2},'FaceColor','None','LineWidth',2,'LineStyle','--','HandleVisibility','off');
% 
% X = [8.51,11.5,11.5,8.51];
% patch('XData',X,'YData',Y,'EdgeColor',coverColors{3},'FaceColor','None','LineWidth',2,'LineStyle','--','HandleVisibility','off');

xticks([2.5 8 14]);
xticklabels({'Tidal marsh','Mangrove','Tidal flats'});
xlim([-0.5 15.51]);

ylabel('Cover changes (%)');
tg = text(0,8,'Gain','FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
tl = text(0,-25,'Loss','FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
tg.Rotation = 90;
tl.Rotation = 90;
% box off;
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');

% LimitsX = xlim; LimitsY = ylim;
% title('(a) Cover change agent','HorizontalAlignment','left', 'position', [LimitsX(1), LimitsY(2)]);

ttl = title('a','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
% ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

%% -----------------------------------------------------------Plot annual cover changes-----------------------------------------------------------------
ax2 = nexttile(2);
% Load the data
intervals = [1];
for i_i = 1:length(intervals)
    interval = intervals(i_i);
    % interval = 5;
    if interval==1
        years = 1986:interval:2020;
    else
        years = [1986,1990:interval:2020];
    end
    
    countCoverChange = fullfile(pathChangeStatistic,sprintf('tileBasedCoverConverstionCount_%d_Years.mat',interval));
    statistics(i_i).CoverChange = readtable(replace(countCoverChange,'mat','csv'),'ReadVariableNames',true,'VariableNamingRule','preserve');
    statistics(i_i).Interval = interval;
    statistics(i_i).Years = years;
end

%-------------------------------------------------------------------
%-------------------------------------------------------------------


hold on;
for i_i = 1:length(intervals)
    interval = statistics(i_i).Interval;
    accumulatedCoverChange = statistics(i_i).CoverChange;
    years = statistics(i_i).Years(2:end);
    if interval == 1
        lineType = '-';
        markerType = 'o';
        str = '(Annual)';
    else
        lineType = '--';
        markerType = 's';
        str = '(5-year interval)';
    end
    plot(years,accumulatedCoverChange.GainsAdjust,'-','Color',changeColors{1},'lineWidth',2,'DisplayName','Annual gain');
   
    patch([years'; flip(years')], [accumulatedCoverChange.GainsAdjust - accumulatedCoverChange.GainsStd; flip(accumulatedCoverChange.GainsAdjust + accumulatedCoverChange.GainsStd)],...
        'b', 'FaceColor',changeColors{1}, 'FaceAlpha',0.25, 'EdgeColor','none','HandleVisibility','off');
    
    plot(years,-accumulatedCoverChange.LossesAdjust,'-','Color',changeColors{2},'lineWidth',2,'DisplayName','Annual loss');
   
    patch([years'; flip(years')], [-accumulatedCoverChange.LossesAdjust - accumulatedCoverChange.LossesStd; flip(-accumulatedCoverChange.LossesAdjust + accumulatedCoverChange.LossesStd)],...
        'b', 'FaceColor',changeColors{2}, 'FaceAlpha',0.25, 'EdgeColor','none','HandleVisibility','off');
    
    % errorbar(years,accumulatedCoverChange.GainsAdjust,accumulatedCoverChange.GainsStd,'Color',changeColors{1},'LineStyle',lineType,'LineWidth',2,'DisplayName',sprintf('Gain'));
    % errorbar(years,-accumulatedCoverChange.LossesAdjust,accumulatedCoverChange.LossesStd,'Color',changeColors{2},'LineStyle',lineType,'LineWidth',2,'DisplayName',sprintf('Loss'));
%     errorbar(years,accumulatedCoverChange.NetChangeAdjust,accumulatedCoverChange.NetChangeStd,'Color',changeColors{3},'LineStyle',lineType,'LineWidth',3,'DisplayName',sprintf('Net change'));
    %yline(0,'-k','LineWidth',2,'HandleVisibility','off');
    % ylim([1800 2400])
    ylabel('Area (km^2)');
    xlabel('Year');
end

% %% Plot the condition change if needed
% % Load the data
% areaRate = 900/1000000;
% countContionChange = readtable(fullfile(pathConditionChange,sprintf('change_Patch_Statistics_81Tiles.csv')),'ReadVariableNames',true,'VariableNamingRule','preserve');
% % plot(countContionChange.Year(2:end),countContionChange.Changed(2:end)*areaRate,'Color',changeColors{4},'LineStyle',lineType,'LineWidth',3,'DisplayName',sprintf('Gross change'));
% hold on;
% plot(countContionChange.Year,countContionChange.Disturbed*areaRate,'Color',changeColors{5},'LineStyle','--','LineWidth',3,'DisplayName',sprintf('Disturbed'));
% plot(years,accumulatedCoverChange.GainsAdjust-accumulatedCoverChange.LossesAdjust,'Color',changeColors{5},'LineStyle','-','LineWidth',3,'DisplayName',sprintf('Cover changes'));
% title('(b) Annual changes');

% -------------------------------------
xlim([1985 2020]);
ylim([0 350]);
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,.0f';
lgd = legend('Location','Northwest','NumColumns',1,'Box','off');
lgd.ItemTokenSize = [50,15];
box on;

ttl = title('b','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.12; % use negative values (ie, -0.1) to move further left
% ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
%-------------------------------------------------------------------
%-------------------------------------------------------------------



%% ----------------------------------------------------------Plot subclass annual cover areas----------------------------------------------------------
% Load data
areaEstimates = readtable(fullfile(pathMapStatistics,'areaAdjustment.csv'),'VariableNamingRule','preserve','ReadVariableNames',true);

% Select the cover statistics every several years
if sum(ismember(intervals,5))
    intervalYr = [1986,1990:5:2020];
    areaEstimatesIntervals = areaEstimates(ismember(areaEstimates.Year,intervalYr),:);
end
yrs = areaEstimates.Year;
yrLength = length(areaEstimates.Year);
%-------------------------------------------------------------------
%-------------------------------------------------------------------
ax3 = nexttile;
hold on;
% shadedErrorBar(yrs,areaEstimates.adj_Tidalmarsh,areaEstimates.std_Tidalmarsh,'lineProps',{'-s','Color',coverColors{1},'lineWidth',2});

plot(yrs,areaEstimates.adj_Tidalmarsh,'-','Color',coverColors{1},'lineWidth',2,'HandleVisibility','off')
   
patch([yrs; flip(yrs)], [areaEstimates.adj_Tidalmarsh - areaEstimates.std_Tidalmarsh; flip(areaEstimates.adj_Tidalmarsh + areaEstimates.std_Tidalmarsh)],...
    'b', 'FaceColor',coverColors{1}, 'FaceAlpha',0.25, 'EdgeColor','none','HandleVisibility','off');
    
if sum(ismember(intervals,5))
    plot(intervalYr,areaEstimatesIntervals.adj_Tidalmarsh,'--o','Color','k','MarkerSize',12,'lineWidth',2);
end
xlabel('Year');
ylabel('Area (km^2)');
ylim([15000 20000]);
% ylim([0 20000]);
xlim([1985 2020]);
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,.0f';
if sum(ismember(5,intervals))
    lgd = legend({'Annual','5-year interval'},'Location','southwest');
    lgd.ItemTokenSize = [50,15];
end
box on;
ttl = title('c','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
% ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

yLoc = (max(ylim)-min(ylim))*0.95 + min(ylim);
xLoc = (max(xlim)-min(xlim))*0.5 + min(xlim);
text(xLoc,yLoc,sprintf('%s','Tidal marsh'),'FontName',fontName,'FontSize',fontSize+4,'HorizontalAlignment','center');

%-------------------------------------------------------------------
%-------------------------------------------------------------------
ax4 = nexttile;
hold on;

plot(yrs,areaEstimates.adj_Tidalflats,'-','Color',coverColors{4},'lineWidth',2,'HandleVisibility','off')
   
patch([yrs; flip(yrs)], [areaEstimates.adj_Tidalflats - areaEstimates.std_Tidalflats; flip(areaEstimates.adj_Tidalflats + areaEstimates.std_Tidalflats)],...
    'b', 'FaceColor',coverColors{4}, 'FaceAlpha',0.25, 'EdgeColor','none','HandleVisibility','off');
    
if sum(ismember(intervals,5))
    plot(intervalYr,areaEstimatesIntervals.adj_Tidalflats,'--o','Color','k','MarkerSize',12,'lineWidth',2);
end
xlabel('Year');
ylabel('Area (km^2)');
ylim([1000 3000]);
xlim([1985 2020]);
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,.0f';

if sum(ismember(5,intervals))
    lgd = legend({'Annual','5-year interval'},'Location','southwest');
    lgd.ItemTokenSize = [50,15];
end

box on;
ttl = title('d','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.12; % use negative values (ie, -0.1) to move further left
% ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

yLoc = (max(ylim)-min(ylim))*0.95 + min(ylim);
xLoc = (max(xlim)-min(xlim))*0.5 + min(xlim);
text(xLoc,yLoc,sprintf('%s','Tidal flats'),'FontName',fontName,'FontSize',fontSize+4,'HorizontalAlignment','center');


%-------------------------------------------------------------------
%-------------------------------------------------------------------
ax5 = nexttile(5);
hold on;

plot(yrs,areaEstimates.adj_Mangrove,'-','Color',coverColors{2},'lineWidth',2,'HandleVisibility','off')
   
patch([yrs; flip(yrs)], [areaEstimates.adj_Mangrove - areaEstimates.std_Mangrove; flip(areaEstimates.adj_Mangrove + areaEstimates.std_Mangrove)],...
    'b', 'FaceColor',coverColors{2}, 'FaceAlpha',0.25, 'EdgeColor','none','HandleVisibility','off');


if sum(ismember(intervals,5))
    plot(intervalYr,areaEstimatesIntervals.adj_Mangrove,'--o','Color','k','MarkerSize',12,'lineWidth',2);
end
xlabel('Year');
ylabel('Area (km^2)');
ylim([2000 2400]);
xlim([1985 2020]);
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,.0f';
if sum(ismember(5,intervals))
    lgd = legend({'Annual','5-year interval'},'Location','southwest');
    lgd.ItemTokenSize = [50,15];
end

box on;

ttl = title('e','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
% ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

yLoc = (max(ylim)-min(ylim))*0.95 + min(ylim);
xLoc = (max(xlim)-min(xlim))*0.5 + min(xlim);
text(xLoc,yLoc,sprintf('%s','Mangrove'),'FontName',fontName,'FontSize',fontSize+4,'HorizontalAlignment','center');

%-------------------------------------------------------------------
%-------------------------------------------------------------------
ax6 = nexttile(6);
hold on;

plot(yrs,areaEstimates.adj_Dieback,'-','Color',coverColors{3},'lineWidth',2,'HandleVisibility','off')
   
patch([yrs; flip(yrs)], [areaEstimates.adj_Dieback - areaEstimates.std_Dieback*10; flip(areaEstimates.adj_Dieback + areaEstimates.std_Dieback*10)],...
    'b', 'FaceColor',coverColors{3}, 'FaceAlpha',0.25, 'EdgeColor','none','HandleVisibility','off');

if sum(ismember(intervals,5))
    plot(intervalYr,areaEstimatesIntervals.adj_Dieback,'--o','Color','k','MarkerSize',12,'lineWidth',2);
end
xlabel('Year');
ylabel('Area (km^2)');
ylim([0 210]);
xlim([1985 2020]);
ylim([0 260]);
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.TickLabelFormat = '%,.0f';
if sum(ismember(5,intervals))
    lgd = legend({'Annual','5-year interval'},'Location','southwest');
    lgd.ItemTokenSize = [50,15];
end

box on;

ttl = title('f','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.12; % use negative values (ie, -0.1) to move further left
% ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  

yLoc = (max(ylim)-min(ylim))*0.95 + min(ylim);
xLoc = (max(xlim)-min(xlim))*0.5 + min(xlim);
text(xLoc,yLoc,sprintf('%s','Mangrove dieback'),'FontName',fontName,'FontSize',fontSize+4,'HorizontalAlignment','center');

exportgraphics(figCover,fullfile(pathFigure,'Figure2_Trend.jpg'),'Resolution',300);
saveas(figCover,fullfile(pathFigure,'Figure2_Trend.svg'));