%Aggregate the cover maps into CONUS wide with specific window size;
%default value is 10*10; 

clc
clear
close all

fontName = 'Arial';
fontSize = 14;
c = parula(36);
markerColor = {'#009e73','#e69f00'};
trendColor = {'#0072b2','#d55e00'};

folderpath_mfile = fileparts(fileparts(mfilename('fullpath')));
addpath(folderpath_mfile);
addpath(fullfile(folderpath_mfile, 'Packages/MTools'));
addpath(fullfile(folderpath_mfile,'Packages/ktaub/'));

%% Output location
pathFigure = fullfile(fileparts(mfilename('fullpath')),'Figures');
years = 1986:2020;

% Accuracy assessment info to conduct area adjustment
folderCoverAnalyses = 'CoverAnalyses';
pathCoverAnalyses = fullfile(fileparts(fileparts(globalsets.PathDECODE)),folderCoverAnalyses);

% Mangrove statistics at latitude direction
pathMangroveMigration = fullfile(pathCoverAnalyses,'mangroveMigration');
statisticAlongLat = readtable(fullfile(pathMangroveMigration,'MigrationAlongLatitude.csv'),'ReadVariableNames',true,'VariableNamingRule','preserve');
statisticAlongLatStd = readtable(fullfile(pathMangroveMigration,'MigrationAlongLatitudeSTD.csv'),'ReadVariableNames',true,'VariableNamingRule','preserve');

% Convert table to matrix
latsUnique = statisticAlongLat.Lat;
statisticAlongLat = statisticAlongLat{:,2:end};
statisticAlongLatStd = statisticAlongLatStd{:,2:end};

%% Load the good practice results for adjustment to calculate area bias 
pathMapStatistics = fullfile(fileparts(fileparts(globalsets.PathDECODE)),folderCoverAnalyses,'MapStatistics');
confusionMatrixTB = readtable(fullfile(pathMapStatistics,'goodPractice_CoverMap.csv'),'ReadVariableNames',true,'VariableNamingRule','preserve');
realArea = confusionMatrixTB{ismember(confusionMatrixTB.Predicted,'realArea [%]'),2:end};
realAreaStd = confusionMatrixTB{ismember(confusionMatrixTB.Predicted,'realAreaStd [%]'),2:end};

%% calculate the trend along latitude direction
sens = zeros(length(latsUnique),9);
temporalWindowSize = 10; % Trend calculate every 10 years
areaMWOfLat = zeros(length(latsUnique),length(years));

for i_lat = 1:length(latsUnique)
    datain = [years',statisticAlongLat(i_lat,:)'];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.1, 0);
    sens(i_lat,1) = h;
    sens(i_lat,2) = sen;
    sens(i_lat,3) = CIlower;    
    sens(i_lat,4) = CIupper;

    %% moving window to calculate the accelerated rate
    speedMW = zeros(length(years),1); 
    areaMW = zeros(length(years),1);
    for i_w = 1:length(years)-temporalWindowSize+1
        yrSpan = years(i_w:i_w+temporalWindowSize-1);
        idxSpan = ismember(years,yrSpan);        
        datainInWindow = datain(idxSpan,:);
        %         areaMW(i_w) = (datainInWindow(end,2)-datainInWindow(1,2))/windowSize; 
        % If use the slope within the window
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datainInWindow, 0.1, 0); 
        areaMW(i_w) = sen;     
    end
    
    areaMWOfLat(i_lat,:) = areaMW;
    
    datain = [years',areaMW];
    datain = datain((1:length(years)-temporalWindowSize+1),:);
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.1, 0);
    sens(i_lat,6) = h;
    sens(i_lat,7) = sen;
    sens(i_lat,8) = CIlower;  
    sens(i_lat,9) = CIupper;
end

trendAlongLat = [latsUnique,sens];
trendAlongLat = array2table(trendAlongLat,'VariableNames',{'Lat','SignificantChange','SenChange','SenChangeLower','SenChangeUpper',...                            
                            'Space','SignificantAreaGrowth','SenAreaGrowth','SenAreaGrowthLower','SenAreaGrowthUpper'});
statisticAlongLat = array2table(statisticAlongLat,'VariableNames',strcat('y',string(years)));
trendAlongLat = [trendAlongLat,statisticAlongLat];
% writetable(trendAlongLat,fullfile());

latitudeLimits = [25 30];
idxInLimit = find(trendAlongLat.Lat>=latitudeLimits(1)&trendAlongLat.Lat<=latitudeLimits(end));
trendAlongLatToPlot = trendAlongLat(idxInLimit,:);

%% Begin to plot
%% _______________________________________________________Plot the state-based results__________________________________________________________________
figMangrove = figure('Position',[1 1 1600 1000]);
t = tiledlayout(2,2,"TileSpacing","loose",'Padding','loose');

%% ------ Plot the rate of change compared to the value in 2020 -----------
ax1 = nexttile(1,[1,1]);
hold on;
latsToPlot = trendAlongLatToPlot.Lat;
mangroveAnnualDistribution = trendAlongLatToPlot{:,11:end};

averageArea = mangroveAnnualDistribution(:,end);
averageAreaStd = averageArea*realAreaStd(2)./realArea(2);

% averageArea = mean(tidalMarshAnnualDistribution,2);
% averageAreaStd = mean(marshElevation.(sprintf('Std%s',coast)),2);
errorbar(latsToPlot,averageArea,averageAreaStd,...
        'Marker','.','color','k','LineStyle','-','LineWidth',1.5,'CapSize',6,'DisplayName','Gain/loss trend');
ylabel('Area (km^2)');
xlim([latitudeLimits(1) latitudeLimits(2)]);
xticks([latitudeLimits(1):1:latitudeLimits(2)]);
xticklabels(string([latitudeLimits(1):1:latitudeLimits(2)]));
xlabel('Latitude (^o)');
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal'); 
box on;

% Add subfigure series number
ttl = title('a','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
ttl.Position(2) = 0.98; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';   


ax1 = nexttile(2,[1,1]);
hold on;
variationArea = mangroveAnnualDistribution - averageArea;
for i_y = 1:length(years)
    plot(latsToPlot,variationArea(:,i_y),'-','Color',c(i_y,:),...
        'MarkerSize',2,'MarkerFaceColor',c(i_y,:),'lineWidth',1,'HandleVisibility','off');
end

ylabel('Area (km^2)');
xlim([latitudeLimits(1) latitudeLimits(2)]);
xticks([latitudeLimits(1):1:latitudeLimits(2)]);
xticklabels(string([latitudeLimits(1):1:latitudeLimits(2)]));
xlabel('Latitude (^o)');
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal'); 

% Add subfigure series number
ttl = title('b','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
ttl.Position(2) = 0.98; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left'; 

c = parula(36);
colormap(c);
cb = colorbar('FontSize',fontSize,'FontName','Arial');
% set(gca,'color', [0.8 0.8 0.8]);
cb.Ticks = [0.02,0.5,0.98];
cb.TickLabels = {'1986','2003','2020'};    
box on
% 

%% -----------------------------------------Plot the annual gain/loss rate-------------------------------
ax = nexttile(3);

yline(0,'LineWidth',2,'LineStyle',':','HandleVisibility','off');
hold on;
errorbar(trendAlongLatToPlot.Lat,trendAlongLatToPlot.SenChange.*trendAlongLatToPlot.SignificantChange,...
        (trendAlongLatToPlot.SenChange-trendAlongLatToPlot.SenChangeLower).*trendAlongLatToPlot.SignificantChange,...
        (trendAlongLatToPlot.SenChangeUpper-trendAlongLatToPlot.SenChange).*trendAlongLatToPlot.SignificantChange,...
        'Marker','o','color','k','LineStyle','-','LineWidth',1.5,'CapSize',12,'DisplayName','Gain/loss trend');

xlabel('Latitude (^o)');
ylabel('Trend (km^2/yr^-^1)');       

leg = legend('Location','Northeast','Box','off','FontName','Arial','FontSize',fontSize,'FontWeight','normal');
leg.ItemTokenSize = [60,25];
ylim([-4 4]);

set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal'); 

xlim([25 30]);
box on;
ttl = title('c','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
ttl.Position(2) = 0.98; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left'; 


%% ----------------------------------------- Strength      -----------------------------------------
ax2 = nexttile(4);
MWSens = trendAlongLatToPlot{:,{'Lat','SignificantAreaGrowth','SenAreaGrowth','SenAreaGrowthLower','SenAreaGrowthUpper'}};
rateSen = trendAlongLatToPlot{:,{'Lat','SignificantChange','SenChange','SenChangeLower','SenChangeUpper'}};

hold on
errorbar(trendAlongLatToPlot.Lat,trendAlongLatToPlot.SenAreaGrowth.*trendAlongLatToPlot.SignificantAreaGrowth,...
        (trendAlongLatToPlot.SenAreaGrowth-trendAlongLatToPlot.SenAreaGrowthLower).*trendAlongLatToPlot.SignificantAreaGrowth,...
        (trendAlongLatToPlot.SenAreaGrowthUpper-trendAlongLatToPlot.SenAreaGrowth).*trendAlongLatToPlot.SignificantAreaGrowth,...
        'Marker','o','color','k','LineStyle','-','LineWidth',1.5,'CapSize',12,'DisplayName','Strength of Trend');
yline(0,'LineWidth',2,'LineStyle',':','HandleVisibility','off');

xlabel('Latitude (^o)');
ylabel('Strength (km^2/yr^-^2)');       
% leg = legend('Location','Northwest','Box','off','FontName',fontName,'FontSize',fontSize,'FontWeight','normal');
% leg.ItemTokenSize = [60,25];
set(gca,'FontName',fontName,'FontSize',fontSize,'FontWeight','normal'); 
ylim([-0.2 0.2]);
xlim([25 30]);
box on;
ttl = title('d','FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
% ax2.TitleHorizontalAlignment ='left';
ttl.Units = 'Normalize'; 
ttl.Position(1) = -0.15; % use negative values (ie, -0.1) to move further left
ttl.Position(2) = 0.98; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left'; 
lgd = legend('Location','northeast','FontSize',fontSize,'FontName','Arial','Box','off');
anArrow = annotation('arrow') ;
anArrow.Parent = gca;
anArrow.Position = [1, 0.2, 0, 0.1] ;


% Create arrow
annotation('arrow',[0.9375 0.9375],...
    [0.295000000000001 0.443000000000001],'LineWidth',2);

% Create arrow
annotation('arrow',[0.939375 0.93875],[0.27 0.126],'LineWidth',2);

% Create textbox
annotation('textbox',...
    [0.980625000000003 0.138 0.10218749714084 0.0629999999999998],...
    'String',{'    Dec- gain','OR Acc- loss'},...
    'Rotation',90,...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Arial',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.980000000000003 0.297000000000001 0.10218749714084 0.063],...
    'String',{'     Acc- gain','OR Dec- loss'},...
    'Rotation',90,...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Arial',...
    'FitBoxToText','off');

exportgraphics(figMangrove,fullfile(pathFigure,sprintf('Figure4_Mangrove.jpg')),'Resolution',300);
saveas(figMangrove,fullfile(pathFigure,sprintf('Figure4_Mangrove.svg')));

