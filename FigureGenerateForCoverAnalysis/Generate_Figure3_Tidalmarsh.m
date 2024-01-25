%Aggregate the cover maps into CONUS wide with specific window size
%(default value is 10*10)

%% This paper displayed the migration of the tidal marsh along different elevation range
clc
clear
close all

fontName = 'Arial';
fontSize = 12;
c = parula(36);

% calculate the tidal marsh distribution every 0.2 DEM values
DEMRanges = -10:20; % 10 times of the original elevation value
DEMRangesTB = array2table([DEMRanges';max(DEMRanges)+1],'VariableNames',{'DEMRange'});
% Define the landward by using DEM, assuming that upland direction has higher altitude values

folderpath_mfile = fileparts(fileparts(mfilename('fullpath')));
addpath(folderpath_mfile);
addpath(fullfile(folderpath_mfile, 'Packages/MTools')); % Good practice to adjust the area estimates
addpath(fullfile(folderpath_mfile,'Packages/ktaub/')); % Sen's slope estimator

%% Output location
pathFigure = fullfile(fileparts(mfilename('fullpath')),'Figures');

years = 1986:2020;
% Forlder to save the results
folderCoverAnalyses = 'CoverAnalyses';

pathCoverAnalyses = fullfile(fileparts(fileparts(globalsets.PathDECODE)),folderCoverAnalyses);

pathMarshMigration = fullfile(pathCoverAnalyses,'tidalMarshElevation');
if ~exist(pathMarshMigration,'dir')
    mkdir(pathMarshMigration);
end

%% Load the good practice results for adjustment
pathMapValidation = fullfile(fileparts(fileparts(globalsets.PathDECODE)),folderCoverAnalyses,'MapValidation');
pathMapStatistics = fullfile(fileparts(fileparts(globalsets.PathDECODE)),folderCoverAnalyses,'MapStatistics');
zScore = 1.645; %1.96
areaRate = 900/1000000; %km2

%% Good Practice adjusted weights
confusionMatrixTB = readtable(fullfile(pathMapStatistics,'goodPractice_CoverMap.csv'),'ReadVariableNames',true,'VariableNamingRule','preserve');
realArea = confusionMatrixTB{ismember(confusionMatrixTB.Predicted,'realArea [%]'),2:end};
realAreaStd = confusionMatrixTB{ismember(confusionMatrixTB.Predicted,'realAreaStd [%]'),2:end};

%---------------------------------------------- Completete the adjustment --------------------------------

% Load the coast information
marshElevation = struct();
coasts = {'CONUS','South','East','West'};
for i_c = 1:length(coasts)
    marshElevation.(sprintf('count%s',coasts{i_c})) = readtable(fullfile(pathMarshMigration,sprintf('count%s_TidalMarsh.csv',coasts{i_c})));
end

%% Conduct error adjustment
% Either mangrove or dieback is involved
targetCovers = num2cell(1:6);

[Tidalmarsh,Mangrove,Dieback,Tidalflats,Openwater,Others] = deal(targetCovers{:});
coverTypes = {'TidalMarsh','Mangrove','Dieback','TidalFlats','Water','Others'};
patchMap = 'cover_Patch'; % map types % Cover maps

TidalmarshCoefficient = realArea(Tidalmarsh);
TidalmarshStdCoefficient = realAreaStd(Tidalmarsh);

%% Moving window to measure the decline speed change

for i_c = 1:length(coasts)
    coast = coasts{i_c};
    countMarshValue = round(marshElevation.(sprintf('count%s',coast)){:,2:end}*TidalmarshCoefficient)*areaRate;
    stdMarshValue= marshElevation.(sprintf('count%s',coast)){:,2:end}*TidalmarshStdCoefficient*areaRate;
    %% Gain or loss 
    RateSen = zeros(size(countMarshValue,1),4);
    for i_e = 2:length(DEMRanges) 
        % Exclude the extreme low elevation and high elevation ranges to
        % eliminate the commission error
        datain = [years',countMarshValue(i_e,:)'];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.1, 0);
        RateSen(i_e,:) = [h, sen, CIlower, CIupper];
    end
    RateSen = [DEMRangesTB,array2table(RateSen,'VariableNames',{'significant','sen','lower','upper'})];
    
    %% Strength of the trend
    MWArea = ((countMarshValue(:,10:end)-countMarshValue(:,1:end-9))/10);
    MWYr = years(5:end-5);
    MWSens = zeros(size(MWArea,1),4);
    for i_e = 2:length(DEMRanges) 
        % Exclude the extreme low elevation and high elevation ranges to
        % eliminate the commission error
        datain = [MWYr',MWArea(i_e,:)'];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] = ktaub(datain, 0.1, 0);
        MWSens(i_e,:) = [h, sen, CIlower, CIupper];
    end
    MWSens = [DEMRangesTB,array2table(MWSens,'VariableNames',{'significant','sen','lower','upper'})];
    MWArea = [DEMRangesTB,array2table(MWArea,'VariableNames',strcat('y',string(MWYr)))];
     
    marshElevation.(sprintf('Adjust%s',coast)) = countMarshValue;
    marshElevation.(sprintf('Std%s',coast)) = stdMarshValue;
    marshElevation.(sprintf('MW%s',coast)) = MWArea;
    marshElevation.(sprintf('Sen%s',coast)) = MWSens;
    marshElevation.(sprintf('Rate%s',coast)) = RateSen;
end

%% ------------------------------------------- Plot the results with histogram ---------------------------------------------------------------------
nHeight = 1;
figMarsh = figure('position',[1,1, 1800, 1600],'Visible','on');
t = tiledlayout(4,4,"TileSpacing",'loose','Padding','compact');

% regions = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};
coastalNames = {'Conterminous US','Gulf of Mexico','Atlantic Coast','Pacific Coast'};

for i_c =  1:length(coasts)
    %% Plot the annual areas of Mean from 1986 to 2020
    nexttile(i_c);
    coast = coasts{i_c};
    colororder({'k','k'});
    hold on
    %% -----------------------------------------------------------Plot annual marsh distribution-----------------------------------------------------------
    tidalMarshAnnualDistribution = marshElevation.(sprintf('Adjust%s',coast));    
    averageArea = tidalMarshAnnualDistribution(:,end);
    averageAreaStd = marshElevation.(sprintf('Std%s',coast))(:,end);
  
    errorbar(table2array(DEMRangesTB),averageArea,averageAreaStd,...
            'Marker','.','color','k','LineStyle','-','LineWidth',1.5,'CapSize',6,'DisplayName','Gain/loss trend');
    
    
    
    ylabel('Area (km^2)');
%     set(gca,'color', [0.8 0.8 0.8]);
    xticks([-10:5:20]);
    xticklabels({'-1','-0.5','0','0.5','1','1.5','2'});
    xlim([-10 20]);
    set(gca,'FontName',fontName,'fontSize',fontSize,'fontWeight','normal');
     
    ttl = title(sprintf('%s',char(96+i_c)),'FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
    % ax2.TitleHorizontalAlignment ='left';
    ttl.Units = 'Normalize'; 
    if i_c<3
        ttl.Position(1) = -0.2; % use negative values (ie, -0.1) to move further left
    else
        ttl.Position(1) = -0.23; % use negative values (ie, -0.1) to move further left
    end
     %     ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left';  

    text(0.5,1.05,sprintf('%s',coastalNames{i_c}),'FontName',fontName,'FontSize',fontSize,'HorizontalAlignment','center','Units','normalized');
  
    xlabel('Elevation (m)');    
    box on;
  
    %% Plot the annual changes
    nexttile(i_c+4);    
    hold on;
    variationArea = tidalMarshAnnualDistribution - averageArea;
    for i_y = 1:length(years)
        plot(table2array(DEMRangesTB),variationArea(:,i_y),'-','Color',c(i_y,:),...
            'MarkerSize',2,'MarkerFaceColor',c(i_y,:),'lineWidth',1.5,'HandleVisibility','off');
    end
    ylabel('Area (km^2)');
    xticks([-10:5:20]);
    xticklabels({'-1','-0.5','0','0.5','1','1.5','2'});
    xlim([-10 20]);
    set(gca,'FontName',fontName,'fontSize',fontSize,'fontWeight','normal');
     
    ttl = title(sprintf('%s',char(96+i_c+4)),'FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
    % ax2.TitleHorizontalAlignment ='left';
    ttl.Units = 'Normalize'; 
    if i_c<3
        ttl.Position(1) = -0.2; % use negative values (ie, -0.1) to move further left
    else
        ttl.Position(1) = -0.23; % use negative values (ie, -0.1) to move further left
    end
   
%     ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left';  
    xlabel('Elevation (m)');    
    box on;

    if i_c == length(coasts)
        colormap(c);
        cb = colorbar('Location','eastoutside','FontSize',fontSize,'FontName','Arial');
        cb.Ticks = [0.02,0.5,0.98];
        cb.TickLabels = {'1986','2003','2020'};    
       
    end

    %% -------------------------------------------------Plot the gain./loss rate---------------------------------------------------------------------
    nexttile (i_c+8);
    rateSen = marshElevation.(sprintf('Rate%s',coast));
    hold on;
    yline(0,'LineWidth',2,'LineStyle',':','HandleVisibility','off');
    
    errorbar(rateSen.DEMRange,rateSen.sen.*rateSen.significant,...
            (rateSen.sen-rateSen.lower).*rateSen.significant,...
            (rateSen.upper-rateSen.sen).*rateSen.significant,...
            'Marker','o','color','k','LineStyle','-','LineWidth',1.5,'CapSize',6,'DisplayName','Gain/loss trend');
    
    
    xlabel('Elevation (m)');  
    ylabel('Trend (km^2/yr^-^1)');
                 
    xlim([-10 20]);
    xticks([-10:5:20]);
    xticklabels({'-1','-0.5','0','0.5','1','1.5','2'});
    set(gca,'FontName',fontName,'fontSize',fontSize,'fontWeight','normal');
    
    ttl = title(sprintf('%s',char(96+i_c+8)),'FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
    % ax2.TitleHorizontalAlignment ='left';
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = -0.2; % use negative values (ie, -0.1) to move further left    
    ttl.HorizontalAlignment = 'left';  
    box on;

    %% -------------------------------------------------Plot the accelerated/decelerated---------------------------------------------------------------------
    nexttile (i_c+12);
    hold on;    
    rateSen = table2array(marshElevation.(sprintf('Rate%s',coast)));
    MWSens = table2array(marshElevation.(sprintf('Sen%s',coast)));
    MWSens(:,3:5) = MWSens(:,3:5);
    errorbar(MWSens(2:end-1,1),MWSens(2:end-1,2).*MWSens(2:end-1,3),...
        MWSens(2:end-1,2).*(MWSens(2:end-1,3)-MWSens(2:end-1,4)),MWSens(2:end-1,2).*(MWSens(2:end-1,3)-MWSens(2:end-1,5)),...
        'Marker','o','color','k','LineStyle','-','LineWidth',1.5,'CapSize',6,'DisplayName','Strength of trend');
    
    yline(0,'LineStyle','--','Color','k','HandleVisibility','off');

    minValue = min(MWSens(:,3:end),[],'all');
    maxValue = max(MWSens(:,3:end),[],'all');
    yrange = max(abs(minValue),abs(maxValue));

    xlabel('Elevation (m)');  
    ylabel('Strength (km^2/yr^-^2)');
%     set(gca,'color', [0.8 0.8 0.8]);
    
    ylim([-yrange yrange]);
             
    xlim([-10 20]);
    xticks([-10:5:20]);
    xticklabels({'-1','-0.5','0','0.5','1','1.5','2'});
    set(gca,'FontName',fontName,'fontSize',fontSize,'fontWeight','normal');
    
    ttl = title(sprintf('%s',char(96+i_c+12)),'FontName',fontName,'FontSize',fontSize+4,'FontWeight','bold');
    % ax2.TitleHorizontalAlignment ='left';
    ttl.Units = 'Normalize'; 
    
    ttl.Position(1) = -0.2; % use negative values (ie, -0.1) to move further left
%     ttl.Position(2) = 0.97; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left';  
    box on;
end

% Create textbox
annotation('textbox',...
    [0.98610864334629 0.156393851400643 0.076111111111111 0.0322108338116553],...
    'String',{'    Acc- gain','OR Dec- loss'},...
    'Rotation',90,...
    'LineStyle','none',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');

% Create arrow
annotation('arrow',[0.954441668209019 0.954444444444444],...
    [0.155686383601757 0.236456808199122]);

% Create arrow
annotation('arrow',[0.955552470849528 0.955],...
    [0.141026061493411 0.0556368960468521]);

% Create textbox
annotation('textbox',...
    [0.98499660682337 0.0582972188969832 0.0761111111111108 0.0322108338116553],...
    'String',{'    Dec- gain','OR Acc- loss'},...
    'Rotation',90,...
    'LineStyle','none',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off');

exportgraphics(figMarsh,fullfile(pathFigure,sprintf('Figure3_TidalMarsh.jpg')),'Resolution',300);
saveas(figMarsh,fullfile(pathFigure,sprintf('Figure3_TidalMarsh.svg')));


