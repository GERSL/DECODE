path = fileparts(fileparts(mfilename('fullpath'))); 
pathStation = fullfile(path,'NOAA_CO-OPS_Present_and_Historical_Stations'); %path of shapefile
region = 'Alaska';
path_Tide = fullfile(path,'TideDownload',['Stations_',region]);

pathPre = fullfile(path,'TideDownload',['tidePre_',region]);
tideHour = 10; %The time of Landsat acquisition when to extract the tidal water level

stations = dir(fullfile(pathStation,'*shp'));
%stations = shaperead(fullfile(pathStation,stations.name)); % shapefile of the whole guage station

states = [];
% get station list according to the shapefile
% % Download observed hourly tide information and predicted information every
% % 6 minutes
% stationDownloaded = downloadTideInfoFromCOOPS(path,path_Tide,stations,states,(1984:2022),0);
% filterPredictWL(path_Tide,pathPre,tideHour);%extract the predicted water daily water level at 10 am
% % level by using the predicted values
% stationDownloaded = readmatrix('stationListDownloaded.csv');
% integrateWLInfo(pathPre,stationDownloaded,(1984:2022),['dailyWL_latlon_',region,'.csv']); %merge predicted into one csv
commandStr = sprintf(['module load gcc/5.4.0-alt sqlite/3.18.0 tcl/8.6.6.8606 zlib/1.2.11 java/1.8.0_162 mpi/openmpi/3.1.3 python/3.6.1;',...
                'module load proj/4.9.2 geos/3.5.0 mpi/openmpi/2.1.0 libjpeg-turbo/1.1.90 hdf5/1.8.19 netcdf/4.4.1.1-ompi-2.1.0 expat/2.2.0 gdal/2.2.1;',...
            '/apps2/python/3.6.1/bin/python3 /home/xiy19029/DECODE_github/prepareWL/interpolateWLInfo.py']);

system(commandStr);
function stationListFlag = downloadTideInfoFromCOOPS(path,folder_Tide,stations,states,years,i_plot)
    %%fileter the point list to be downloaded
    if isempty(states)
        lat = [stations.Y]';
        lon = [stations.X]';
        ids = {stations.id};
        ids = str2double(ids');
        stationList = [ids lat lon];
    else
        polygon1_y = states(1).X;     
        polygon1_x = states(1).Y;

        polygon1_x = polygon1_x.';
        polygon1_y = polygon1_y.';

        lat = [stations.Y];
        lon = [stations.X];
        ids = {stations.id};
        [in,on] = inpolygon(lon,lat,polygon1_y,polygon1_x);        % Logical Matrix
        inon = in | on;                                             % Combine ‘in’ And ‘on’
        idx = find(inon(:));                                     % Linear Indices Of ‘inon’ Points
        latcoord = lat(idx);                                        % X-Coordinates Of ‘inon’ Points
        loncoord = lon(idx);                                        % Y-Coordinates Of ‘inon’ Points
        ids = string(ids(idx));
        latcoord = latcoord';
        loncoord = loncoord';
        ids = str2double(ids');
        stationList = [ids latcoord loncoord];
        
        if i_plot >0
            figure(1)
            plot(lon, lat, '.')                                        % Plot All Points
            hold on
            plot(polygon1_y, polygon1_x, '.')                          % Plot Polygon
            plot(loncoord,latcoord, '*')                             % Overplot ‘inon’ Points
        %     hold off
        end
    end
    

    pathPreTide = fullfile(folder_Tide);%location to download the predicted tidal water level info
    
    if ~exist(pathPreTide) 
         mkdir(pathPreTide)         
    end
    
    flag = ones(length(stationList),1);
    if isfile(fullfile(pathPreTide,'stationListDownloaded.csv'))
        stationListFlag = readmatrix(fullfile(pathPreTide,'stationListDownloaded.csv'));
    else
        stationListFlag = [stationList,flag];
    end
    
    api = 'https://tidesandcurrents.noaa.gov/api/';
    % download water level from https://tidesandcurrents.noaa.gov/api-helper/url-generator.html
    %% download the predicted values
    for i = 1: length(stationList)
        id = num2str(stationList(i));
        if ~exist(fullfile(pathPreTide,id),'dir') 
            mkdir(fullfile(pathPreTide,id))         % ??????????????????‘Figure’
        else
            if exist(fullfile(pathPreTide,id,[num2str(years(end)),'.csv']),'file')
                continue
            end
        end
        fprintf('\n To download the station %d/%d - %s\n',i,length(stationList),id);
        for year = years % every nine years because the limitation of maximum 3650 days per download
            fprintf('%d - ',year);
            url = [api,'datagetter?begin_date=',num2str(year),'0101&end_date=',num2str(year),'1231&station=',id,'&product=predictions&datum=STND&time_zone=gmt&units=english&format=csv'];
            filename = fullfile(pathPreTide,id,[num2str(year),'.csv']);
            options = weboptions('Timeout',Inf);
            outfilename = websave(filename,url,options);
            if(year == years(1))
                s = dir(filename);
                a = s.bytes;
                if (a < 100)
                    % exclude the old stations
                    stationListFlag(i,end)=0;
                    rmdir(fullfile(pathPreTide,id),'s');
                    break;
                end
            end
        end
        writematrix(stationListFlag, 'stationListDownloaded.csv');
    end
    
end
function integrateWLInfo(pathPre,stationList,years,nameOut)
    dir_files=dir(fullfile(pathPre,'*.csv'));
    fileNames= {dir_files.name};
    
    doyWLInfo = [];
    lat_lon = [];
    for i = 1:length(fileNames)
        i
        fileName = char(fileNames(i));
        id = str2num(fileName(1:7));
        station = stationList(find(stationList(:,1)==id),:);
        lat = station(1,2);
        lon = station(1,3);
        wlInfo = readtable(fullfile(pathPre,fileName),'PreserveVariableNames',true);
        wlTemp = table2array(wlInfo(:,2));
        if ~isa(wlTemp(1),'double')
            fprintf('Invalid tide information for %d--%d,\n',i,id);
            continue
        end
        if i>1 && length(wlTemp)~=size(doyWLInfo,2)
            fprintf('Missing part of tide information for %d--%d,\n',i,id);
            continue
        end
        doyWLInfo = [doyWLInfo;wlTemp'];
        lat_lon = [lat_lon;id lat lon];
    end
    
    dailyList = datetime(years(1),1,1):datetime(years(end),12,31);
    [YList,~,~] = ymd(dailyList);
    doyList = day(dailyList,'dayofyear');
    doyList = YList*1000+doyList;
    nameDoy = cellstr(num2str(doyList'));
    headerDoy = ['id' 'lat' 'lon' nameDoy'];
    latlonWL = [lat_lon doyWLInfo];
    latlonWLTb = array2table(latlonWL,'VariableNames',headerDoy);
    writetable(latlonWLTb,nameOut);
end
function filterPredictWL(path_Tide,pathPre,tideHour)
    path_Tide = fullfile(path_Tide);
    if ~exist(pathPre) 
        mkdir(pathPre)         
    end

    stationsPre = dir(fullfile(path_Tide,'*'));
    tic
    for i = 3:length(stationsPre)
        
        id = stationsPre(i).name;
        preWLs = dir(fullfile(path_Tide,id,'*.csv'));
        if isempty(preWLs)
            continue
        end
        if exist(fullfile(pathPre,[id,'_10h.csv']))
            continue
        end
%         if str2num(id(1:7))~= 8417227
%             continue
%         end
        hourlyTide = [];
        for j = 1:length(preWLs)
            T = readtable(fullfile(path_Tide,id,preWLs(j).name),'PreserveVariableNames',true);
            hourlyTide = [hourlyTide;T];
        end
        dateTide = table2array(hourlyTide(:,1));
        numObserve = length(dateTide);
        [y,mt,d] = ymd(dateTide);
        doy = day(dateTide, 'dayofyear');
        doy = y*1000+doy;

        [h,mi,s] = hms(dateTide);
        index = find(h==tideHour & mi ==0);
        tideInfo = hourlyTide(index,:);
        doy = doy(index);
        doy = array2table(doy);
        writetable([tideInfo doy], fullfile(pathPre,[id,'_10h.csv']));
        fprintf('Having completed %d station - %s with in %0.2f mins\n',i-2,id,toc/60);
    end
end

