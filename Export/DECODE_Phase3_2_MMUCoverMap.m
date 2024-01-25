function DECODE_Phase3_2_MMUCoverMap(pathMap,maskImage,yr,varargin)    
   
    % optional
    
    p = inputParser;
    addParameter(p,'mmu', 5); % minimum mapping unit -->pixel numbers 900m2 * mmu(5) = 4500 m2 ~ 1 arce
    addParameter(p,'conn', 4); % Four directional or eight directional connection; 
    addParameter(p,'types',[6,5,4,1,3,2]); % Type number preferred to eliminate
    addParameter(p,'isRemoveIsolated',true); % whether to remove the isolated objects
    parse(p,varargin{:});
    MMU = p.Results.mmu;
    conn = p.Results.conn;
    types= p.Results.types;
    isRemoveIsolated = p.Results.isRemoveIsolated;
    % Load the wetland mask
    patchSize = 100;
    wetlandLayer = GRIDobj(maskImage);
    maskWetland = wetlandLayer.Z==1; % wetland

    % Load the classification image
    path_covermap = fullfile(pathMap,['cover_',num2str(yr),'.tif']);
    path_seivemap = fullfile(pathMap,['cover_MMU_',num2str(yr),'.tif']);
    
    if isfile(path_seivemap)
        fprintf('Skip %s\n',path_seivemap);
        return
    end

    mapValue = imread(path_covermap);
    
    %% Merge the subclasses (5-9) of others into others (5)
%     mapValue(mapValue>=5&mapValue<=9) = 6;
%     mapValue(mapValue==10) = 5; % Die back
    mapValue(mapValue>=6&mapValue<=10) = 6;
    fprintf('      Sieve the cover map: %s \n',path_covermap);
    %% Iteratively sieve each type not larger than MMU
    for iType = types
        ROIType = mapValue;
        ROIType(ROIType~=iType)=0;
        
        CC = bwconncomp(ROIType, conn);
        Stat = regionprops(CC, 'area','BoundingBox','SubarrayIdx');
        L = labelmatrix(CC);
        SmallIds = find([Stat.Area] < MMU);

        SmallObj = ismember(L, SmallIds);
    %         coverMap.Z = uint16(uint16(SmallObj).*L);
    %         GRIDobj2geotiff(coverMap, fullfile('G:\My Drive\DECODE\NEOutMap\Test','cover_small1_2020.tif'));

        bounding = [Stat.BoundingBox];
        bounding = reshape(bounding,[4,length(Stat)]);
        bounding = bounding';
        bounding(:,1:2) = floor(bounding(:,1:2));
        bounding(:,3:4) = bounding(:,3:4)+2;
        % Define the type based on the boundary pixels
        dilatedI = imdilate(SmallObj,ones(5,5));
        edgePixels = mapValue;
%         edgePixels(edgePixels==9)=0; % Dont modify to others
        edgePixels(dilatedI<1)=0;
        edgePixels(SmallObj>0)=0;
        tic
        for iObj = 1:length(Stat)
            if Stat(iObj).Area>=MMU
                Stat(iObj).Type = iType;
            else
                Icropped = imcrop(edgePixels,bounding(iObj,:));
                Icropped(Icropped==0)=[];
                if isempty(Icropped)
                    Stat(iObj).Type = iType; % Dont modify if no boundary pixel
%                     if removeIsolated
%                         Stat(iObj).Type = 9; % modify to background/remove it
%                     else
%                         Stat(iObj).Type = iType; % remain the small objects
%                     end
                else
                    Stat(iObj).Type  = mode(Icropped(:)); % seive into the neighbor type
                end
            end
        end
        if isempty(Stat)
            continue
        end
        typeNew = unique([Stat(:).Type]);
    %   typeNew(typeNew==iType)=[];
        for iNew = 1:length(typeNew)
            typeIds = find([Stat.Type]==typeNew(iNew));
            mapValue(ismember(L,typeIds)) = typeNew(iNew);
        end
        fprintf('      %0.0f seconds to process the type-%d \n',toc,iType);
    end
    
%      %% Remove the isolated patches
%      if isolatedDistance>0       
%         % Remove the middle size objects
%         ROIType = mapValue;
%         ROIType(~ismember(ROIType,[1,2,3]))=0;
%         
%         CC = bwconncomp(ROIType, 8);
%         Stat = regionprops(CC, 'area');
%         L = labelmatrix(CC);
%         patchIds = find([Stat.Area] < patchSize);
%         MiddleSizeObj = ismember(L, patchIds);
%         % Build the ocean water maps
%         background = mapValue; 
%         background(MiddleSizeObj) = 0;
%         patchDistance = bwdist(background);
%         
%         for i_p = 1:length(patchIds)
%             patchId = patchIds(i_p);    
%             patchObj = ismember(L, patchId);
%             if min(patchDistance(patchObj))>isolatedDistance
%                 maskValue(patchObj) = 0; % Minimum distance to the others
%             end
%         end
%         fprintf('Having removed middle sized isolated objects\n');
%     else
%         fprintf('Do not remove middle sized isolated objects\n');
%      end
    
    mapValue(~maskWetland) = 0; % Nonwetland region
%     mapValue(mapValue==255) = 20; % Unclassified
    outPutMap = wetlandLayer;
    outPutMap.Z = uint8(mapValue);
    GRIDobj2geotiff(outPutMap, path_seivemap);
    fprintf('    End of the process and save the processed map: %s \n',path_seivemap);
end
