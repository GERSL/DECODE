function DECODE_Phase3S_removeIsolatedTidalWetlandOverestimation(coverImage,patchImage,maskImage,varargin)    
   
    % optional    
    p = inputParser;
    addParameter(p,'mmu', 5); % minimum mapping unit -->pixel numbers 900m2 * mmu(5) = 4500 m2 ~ 1 arce
    addParameter(p,'conn', 4); % Four directional or eight directional connection; 
    addParameter(p,'types',[1,2,4]); % Type number preferred to eliminate
    addParameter(p,'patchSize',100); % whether to remove the isolated objects
    parse(p,varargin{:});
    MMU = p.Results.mmu;
    conn = p.Results.conn;
    types= p.Results.types;
    patchSize = p.Results.patchSize;
    
    wetlandLayer = GRIDobj(maskImage);
    maskWetland = wetlandLayer.Z==1; % wetland
    
%     isolatedDistance = 20;
        
    mapValue = imread(coverImage);
    patchImageValue = uint8(ones(size(mapValue))*6); % Default value as others
    
    openWater = mapValue==5; % Define the open water and non-open water
    patchImageValue(openWater) = 5;
    CC = bwconncomp(openWater, conn);
    Stat = regionprops(CC, 'area','BoundingBox','SubarrayIdx');
    L = labelmatrix(CC);
    SmallIds = find([Stat.Area] < patchSize);
    SmallObj = ismember(L, SmallIds);
    openWater(SmallObj) = 0;
    %Export Open water
    outPutMap = wetlandLayer;
    outPutMap.Z = uint8(openWater);
    GRIDobj2geotiff(outPutMap, replace(patchImage,'cover_Patch','LargeWater'));
    % Distance to the nearest water
    distanceToWater = round(bwdist(openWater));
    
    fprintf('      Begin to remove the isolated patch: %s (%d mins) \n',coverImage,round(toc/60));
    %% Iteratively sieve each type not larger than MMU
    for iType = types
        ROIType = mapValue==iType;
        CC = bwconncomp(ROIType, conn);
        Stat = regionprops(CC, 'area','BoundingBox','SubarrayIdx');
        
        [Stat(:).Distance] = deal(0); %Default value set as close to water
        
        L = labelmatrix(CC);
        SmallIds = find([Stat.Area] < patchSize);
      
        for i_p = 1:length(SmallIds)
            patchId = SmallIds(i_p);    
            patchObj = ismember(L, patchId);
            Stat(patchId).Distance = min(distanceToWater(patchObj));     
            
            if mod(i_p,100)==0
                fprintf('%05d/%05d (%d mins)\n',i_p,length(SmallIds),round(toc/60));
            end
        end
        IsolatedIds = find([Stat.Distance]-[Stat.Area]>0);
        IsolatedObj = ismember(L, IsolatedIds);
%         % Export the isolated objects
%         if ~isempty(Stat) && ~isempty(IsolatedIds)
%             outPutMap.Z = uint16(uint16(IsolatedObj).*uint16(L));
%             GRIDobj2geotiff(outPutMap, replace(patchImage,'cover_Patch',['Isolated_T',num2str(iType)]));
%         end
        % Update the cover maps
        ROIType(IsolatedObj) = 0;
        patchImageValue(ROIType) = iType;
          
        fprintf('      %0.0f mins to process the type- %d with %d objects\n',round(toc/60),iType,length(Stat));
    end
    
    %% Add the mangrove dieback
    patchImageValue(mapValue==3) = 3;

    %% Mask and export
    patchImageValue(~maskWetland) = 0; % Nonwetland region
    outPutMap = wetlandLayer;
    outPutMap.Z = uint8(patchImageValue);
    GRIDobj2geotiff(outPutMap, patchImage);
    fprintf('    End of the process and save the processed map: %s \n',patchImage);
end
