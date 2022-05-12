function seiveSmallObjects(pathMap,maskImage,year,varargin)    
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'GRIDobj'));
    % optional
    
    p = inputParser;
    addParameter(p,'mmu', 5); % minimum mapping unit -->pixel numbers 900m2 * mmu(5) = 4500 m2 ~ 1 arce
    addParameter(p,'conn', 8); % Four directional or eight directional connection; 
    addParameter(p,'types',9); % Type number preferred to eliminate
    addParameter(p,'removeIsolated',false); % whether to remove the isolated small objects
    parse(p,varargin{:});
    mmu = p.Results.mmu;
    conn = p.Results.conn;
    types= p.Results.types;
    removeIsolated = p.Results.removeIsolated;
    % Load the wetland mask
    
    wetlandLayer = GRIDobj(maskImage);
    maskWetland = wetlandLayer.Z==1; % wetland

    % Load the classification image
    path_covermap = fullfile(pathMap,['cover_',num2str(year),'.tif']);
    path_seivemap = fullfile(pathMap,['cover_seive_',num2str(year),'.tif']);
    coverMap = GRIDobj(path_covermap);
    mapValue = coverMap.Z;
    
    fprintf('Begin to sieve the cover map: %s \n',path_covermap);
    %% Iteratively sieve each type not larger than MMU

    ROIType = mapValue;
    ROIType(ROIType~=types)=0;

    CC = bwconncomp(ROIType, conn);
    Stat = regionprops(CC, 'area','BoundingBox','SubarrayIdx');
    L = labelmatrix(CC);
    SmallIds = find([Stat.Area] < mmu);

    SmallObj = ismember(L, SmallIds);

    bounding = [Stat.BoundingBox];
    bounding = reshape(bounding,[4,length(Stat)]);
    bounding = bounding';

    dilatedI = imdilate(SmallObj,ones(5,5));
    edgePixels = mapValue;
    edgePixels(edgePixels==255)=0;
    edgePixels(dilatedI<1)=0;
    edgePixels(SmallObj>0)=0;
    tic
    for iObj = 1:length(Stat)
        if Stat(iObj).Area>=mmu
            Stat(iObj).Type = types;
        else
            Icropped = imcrop(edgePixels,bounding(iObj,:));
            Icropped(Icropped==0)=[];
            if isempty(Icropped)
                if removeIsolated
                    Stat(iObj).Type = 0; % modify to background/remove it
                else
                    Stat(iObj).Type = types; % remain the small objects
                end
            else
                Stat(iObj).Type  = mode(Icropped(:)); % seive into the neighbor type
            end
        end
    end
    if isempty(Stat)
       mapValue(~maskWetland) = 0;
       coverMap.Z = uint8(mapValue);
       GRIDobj2geotiff(coverMap, path_seivemap);
       fprintf('End of the process and no need to conduct seive: %s \n',path_seivemap);
       return
    end
    
    typeNew = unique([Stat(:).Type]);

    for iNew = 1:length(typeNew)
        typeIds = find([Stat.Type]==typeNew(iNew));
        mapValue(ismember(L,typeIds)) = typeNew(iNew);
    end
    fprintf('Having spent %0.0f seconds to process the type-%d \n',toc,types);

    mapValue(~maskWetland) = 0;
    coverMap.Z = uint8(mapValue);
    GRIDobj2geotiff(coverMap, path_seivemap);
    fprintf('End of the process and save the processed map: %s \n',path_seivemap);
end
