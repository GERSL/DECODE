function DECODE_Phase3_3_ObjectBasedChangeMap(pathMap,maskImage,year,varargin)
    p = inputParser;
    addParameter(p,'mmu', 5); % minimum mapping unit -->pixel numbers 900m2 * mmu(5) = 4500 m2 ~ 1 arce
    addParameter(p,'conn', 4); % Four directional or eight directional connection; 
    
    parse(p,varargin{:});
    MMU = p.Results.mmu;
    conn = p.Results.conn;

    % Load the classification image
    path_disturbancemap = fullfile(pathMap,['disturbance_',num2str(year),'.tif']);
    path_objectmap = fullfile(pathMap,['disturbance_MMU_',num2str(year),'.tif']);
    if isfile(path_objectmap)
        fprintf('Skip %s\n',path_objectmap);
        return
    end
    
    disturbanceImage = imread(path_disturbancemap);
    binaryChangeImage = disturbanceImage;
    %% Exclude the gradual disturbance map in the current phase
    binaryChangeImage(binaryChangeImage>=0&binaryChangeImage<=3) = 1; % 0-> non tidal wetland disturbance; 1-> pulse cover change; 2-> pulse condition change 3-> press cover change 4-> press condition change
    binaryChangeImage(binaryChangeImage~=1) = 0; % Binary map
    
    %% detect the small object
    CC = bwconncomp(binaryChangeImage, 8); % eight-connected; either altered to 4-connected
    Stat = regionprops(CC, 'area','BoundingBox','SubarrayIdx');
    L = labelmatrix(CC);
    largePatchIds = find([Stat.Area] >= MMU);
    largePatch = ismember(L, largePatchIds);    
    disturbanceImage(~largePatch) = 0; % Stable
    
    % Export the object_based_change_map
    maskLayer = GRIDobj(maskImage);
    disturbanceImage(maskLayer.Z~=1) = 0; % 0 is either stable or background
    maskLayer.Z = uint8(disturbanceImage);
    GRIDobj2geotiff(maskLayer,path_objectmap);
end

