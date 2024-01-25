function maximumTidalWetlandDistribution(pathMap,path_maximumMap,maskImage,years) 
    % Load the wetland mask    
    wetlandLayer = GRIDobj(maskImage);
   
    % Export the patch map after seiving process with MMU and distance to open water

    for yr = years
        path_patchmap = imread(fullfile(pathMap,['cover_Patch_',num2str(yr),'.tif']));
        path_patchmap(path_patchmap>0&path_patchmap<=4) = 1; % binary existence of tidal wetlands
        path_patchmap(path_patchmap~=1) = 0;
        if ~exist('maximumExtent','var')
            maximumExtent = path_patchmap;
        else
            maximumExtent = maximumExtent|path_patchmap;
        end
    end
    
    wetlandLayer.Z = uint8(maximumExtent);
    GRIDobj2geotiff(wetlandLayer,path_maximumMap);
end