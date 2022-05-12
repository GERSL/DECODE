function outputImagePath = clipRasterRIOWarp(inputImagePath, imgPathLike, outputImagePath)
%CLIPRASTERRIOWARP is to use Python function to clip raster to same extent
%with same resolution
    try
        commandStr = sprintf(['module load gcc/5.4.0-alt sqlite/3.18.0 tcl/8.6.6.8606 zlib/1.2.11 java/1.8.0_162 mpi/openmpi/3.1.3 python/3.6.8;',...
            '/apps2/python/3.6.8/bin/python3 /home/xiy19029/DECODE_github/prepareWL/cust_warp.py %s %s --like %s'],...
            inputImagePath, outputImagePath, imgPathLike);
%         commandStr = sprintf('. $HOME/conda/etc/profile.d/conda.sh; conda activate; rio warp %s %s --like %s; conda deactivate;',...
%             inputImagePath, outputImagePath, imgPathLike);
        system(commandStr);
    catch
        outputImagePath = [];
        warning('Error: Clipping raster for %s \n', inputImagePath);
    end
map

end

