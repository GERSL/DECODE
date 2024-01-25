function ARDTiles = divideZones(ARDTiles,region)
    Tiles = char(ARDTiles);
    h = str2num(Tiles(:,2:4));
    v = str2num(Tiles(:,end-2:end));
    switch region
        case 'Northwest'
            idx = h<10 & v<=3;
        case 'West'
            idx = h<10 & v>3 & v<=8;    
        case 'Southwest'
            idx = h<10 & v>8;
        case 'Texas' 
            idx = h>10 & h<=17 & v>15;
        case 'WestLouisiana'
            idx = h>17 & h<=19 & v>15;    
        case 'LouisianaDelta'
            idx = h>19 & h<=21 & v>15;
        case 'South' 
            idx = h>21 & v==16;  
        case 'NorthFlorida'
            idx = h>24 & v==17;    
        case 'CentralFlorida'
            idx = h>24 & v==18;    
        case 'SouthFlorida'
            idx = h>24 & v>18; 
        case 'Southeast'
            idx = h>24 & v<=15 & v>=11;   
        case 'East'
            idx = h>24 & v<11 & v>=7;   
        case 'Northeast'
            idx = h>24 & v<7;          
    end
    ARDTiles = ARDTiles(idx,:);
end