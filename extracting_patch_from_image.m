ETH=imread('ETH_example_3D_model/dense3.nvm.cmvs/00/visualize/00000000.jpg');


% EXAMPLE FEATURE MATCH
% 5.44120297374 0.586922990216 4.70341271136        XYZ
% 173 156 136                                       RGB
% 2                                                 number of appearances
% 0 1578 387.673095703 286.252929688                imIndex,FeIndex,xy
% 1 1058 519.885009766 201.884277344 

% LEGEND
% <Point>  = <XYZ> <RGB> <number of measurements> <List of Measurements>
% <Measurement> = <Image index> <Feature Index> <xy>

ETH(388,286,1)=255;
ETH(388,286,2)=255;
ETH(388,286,3)=255;

imshow(ETH);