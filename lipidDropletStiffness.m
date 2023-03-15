% Lipid Droplet Stiffness Estimates from Compression
% Abigail Loneker, Wells Lab, UPenn
% March 2023

% Program for estimating lipid droplet stiffness based on the Eshelby
% equation. Parameters need to be modified to match experimental set up and
% aquired images including pixSize, gelStiffness, and gelStrain

clc;
clear;
close;


pathstem = "C:\Users\Abby\Downloads\2023_01_10_Lipid Droplets in PAA\2023_01_10_Lipid Droplets in PAA\Titchenell Collab\HFD 3"; % Change to match where images are stored
nuclearFiles = dir(fullfile(pathstem,'*.tif')); % Images should be .tif stacks with the uncompressed image first, bleach correction applied in imageJ
minSize = 5; % minimum size nucleus area in pixels
maxSize = 5000; % maximum size nucleus in pixels
pixSize = 562.56/1920; % pixel size in microns
gelStiffness = 10000; % Pa



for i=1:length(nuclearFiles)

    
    filestem = nuclearFiles(i).name;
    image = imread(fullfile(pathstem,filestem));

    info = imfinfo(fullfile(pathstem,filestem));
    numberOfImages = length(info);

    filename = extractBefore(filestem, '.tif');

    for k = 1:numberOfImages
        currentImage = imread(fullfile(pathstem,filestem), k, 'Info', info);
        imageStack(:,:,k) = currentImage; 
    end
    
    stiffnessTable = table;

    for k = 1:numberOfImages-1
        clearvars -except i image filestem info numberOfImages filename imageStack k match pathstem nuclearFiles minSize maxSize pixSize gelStiffness match currentImage
        lipid1 = imageStack(:,:,2);
        lipid2 = imageStack(:,:,k+1); 
        
        % Image registration: matches up fields of view to account for
        % drift
        imshowpair(lipid1, lipid2)
        [optimizer,metric] = imregconfig("multimodal");
        lipid2_reg = imregister(lipid2, lipid1, "rigid", optimizer, metric);
        imshowpair(lipid1, lipid2_reg);
        
        % Segmentation: gets label matrix for each image
        [lipid1numbered] = getLDlabel(lipid1, minSize, maxSize);
        [lipid2numbered] = getLDlabel(lipid2_reg, minSize, maxSize);

        % Object matching: filters lipid droplets and matches specific
        % droplets before and after compression
        overlap = lipid1numbered & lipid2numbered;
        pairs = [lipid1numbered(overlap), lipid2numbered(overlap)];
        pairs = unique(pairs, 'rows');
        pairs = pairs(all(pairs,2),:);
        [~,uidx] = unique(pairs(:,1),'stable');
        pairs = pairs(uidx,:);
        [~,uidx] = unique(pairs(:,2),'stable');
        pairs_clean = pairs(uidx,:);
        lipid2renum = zeros(size(currentImage));
        lipid1renum = zeros(size(currentImage));
        for r = 1:size(pairs_clean,1)
            lipid2renum(lipid2numbered == pairs_clean(r,2)) = r;
            lipid1renum(lipid1numbered == pairs_clean(r,1)) = r;
        end
        
        if k == 1
        outputfileName1 = ['Masks\LDmask_',filestem,'_',num2str(i),'.tif'];
        outputPath1 = fullfile(pathstem,outputfileName1);
        imwrite(lipid1renum, outputPath1, 'WriteMode', 'overwrite','Compression','none');
        outputfileName2 = ['Masks\LDmask_',filestem,'_',num2str(i),'_', num2str(k) 'kPa.tif'];
        outputPath2 = fullfile(pathstem,outputfileName2);
        imwrite(lipid2renum, outputPath2, 'WriteMode', 'overwrite','Compression','none');
        else
        outputfileName = ['Masks\LDmask_',filestem,'_',num2str(i),'_', num2str(k) 'kPa.tif'];
        outputPath = fullfile(pathstem,outputfileName);
        imwrite(lipid2renum, outputPath, 'WriteMode', 'overwrite','Compression','none');
        end
        
        % Shape analysis of lipid droplets
        lipid1stats = getLDstats(lipid1renum);
        lipid2stats = getLDstats(lipid2renum);
        
        % Gel strain calculated based on weight applied 
        gelStrain = 614.532287*k/10000; % Change depending on weight applied
        gelStress = gelStiffness * gelStrain;
        
        areaMicrons2 = lipid2stats.Area.*pixSize^2;
        areaMicrons1 = lipid1stats.Area.*pixSize^2;
        
        radii2 = lipid2stats.EquivDiameter/2;
        radii1 = lipid1stats.EquivDiameter/2;
        radiiMicrons2 = lipid2stats.EquivDiameter/2*pixSize;
        radiiMicrons1 = lipid1stats.EquivDiameter/2*pixSize;
        deltaR = radii1-((radii1.^3)./(radii2.^2));
        LDstrain = deltaR./radii1;
        LDstiffness2 = ((gelStrain./LDstrain.*5.*gelStiffness)-(3.*gelStiffness))/2;
        
        lipid1stats = addvars(lipid1stats, areaMicrons1, radiiMicrons1);
        lipid2stats = addvars(lipid2stats, areaMicrons2, radiiMicrons2, LDstiffness2);

        outputfilename = ['LDmeasurements_', num2str(filename), '_', num2str(k),'.xls'];
        outputPath = fullfile(pathstem, outputfilename);
        
        writetable(lipid1stats,outputPath,'Sheet','Lipid1');
        writetable(lipid2stats,outputPath,'Sheet','Lipid2');
             
    end


end

function LDnumbered = getLDlabel(image, minSize, maxSize)

        smoothedlipid1 = imgaussfilt(image,2);
        T = adaptthresh(smoothedlipid1, 0.1);
        BW = imbinarize(smoothedlipid1,T);
        lipidMask1 = imfill(BW,'holes');
        lipidMaskSeg1 = imclearborder(lipidMask1);
        lipidMaskSeg1 = watershedSegmentation(lipidMaskSeg1);
        lipidMaskSeg1 = bwconncomp(lipidMaskSeg1);
        L = labelmatrix(lipidMaskSeg1);
    
        stats = regionprops(L,'Area',"Circularity");
        idx = find(([stats.Area] > minSize) & ([stats.Area]< maxSize) & ([stats.Circularity] > 0.8)); 
        L2 = ismember(L,idx);
    
        LDnumbered = bwconncomp(L2);
        LDnumbered = labelmatrix(LDnumbered);
end

function [stats] = getLDstats(LDmask)
        
        
        stats = regionprops('table', LDmask, 'Area', 'Centroid','Perimeter', 'EquivDiameter' ,'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'ConvexArea' );
    
        areas = stats.Area;
        perims = stats.Perimeter;
        circularity = 4*pi*areas./(perims.^2);
        majorAxis = stats.MajorAxisLength;
        minorAxis = stats.MinorAxisLength;
        aspectRatio = majorAxis./minorAxis;
        roundness = 4*areas./(pi*majorAxis.^2);
        convexAreas = stats.ConvexArea;
        solidity = areas./convexAreas;
        stats = addvars(stats, circularity, aspectRatio, roundness, solidity);
end

