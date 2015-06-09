% INFUSION ANALYSIS processes image series of infusion. The workflow is
% the following: 
%
% - load specific directory of infusion images
% - crop the stack according to the last frame (draw rectangle and double
%   click the area)
% - merge cap of infusion object subtracted by needle object using
%   morphological close 
% - binarize the stack with absolute threshold level (range 0:1)
% - reduce noise by deleting small connected component objects
% - detect the edge of infusion wave front using Canny edge detector
% - find the center of needle tip, the center of mass of infusion at each
%   frame, and compute the infusion rates (averaged rate, rates in N,S,E,W
%   directions)
% - save output in matlab file 'Data.mat'
% - visualize data
%
%   figure 1 - coordinates of needle tip and center of mass, and infusion
%              wave fronts at sampled frame 
%   figure 2 - averaged infusion rate over time
% 
% @param
%   pathName        - directory of image series
%   sizeImClose     - kernel size of morphological close operator
%   levelThreshold  - level of thresholding
%   secPerFrame     - second per frame
%   mmPerPixel      - minimeter per pixel
%
% @return
%   ctrInfusion
%   ctrNeedleTip
%   edgeInfusion
%   radiusInfusion
%   rateInfusion
% 
% @author 
%   Yiming Kang, Cornell University
% @version 
%   01/15/2014

clear; close all; clc;

%% PREPROCESSING
% load directory and convert all rgb png to 8-bit grayscale
pathName=uigetdir('*.*'); % @param
fileNameList=dir([pathName,'/','*.png']);
numFrames=size(fileNameList,1);

handle=waitbar(0/numFrames,['Loading Images: 0/',num2str(numFrames)]);
im0=cell(numFrames,1);
for i=1:numFrames
    if mod(i,20)==0
        waitbar(i/numFrames,handle,['Loading Images: ',num2str(i),'/',num2str(numFrames)]);
    end
    im0{i}=imread([pathName,'/',fileNameList(i).name],'png');
end
delete(handle);

% choose region of interest 
figure('name','Crop Image');
imshow(im0{numFrames});
pos=uint16(wait(imrect)); close; 
pause(0.5);

% subtract frames and close needle gap
im=im0;
im{1}=imcomplement(rgb2gray(im{1}(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),:)));
sizeImClose=10; % @param
se=strel('disk',sizeImClose);
handle=waitbar(0,['Cropping Images: 0/',num2str(numFrames)]);
for i=2:numFrames
    if mod(i,50)==0
        waitbar(i/numFrames,handle,['Cropping Images: ',num2str(i),'/',num2str(numFrames)]);
    end
    imCurrent=imcomplement(rgb2gray(im{i}(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),:)))-im{1};
    imCurrent=imclose(imCurrent,se);
    im{i}=imCurrent;
end
delete(handle);

%% SEGMENTATION
% binarize infusion and detect edge front
bwIm=cell(numFrames,1);
edgeInfusion=cell(numFrames,1);
levelThreshold=0.25; % @param
handle=waitbar(0/numFrames,['Processing Images: 0/',num2str(numFrames)]);
bwIm{1}=im2bw(im{1},0.5);
for i=2:numFrames
    if mod(i,20)==0
        waitbar(i/numFrames,handle,['Processing Images: ',num2str(i),'/',num2str(numFrames)]);
    end
    % binarize image
    bwCurrent=im2bw(im{i},levelThreshold);
    bwCurrent(1,:)=0; bwCurrent(end,:)=0;
    bwCurrent(:,1)=0; bwCurrent(:,end)=0;
    
    % delete small connected components
    CC=bwconncomp(bwCurrent,8);
    idxMaxObj=1;
    for j=1:CC.NumObjects-1
        if length(CC.PixelIdxList{idxMaxObj})<length(CC.PixelIdxList{j+1})
            idxMaxObj=j+1;
        end
    end
    for j=1:CC.NumObjects
        if j~=idxMaxObj
            bwCurrent(CC.PixelIdxList{j})=0;
        end
    end
    bwIm{i}=bwCurrent;
    
    % detect edge front
    edgeCurrent=logical(edge(bwCurrent,'canny'));
    [idxI,idxJ]=find(edgeCurrent);
    if ~isempty(idxI)
        edgeInfusion{i}(1:size(idxI),1)=uint16(idxI)+pos(2);
        edgeInfusion{i}(1:size(idxJ),2)=uint16(idxJ)+pos(1);
    end
end
delete(handle);

%% MEASUREMENT
% find tip of needle 
CC=bwconncomp(bwIm{1},8);
idxMaxObj=1;
for i=1:CC.NumObjects-1
    if length(CC.PixelIdxList{idxMaxObj})<length(CC.PixelIdxList{i+1})
        idxMaxObj=i+1;
    end
end
[idxI,idxJ]=ind2sub(CC.ImageSize,CC.PixelIdxList{idxMaxObj});
ITip=max(idxI);
listIdx=find(idxI==ITip);
JTipArray=zeros(length(listIdx),1);
for i=1:length(listIdx)
    JTipArray(i)=idxJ(listIdx(i));
end
JTipArray=sort(JTipArray);
ctrNeedleTip=[ITip+pos(2),(JTipArray(end)+JTipArray(1))/2+pos(1)];

% find center of mass of infusion
ctrInfusion=cell(numFrames,1);
for i=2:numFrames
    [currI,currJ]=find(bwIm{i});
    ctrInfusion{i}=[sum(currI)/length(currI)+pos(2),sum(currJ)/length(currI)+pos(1)];
end

% compute infusion front propagation
radiusInfusion=zeros(numFrames,5);
for i=1:numFrames
    if ~isempty(edgeInfusion{i})
        % average radius
        radiusCurrent=pdist2(double(ctrInfusion{i}),double(edgeInfusion{i}));
        radiusInfusion(i,1)=mean(radiusCurrent);
        % north and south directions
        idxList=find(edgeInfusion{i}(:,2)==ctrInfusion{i}(2));
        distList=double(ctrInfusion{i}(1))-double(edgeInfusion{i}(idxList,1));
        radiusInfusion(i,2)=max(distList);
        radiusInfusion(i,3)=-min(distList);
        % east and south directions
        idxList=find(edgeInfusion{i}(:,1)==ctrInfusion{i}(1));
        distList=double(ctrInfusion{i}(2))-double(edgeInfusion{i}(idxList,2));
        radiusInfusion(i,4)=-min(distList);
        radiusInfusion(i,5)=max(distList);
    end
end
rateInfusion=zeros(numFrames,5);
for i=2:numFrames-1
    for j=1:5
        rateInfusion(i,j)=radiusInfusion(i+1,j)-radiusInfusion(i,j);
    end
end

%% SAVE DATA
save([pathName,'/Data.mat'],'ctrNeedleTip','ctrInfusion','edgeInfusion',...
    'radiusInfusion','rateInfusion');

%% DATA VISUALIZATION
numGrad=10;
figure;
imshow(im0{numFrames}); hold on;
scatter(ctrNeedleTip(2),ctrNeedleTip(1),'r','+'); hold on;
scatter(ctrInfusion{numFrames}(2),ctrInfusion{numFrames}(1),'g','+'); hold on;
for i=1:numGrad
    j=round(numFrames/numGrad*i);
    handle=scatter(edgeInfusion{j}(:,2),edgeInfusion{j}(:,1),1,'fill'); 
    set(handle,'MarkerFaceColor',[i/numGrad,i/numGrad,0]);
    hold on;
end
hold off;
legend('Center of Needle Tip','Center of Mass of Infusion');

secPerFrame=5; % @param
mmPerPixel=1; % @param
widthSampling=10;
timeSample=secPerFrame*widthSampling*(1:floor((numFrames)/widthSampling));
rateSample=zeros(floor(numFrames/widthSampling),5);
for i=1:floor((numFrames)/widthSampling)
    for j=1:5
        rateSample(i,j)=mean(rateInfusion((i-1)*widthSampling+1:i*widthSampling,j));
    end
end
figure;
plot(timeSample,(rateSample(:,1))','Color',[0,0,0],'LineWidth',1.5);
legend('Average Infusion');
xlabel('Time Span (sec)'); ylabel('Infusion Rate (mm/sec)');
title('Infusion Rate Plot');
% optional to plot infusion rate in N,S,E,W directions
for i=2:5
    hold on; plot(timeSample,(rateSample(:,i))','Color',rand(1,3),'LineWidth',1);
end
hold off;
