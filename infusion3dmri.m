function [obj,im]=infusion3dmri(fname)
% INFUSION MRI ANALYSIS processes a stack of infusion, to analyze the
% infusion profile based on two sets of ground truth.
%
% @param
%   fName           - file name of tif stack
%   offset          - image intensity offset
%   threshold       - threshold of objects
%   threshold2      - threshold of gel
%   sizeObj         - size of small connected object to eliminate
%   erodeDist       - erode pixel distance of gel
%   concGroundTruth - ground truth concentrations
%
% @return
%   obj             - comprehensive information about objects and gel
%   fitP            - fit parameters (a,b) as y = (a)x + (b)
%   rsq             - R square
%
% @author
%   Yiming Kang, Cornell University
% @version
%   04/02/2014

% import tif image stack
% fname='/Users/schafferlab/Documents/MATLAB/yk_Data/JohnFoo/Infusion MRI/Stack.tif';
% fname='/Users/schafferlab/Documents/MATLAB/yk_Data/JohnFoo/Infusion MRI/Exam 195 Stacks/Exam 195 Series 15.tif';
offset=-32768;

infoim=imfinfo(fname);
if strcmp(infoim(1).PhotometricInterpretation,'BlackIsZero') == 1
    im=NaN(infoim(1).Width,infoim(1).Height,length(infoim));
    for i=1:length(infoim)
        im(:,:,i)=imread(fname,'Index',i);
    end
end
im=im+offset;

% identify objects and pixel intensity profile
threshold=30;
sizeObj=70;
erodeDist=5;

bw=im>threshold;
obj=bwconncomp(bw);
for i=1:obj.NumObjects
    if (length(obj.PixelIdxList{i})<sizeObj)
        obj.PixelIdxList{i}=[];
    end
end
obj.PixelIdxList=obj.PixelIdxList(~cellfun('isempty',obj.PixelIdxList));
obj.NumObjects=length(obj.PixelIdxList);

if (obj.NumObjects<3)
    fprintf('No infusion identified in:\n%s\n',fname);
elseif (obj.NumObjects>3)
    for i=1:obj.NumObjects
        if (length(obj.PixelIdxList{i})<sizeObj*2)
            obj.PixelIdxList{i}=[];
        end
    end
    obj.PixelIdxList=obj.PixelIdxList(~cellfun('isempty',obj.PixelIdxList));
    obj.NumObjects=length(obj.PixelIdxList);
end

obj.CenterPos=cell(1,obj.NumObjects);
obj.PixelIntList=cell(1,obj.NumObjects);
obj.PixelIntMean=zeros(1,obj.NumObjects);
obj.PixelIntStd=zeros(1,obj.NumObjects);
for i=1:obj.NumObjects
    [posX,posY,posZ]=ind2sub(obj.ImageSize,obj.PixelIdxList{i});
    obj.CenterPos{i}=round([mean(posX),mean(posY),mean(posZ)]);
    obj.PixelIntList{i}=im(obj.PixelIdxList{i});
    obj.PixelIntMean(i)=mean(obj.PixelIntList{i});
    obj.PixelIntStd(i)=std(obj.PixelIntList{i});
end

% identify gel and profile
threshold2=[-100,15];

bw2=im>threshold2(1)&im<threshold2(2);
bw2Dist=bwdist(imcomplement(bw2));
bw2=bw2Dist>erodeDist;
obj2=bwconncomp(bw2);

obj.GelPixelIdxList{1}=obj2.PixelIdxList{1};
obj.GelPixelIntList{1}=im(obj.GelPixelIdxList{1});
obj.GelPixelIntMean(1)=mean(obj.GelPixelIntList{1});
obj.GelPixelIntStd(1)=std(obj.GelPixelIntList{1});

% correlate intensity profile to concentration
concGroundTruth=[0,1.25,2.5];

posZList=zeros(1,obj.NumObjects);
for i=1:obj.NumObjects
    posZList(i)=obj.CenterPos{i}(3);
end
[foo,bubbleIdxList(1)]=min(posZList);
[foo,bubbleIdxList(2)]=max(posZList);
tmpList=1:obj.NumObjects;
infusionIdx=find(tmpList~=bubbleIdxList(1)&tmpList~=bubbleIdxList(2));
if (obj.PixelIntMean(bubbleIdxList(1))>obj.PixelIntMean(bubbleIdxList(2)))
    bubbleHighIdx=bubbleIdxList(1);
    bubbleLowIdx=bubbleIdxList(2);
else
    bubbleHighIdx=bubbleIdxList(2);
    bubbleLowIdx=bubbleIdxList(1);
end

intInfusion=obj.PixelIntMean(infusionIdx);
intGroundTruth=[obj.GelPixelIntMean(1),obj.PixelIntMean(bubbleLowIdx),obj.PixelIntMean(bubbleHighIdx)];

% calculate linear fit
fitP=polyfit(intGroundTruth,concGroundTruth,1);
fitF=polyval(fitP,intGroundTruth);
ssResid=sum((concGroundTruth-fitF).^2);
ssTotal=(length(concGroundTruth-1)*var(concGroundTruth));
rsq=1-ssResid/ssTotal;

% figure;
% plot(intGroundTruth,concGroundTruth,'x',intGroundTruth,fitF,'-');
% xlabel(upper('pixel intensity')); ylabel(upper('concentration'));
% annotation('textbox',[.35,.75,.1,.1],'String',...
%     sprintf('y = (%.5f) x + (%.5f)\nR^2 = %.5f',fitP(1),fitP(2),rsq));

% identify concentration profile
if (obj.NumObjects<3)
    obj.InfusionCtrSub=[];
    obj.InfusionIdxList=[];
    obj.InfusionIntList=[];
    obj.InfusionProfileX={[],[]};
    obj.InfusionProfileY={[],[]};
    obj.InfusionProfileZ={[],[]};
    obj.InfusionEdgeSub=[];
else
    [foo,tmpListIdx]=max(obj.PixelIntList{infusionIdx});
    [infusionCtr(1),infusionCtr(2),infusionCtr(3)]=ind2sub(obj.ImageSize,obj.PixelIdxList{infusionIdx}(tmpListIdx));
    [infusionSubList(:,1),infusionSubList(:,2),infusionSubList(:,3)]=ind2sub(obj.ImageSize,obj.PixelIdxList{infusionIdx});
    
    tmpIdxX=find(infusionSubList(:,2)==infusionCtr(2)&infusionSubList(:,3)==infusionCtr(3));
    tmpIdxY=find(infusionSubList(:,1)==infusionCtr(1)&infusionSubList(:,3)==infusionCtr(3));
    tmpIdxZ=find(infusionSubList(:,1)==infusionCtr(1)&infusionSubList(:,2)==infusionCtr(2));
    infusionProfileX=fitP(1)*obj.PixelIntList{infusionIdx}(tmpIdxX)+fitP(2);
    infusionProfileY=fitP(1)*obj.PixelIntList{infusionIdx}(tmpIdxY)+fitP(2);
    infusionProfileZ=fitP(1)*obj.PixelIntList{infusionIdx}(tmpIdxZ)+fitP(2);
    [foo,tmpCtrIdxX]=max(infusionProfileX);
    [foo,tmpCtrIdxY]=max(infusionProfileY);
    [foo,tmpCtrIdxZ]=max(infusionProfileZ);
    
    obj.InfusionCtrSub=infusionCtr;
    obj.InfusionIdxList=obj.PixelIdxList{infusionIdx};
    obj.InfusionIntList=obj.PixelIntList{infusionIdx};
    % profile divided into left and right (right includes center)
    obj.InfusionProfileX{1}=infusionProfileX(1:tmpCtrIdxX-1);
    obj.InfusionProfileX{2}=infusionProfileX(tmpCtrIdxX:end);
    obj.InfusionProfileY{1}=infusionProfileY(1:tmpCtrIdxY-1);
    obj.InfusionProfileY{2}=infusionProfileY(tmpCtrIdxY:end);
    obj.InfusionProfileZ{1}=infusionProfileZ(1:tmpCtrIdxZ-1);
    obj.InfusionProfileZ{2}=infusionProfileZ(tmpCtrIdxZ:end);
    
    % store image frame
    bwInfusion=zeros(size(im));
    bwInfusion(obj.InfusionIdxList)=1;
    bwInfusionEdge=logical(edge(bwInfusion(:,:,infusionCtr(3))));
    [obj.InfusionEdgeSub(:,1),obj.InfusionEdgeSub(:,2)]=ind2sub(size(bwInfusionEdge),find(bwInfusionEdge));
    
    % figure;
    % subplot(1,3,1); plot(infusionProfileX); xlabel('Profile X'); axis square;
    % subplot(1,3,2); plot(infusionProfileY); xlabel('Profile Y'); axis square;
    % subplot(1,3,3); plot(infusionProfileZ); xlabel('Profile Z'); axis square;
end