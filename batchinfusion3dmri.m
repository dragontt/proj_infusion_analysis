% BATCH INFUSION MRI ANALYSIS batch processes a stack of infusion, to
% analyze the infusion profile based on two sets of ground truth.
%
% @author
%   Yiming Kang, Cornell University
% @version
%   04/02/2014

close all; clear; clc;

% data information
pname=uigetdir('*.*');
fnameList=dir([pname,'/','*.tif']);
numStacks=size(fnameList,1);

% parse batch analysis from infusion3dmri
handle=waitbar(0/numStacks,['PROCESSING DATA: 0/',num2str(numStacks)]);

infusionVolList=zeros(1,numStacks);
infusionProfileX=cell(1,numStacks);
infusionProfileY=cell(1,numStacks);
infusionProfileZ=cell(1,numStacks);
maxLengthProfileX=[0,0];
maxLengthProfileY=[0,0];
maxLengthProfileZ=[0,0];
infusionEdgeSub=cell(1,numStacks);
for i=1:numStacks
    waitbar(i/numStacks,handle,['PROCESSING DATA: ',num2str(i),'/',num2str(numStacks)]);
    fname=[pname,'/',fnameList(i).name];
    [obj,im]=infusion3dmri(fname);
    infusionVolList(i)=length(obj.InfusionIdxList);
    infusionProfileX{i}=obj.InfusionProfileX;
    infusionProfileY{i}=obj.InfusionProfileY;
    infusionProfileZ{i}=obj.InfusionProfileZ;
    for j=1:2
        if (maxLengthProfileX(j)<length(infusionProfileX{i}{j}))
            maxLengthProfileX(j)=length(infusionProfileX{i}{j});
        end
        if (maxLengthProfileY(j)<length(infusionProfileY{i}{j}))
            maxLengthProfileY(j)=length(infusionProfileY{i}{j});
        end
        if (maxLengthProfileZ(j)<length(infusionProfileZ{i}{j}))
            maxLengthProfileZ(j)=length(infusionProfileZ{i}{j});
        end
    end
    if (~isempty(obj.InfusionEdgeSub))
        infusionEdgeSub{i}=obj.InfusionEdgeSub;
    end
end
delete(handle);

infusionRateList=zeros(1,numStacks);
for i=1:numStacks
    if (i>1)
        infusionRateList(i)=infusionVolList(i)-infusionVolList(i-1);
    end
    infusionProfileX{i}{1}=[zeros(maxLengthProfileX(1)-length(infusionProfileX{i}{1}),1);infusionProfileX{i}{1}];
    infusionProfileX{i}{2}=[infusionProfileX{i}{2};zeros(maxLengthProfileX(2)-length(infusionProfileX{i}{2}),1)];
    infusionProfileY{i}{1}=[zeros(maxLengthProfileY(1)-length(infusionProfileY{i}{1}),1);infusionProfileY{i}{1}];
    infusionProfileY{i}{2}=[infusionProfileY{i}{2};zeros(maxLengthProfileY(2)-length(infusionProfileY{i}{2}),1)];
    infusionProfileZ{i}{1}=[zeros(maxLengthProfileZ(1)-length(infusionProfileZ{i}{1}),1);infusionProfileZ{i}{1}];
    infusionProfileZ{i}{2}=[infusionProfileZ{i}{2};zeros(maxLengthProfileZ(2)-length(infusionProfileZ{i}{2}),1)];
end

% plot figures
figure;
xAxisProfileX=-(maxLengthProfileX(2)-1):maxLengthProfileX(1);
xAxisProfileY=-maxLengthProfileY(1):(maxLengthProfileY(2)-1);
xAxisProfileZ=-maxLengthProfileZ(1):(maxLengthProfileZ(2)-1);
stringTmp=sprintf('FRAME # %d*',1:numStacks);
stringLegend=regexp(stringTmp,'*','split');
subplot(3,1,1);
for i=1:numStacks
    plot(xAxisProfileX,flipud([infusionProfileX{i}{1};infusionProfileX{i}{2}]),'Color',[i/numStacks,0,(numStacks-i)/numStacks]); hold on;
end
hold off; xlabel('PROFILE X'); ylabel('INFUSION CONCETRATION'); legend(stringLegend{:});
subplot(3,1,2);
for i=1:numStacks
    plot(xAxisProfileY,[infusionProfileY{i}{1};infusionProfileY{i}{2}],'Color',[i/numStacks,0,(numStacks-i)/numStacks]); hold on;
end
hold off; xlabel('PROFILE Y'); ylabel('INFUSION CONCETRATION'); legend(stringLegend{:});
subplot(3,1,3);
for i=1:numStacks
    plot(xAxisProfileZ,[infusionProfileZ{i}{1};infusionProfileZ{i}{2}],'Color',[i/numStacks,0,(numStacks-i)/numStacks]); hold on;
end
hold off; xlabel('PROFILE Z'); ylabel('INFUSION CONCETRATION'); legend(stringLegend{:});

figure;
subplot(2,1,1); plot(infusionRateList);
xlabel('VOLUME RATE OVER FRAME'); ylabel('INFUSION VOLUME');
subplot(2,1,2); plot(infusionVolList);
xlabel('VOLUME OVER FRAME'); ylabel('INFUSION VOLUME');

figure;
imagesc(histeq(im(:,:,obj.InfusionCtrSub(3)))); colormap gray; axis image; hold on;
for i=1:numStacks
    if (~isempty(infusionEdgeSub{i}))
        handle=scatter(infusionEdgeSub{i}(:,2),infusionEdgeSub{i}(:,1)); 
        set(handle,'MarkerFaceColor',[i/numStacks,0,(numStacks-i)/numStacks]);hold on;
    end
end
hold off;


