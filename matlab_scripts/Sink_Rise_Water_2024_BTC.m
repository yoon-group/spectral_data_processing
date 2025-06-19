clear all; clc; clf
%% ========================================================================
% 
%            SweeatWater - Rise 2024 dye tracing water samples
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
% The data files naming rule: "trip#_sampler#_bottle#_CorrectionData.txt"
% 
% =========================================================================

campaign = 'SinkRise2024'; % 'SinkRise2024' or 'BearSpring2023'

archiveDir = sprintf('/home/public/dyeTracingData/');

resultsDir = sprintf('%s/processed/%s/',pwd,campaign);
if ~exist(resultsDir,'dir'); mkdir(resultsDir); end

figDir = [resultsDir 'fig/samples/']; 
if ~exist(figDir,'dir'); mkdir(figDir); end

dataDir = sprintf('%s/%s/sampleData/',archiveDir,campaign );

sampleType = 'water'; % either {'water', 'eluent'}

%%  injection times
%   * RWT of 18 liters at 1.8% (1:55 ~ 2:26 pm) 
%   * uranine 19 liters at 4.8% (2:26 – 2:42 pm) 
t_SrBInj = datetime('2024-11-18T14:10','InputFormat','uuuu-MM-dd''T''HH:mm');
t_uranineInj = datetime('2024-11-18T14:32','InputFormat','uuuu-MM-dd''T''HH:mm');

%% dating samples
% trip 1 (Nov. 18)
%    - sampler 1: starting 18:00 PM sampling every 2 hours
%    - sampler 2: starting 18:30 PM sampling every 2 hours
%    - sampler 3: starting 19:00 PM sampling every 2 hours
%    - sampler 4: starting 19:30 PM sampling every 2 hours
% trip 2 (Nov. 20)
%    - sampler 1: starting 16:00 PM sampling every 2 hours
%    - sampler 2: starting 16:30 PM sampling every 2 hours
%    - sampler 3: starting 17:00 PM sampling every 2 hours
%    - sampler 4: starting 17:30 PM sampling every 2 hours
% trip 3 (Nov. 22)
%    - sampler 1: starting 13:00 PM sampling every 3 hours
%    - sampler 2: starting 13:45 PM sampling every 3 hours
%    - sampler 3: starting 14:30 PM sampling every 3 hours
%    - sampler 4: starting 15:15 PM sampling every 3 hours
% trip 4 (Nov. 25)
%    - sampler 1: starting 13:00 PM sampling every 7 hours
%    - sampler 2: starting 14:45 PM sampling every 7 hours
%    - sampler 3: starting 16:30 PM sampling every 7 hours
%    - sampler 4: starting 18:15 PM sampling every 7 hours

nTrip = 4;
nSampler = 4; % number of deployed samplers
nBottle = 24; % number of bottles in each sampler
dates = cell(nTrip,nSampler,nBottle);

for iTrip = 1:nTrip
    switch iTrip
        case 1
            dt_sampler = minutes(30); 
            dt_bottle = hours(2); 
            t0 = datetime('2024-11-18T18:00','InputFormat','uuuu-MM-dd''T''HH:mm');            
        case 2
            dt_sampler = minutes(30); 
            dt_bottle = hours(2); 
            t0 = datetime('2024-11-20T16:00','InputFormat','uuuu-MM-dd''T''HH:mm');            
        case 3
            dt_sampler = minutes(45); 
            dt_bottle = hours(3); 
            t0 = datetime('2024-11-22T13:00','InputFormat','uuuu-MM-dd''T''HH:mm');            
        case 4
            dt_sampler = minutes(105); 
            dt_bottle = hours(7); 
            t0 = datetime('2024-11-25T13:00','InputFormat','uuuu-MM-dd''T''HH:mm');            
    end
    
    
    for iSampler = 1:nSampler
        for iBottle = 1:nBottle
            dates{iTrip,iSampler,iBottle} = t0 + (iSampler-1)*dt_sampler +(iBottle-1)*dt_bottle;
        end
    end
end

%% load previous data
flName = [resultsDir 'fieldSamples.mat'];
if exist(flName,'file')
    load([resultsDir 'fieldSamples.mat'],'sampleDate','area')
    sampleDate_pre = sampleDate; 
    area_pre = area;
    clear area sampleDate
end

%% file name & time recognition
flNameStruct = dir(dataDir);
nSample = 0;
for il = 3:numel(flNameStruct)
    if strcmp(flNameStruct(il).name(end-2:end),'txt')
        if strcmp(flNameStruct(il).name(end-18:end),'_CorrectionData.txt')
            nSample = nSample + 1;
            flNameList{nSample,1} = flNameStruct(il).name;
        end
    end
end


sampleDate = datetime;
for iSample = 1:length(flNameList)
    flName = flNameList{iSample};
    note = textscan(flName, '%s','delimiter','_' ); note = note{1};

    iTrip = str2num(note{1});
    iSampler = str2num(note{2});
    iBottle = str2num(note{3});

    sampleDate(iSample,1) = dates{iTrip,iSampler,iBottle};
end

flNameList_ = cell(size(flNameList));
[sampleDate,ind] = sort(sampleDate);
for iFl = 1:numel(flNameList)
    flNameList_{iFl} = flNameList{ind(iFl)};
end
flNameList = flNameList_;
%% Peak Fitting
area = cell(size(sampleDate));

indNow = [1:length(flNameList)]';
if exist('area_pre','var')
    [~,indPre] = intersect(sampleDate,sampleDate_pre);
    indExist = 0;
    for i = indPre'
        indExist = indExist + 1;
        area{i} = area_pre{indExist};
    end
    
    [~,indNow] = setdiff(sampleDate,sampleDate_pre);
    
end


for iSample = indNow'
    fprintf('%d/%d samples\n',iSample,length(flNameList))
    nIter=20;
    area{iSample} = peakFitter(dataDir,figDir,flNameList{iSample},'water',nIter);
end


save([resultsDir 'fieldSamples.mat'],'sampleDate','area')


%% Breakthrough Curves (direct sampling)
clc; clf

dyeType = {'uranine'};

stdFlName = sprintf('%sSTD_%s.mat',resultsDir,dyeType{1});
load([resultsDir 'fieldSamples.mat'],'sampleDate','area')
load(stdFlName,'slope'); 

area = cell2mat(area);
    
if strcmp(dyeType{1}, 'uranine')
    area = area(:,1);
elseif strcmp(dyeType{1}, 'RWT')
    area = area(:,2);
elseif strcmp(dyeType{1}, 'SrB')
    area = area(:,3);
end

% area = area-area(1); % deduct background fluorescence

cnc = area.*slope;

travelTime = (sampleDate-t_uranineInj)/hours(24);



plot(sampleDate,cnc)


xlabel('dates')
ylabel('concentration [ppb]')



subplot(121); 
plot(sampleDate,cnc,'k-','linewidth',2); hold on

ax = gca;
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

xlabel('dates [day]')
ylabel('concentration [ppb]')
ylim([-0.1 4])


subplot(122); 
plot(log10(travelTime),log10(cnc),'k-','linewidth',2);  hold on


ytick = [0.001 0.01 0.1  1 10 100 ]';
xtick = [0.2 0.5 1 2 4 8 16]';

ax = gca;
ax.XTick = log10(xtick);
ax.YTick = log10(ytick);
ax.XTickLabel = num2str(xtick);
ax.YTickLabel = num2str(ytick);
ax.FontSize = 14;
xlim([log10(0.1) log10(10)])
ylim([log10(0.05) log10(10)])

xlabel('travel time [day]')
ylabel('concentration [ppb]')

