clear all; clc; clf
%% ========================================================================
% 
%                   Standard fluorescein solution samples
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
%     Dilution factors (sample : Square lake water)
%        20x: 0.2 ml : 4 ml  -> 21 x
%        10x: 0.4 ml : 4 ml  -> 11 x
%         5x: 0.8 ml : 4 ml  ->  6 x
%         2x:   2 ml : 2 ml  ->  2 x
% 
%     concentrations of standard samples
%        cnc = [0.00995; 0.0305; 0.1044; 0.3152; 0.97891; 2.87]; % [ppb]
% =========================================================================

codeDir = [pwd '/'];
dataDir = [codeDir 'uranineSTDs/'];
resultsDir = [dataDir 'results/'];
figDir = [dataDir 'fig/']; if ~exist(figDir,'dir'); mkdir(figDir); end

dyeType = 'uranine';
sampleType = 'water';
%% file name recognition
flNameStruct = dir(dataDir);
nSample = 0;
for il = 3:numel(flNameStruct)
    if strcmp(flNameStruct(il).name(end-2:end),'txt')
        % if strcmp(flNameStruct(il).name(end-18:end),'_CorrectionData.txt')
            nSample = nSample + 1;
            flNameList{nSample,1} = flNameStruct(il).name;
        % end
    end
end
flNameList_todo = flNameList;

%% standard solutions
cnc = nan(nSample,1);
for iSample = 1:nSample
    flName = flNameList{iSample}(1:end-4);
    note = textscan(flName, '%s','delimiter','_' ); note = note{1};

    cnc(iSample) = str2double(sprintf('%s.%s',note{2},note{3}));
    if length(note) >= 4
        tmp = textscan(note{5},'%s','delimiter','x'); tmp = tmp{1};
        dilutionFactor(iSample,1) = str2num(tmp{1});
    else
        dilutionFactor(iSample,1) = 1;
    end
end
dilutionFactor(find(dilutionFactor == 20)) = 21;
dilutionFactor(find(dilutionFactor == 10)) = 11;
dilutionFactor(find(dilutionFactor == 5)) = 6;

%% Peak Fitting
area = zeros(length(flNameList_todo),1);
for iSample = 1:length(flNameList_todo)
    area(iSample) = peakFitter(dataDir,figDir,flNameList_todo{iSample},dyeType,sampleType);
end

%% standard data analysis
area = dilutionFactor.*area;


% costF = @(x) sqrt(sum((log10(area) + log10(x) - log10(cnc)).^2));
costF = @(x) sqrt(sum((area .* x - cnc).^2));

options = optimset('PlotFcns',@optimplotfval,'MaxIter',1000);
slope = fminsearch(costF,1e-3,options); 

%%
close all
clf; %scatter(area,cnc); hold on
subplot(121)
plot(area,slope*area,'-'); hold on
scatter(area,cnc,'o')
subplot(122)
loglog(area,slope*area,'-'); hold on
scatter(area,cnc,'o')

title(sprintf('slope: %.3e',slope))
save([resultsDir 'STD.mat'],'slope')