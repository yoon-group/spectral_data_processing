clear all; clc; clf
%% ========================================================================
% 
%                   Standard samples processing
% 
%                                                      Yoon Research Group
% -------------------------------------------------------------------------
%     Dilution factors (sample : Square lake water)
%        20x: 0.2 ml : 4 ml  -> 21 x
%        10x: 0.4 ml : 4 ml  -> 11 x
%         5x: 0.8 ml : 4 ml  ->  6 x
%         2x:   2 ml : 2 ml  ->  2 x
% =========================================================================

campaign = 'SinkRise2024'; % 'SinkRise2024' or 'BearSpring2023'

archiveDir = sprintf('%s/dyeTracingData/',pwd);

resultsDir = sprintf('%s/processed/%s/',pwd,campaign);
if ~exist(resultsDir,'dir'); mkdir(resultsDir); end

figDir = [resultsDir 'fig/STDs/']; 
if ~exist(figDir,'dir'); mkdir(figDir); end

%% file name recognition
sampleType = 'water'; % 'water' or 'eluent'

for dyeType = {'uranine'} % 'SrB', 'uranine', or 'RWT'
    dataDir = sprintf('%s%s/STD_%s/',archiveDir,campaign, dyeType{1} );
    flNameStruct = dir(dataDir);
    nSample = 0;
    for il = 3:numel(flNameStruct)
        if strcmp(flNameStruct(il).name(end-2:end),'txt')
            nSample = nSample + 1;
            flNameList{nSample,1} = flNameStruct(il).name;
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
            tmp = textscan(note{4},'%s','delimiter','x'); tmp = tmp{1};
            dilutionFactor(iSample,1) = str2num(tmp{1});
        else
            dilutionFactor(iSample,1) = 1;
        end
    end
    dilutionFactor(find(dilutionFactor == 20)) = 21;
    dilutionFactor(find(dilutionFactor == 10)) = 11;
    dilutionFactor(find(dilutionFactor == 5)) = 6;
    
    %% Peak Fitting
    area = cell(length(flNameList_todo),1);
    for iSample = 1:length(flNameList_todo)
        nIter=10; if cnc(iSample)>1; nIter=50; end
        area{iSample} = peakFitter(dataDir,figDir,flNameList_todo{iSample},sampleType,nIter);
        fprintf('sample %d/%d \n',iSample,length(flNameList_todo))
    end
    
    %% standard data analysis
    area = cell2mat(area);
    
    if strcmp(dyeType{1}, 'uranine')
        area = area(:,1);
    elseif strcmp(dyeType{1}, 'RWT')
        area = area(:,2);
    elseif strcmp(dyeType{1}, 'SrB')
        area = area(:,3);
    end
    
    area = dilutionFactor.*area;
        
    % costF = @(x) sqrt(sum((log10(area) + log10(x) - log10(cnc)).^2));
    costF = @(x) sqrt(sum((area .* x - cnc).^2));
    options = optimset('PlotFcns',@optimplotfval,'MaxIter',1000);
    slope = fminsearch(costF,1e-3,options); 

    SSR = sum((area*slope - cnc).^2); % sum of squared error
    R2 = 1- SSR / sum((cnc - mean(cnc)).^2); % R^2 value

    flNameSTD = sprintf('%sSTD_%s.mat',resultsDir,dyeType{1} );
    save(flNameSTD,'area','cnc','slope')
    
    %%
    clf; 
    subplot(121)
    plot(area,slope*area,'-'); hold on
    scatter(area,cnc,'o')

    xl = xlim; yl = ylim;
    text(xl(1)+diff(xl)/10,yl(2)-diff(yl)/10,sprintf('y = (%.4e) * x',slope))
    text(xl(1)+diff(xl)/10,yl(2)-diff(yl)/10*2,sprintf('R^2 = %.4f',R2))

    ylabel('concentration [ppb]');
    xlabel('spectral peak area')

    subplot(122)
    loglog(area,slope*area,'-'); hold on
    scatter(area,cnc,'o')
    
    ylabel('concentration [ppb]');
    xlabel('spectral peak area')

    title(sprintf('log-scale'))

    fig = gcf;
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 8.5 11];

    flNameFig = sprintf('%s%s_calibration.png',figDir, dyeType{1});

    print(fig,flNameFig,'-dpng')


end