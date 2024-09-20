function [area] = peakFitter(dataDir,figDir,flName,dyeType,sampleType)    

%% ========================================================================
%                          peakFitting function
%                                                           Seonkyoo Yoon
% -------------------------------------------------------------------------
% This code is to automate the multi-peak fitting process needed for 
% anlayzing the curves produced by a spectrofluorophotometer. The code 
% estimates the peak area correpsponding to a dye tracer. 
% 
% Currently, this code is designed for fluorescein dye which typically has 
% the peak center around 508 nm for water and 514 nm for eluant. For 
% different dye types, the contraint varialbes like A, B, Aeq, Beq, lb, and 
% ub need to be adjusted. To understand the use of the variables, take a 
% look at the MATLAB manual page for the function 'fmincon' which is the 
% optimizer used in this code for the multi-peak fitting.   
% 
% dyeType = currently only 'fluorescein' is avaialbe. 
% sampleType =  either 'water' or 'eluant'
% 
% =========================================================================
    %% data laod
    fidDataFl = fopen([dataDir flName],'r');
    g = textscan(fidDataFl, '%s','delimiter','\n'); fclose(fidDataFl);
    g = g{1};
    
    
    dataStart = 1;
    for iL = 1:length(g)
        if contains(g{iL}, 'Wavelength nm') 
            dataStart = iL + 1;
            break;
        end
    end
    
    iData = 0;
    for iL = dataStart:length(g)
        iData = iData + 1;
        note = textscan(g{iL}, '%s','delimiter',',' ); note = note{1};
        wavelength(iData,1) = str2num(note{1});
        intensity(iData,1) = str2num(note{2});
    end
    
    
    %% contraints according to dye and samples
    if strcmp(dyeType,'NOM') 
        clc
        % nPeak = 9;
        params = [360      390      440                500; ... center
                   300     10*ones(1,2)        10*ones(1,1); ... height
                   100     30*ones(1,2)        30*ones(1,1);         ... hwhm
                   0.8    2*ones(1,2)    30*ones(1,1)];  % shape
        nPeak = size(params,2);
        params = params(:);
        
        costFnc = @(params) sqrt(mean((intensity-peakSum(wavelength,params)).^2)); 
        
        A = [];
        b = [];
        
        fixed = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 ones(1,nPeak-1) ];    % shape
        fixed = fixed(:);
        
        shape = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 2 2   5 ];    % shape
        shape = shape(:);
        
        Aeq = diag(fixed);
        beq = shape;
        
        lb = [320   370   380   450   ; ... center
              5*ones(1,nPeak);    ... height
                0     0     0   0   ;    ... hwhm
              0.4 zeros(1,nPeak-1);];    % shape
        ub = [370   450   430   570   ; ... center
              700*ones(1,nPeak); ... height
              300   100   100   100  ;    ... hwhm
              0.9 100*ones(1,nPeak-1);]; % shape
    elseif strcmp(dyeType,'tracer') 
        clc
        nPeak = 9;
        params = [360 linspace(400, 500, 3) 508 582 linspace(600,740,3); ... center
                   300    10*ones(1,3)      350 380  10*ones(1,3); ... height
                   100 30*ones(1,3) 10 10 30*ones(1,3);         ... hwhm
                   0.8    2*ones(1,3)    10  10  30*ones(1,3)];  % shape
        params = params(:);
        
        costFnc = @(params) sqrt(mean((intensity-peakSum(wavelength,params)).^2)); 
        
        A = [];
        b = [];
        
        fixed = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 ones(1,nPeak-1) ];    % shape
        fixed = fixed(:);
        
        shape = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 2 2 2 10 10 10 5 5];    % shape
        shape = shape(:);
        
        Aeq = diag(fixed);
        beq = shape;
        
        lb = [320   370   380   420   508  580  510   550   550; ... center
              zeros(1,nPeak);    ... height
                0     0     0     0     5     5     0     0     0;    ... hwhm
              0.4 zeros(1,nPeak-1);];    % shape
        ub = [370   450   420   480   509  584  540   780   780; ... center
              700*ones(1,nPeak); ... height
              300   100   100   100    15   15  100   100   100;    ... hwhm
              0.9 100*ones(1,nPeak-1);]; % shape
    elseif strcmp(dyeType,'uranine') && strcmp(sampleType,'water')
        clc
        nPeak = 8;
        params = [360 linspace(400, 500, 3) 508 linspace(600,740,3); ... center
                   300    10*ones(1,3)   350   10*ones(1,3); ... height
                   100 30*ones(1,3) 10 30*ones(1,3);         ... hwhm
                   0.8    2*ones(1,3)    10    30*ones(1,3)];  % shape
        params = params(:);
        
        costFnc = @(params) sqrt(mean((intensity-peakSum(wavelength,params)).^2)); 
        
        A = [];
        b = [];
        
        fixed = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 ones(1,7) ];    % shape
        fixed = fixed(:);
        
        shape = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 2 2 2 10 10 5 5];    % shape
        shape = shape(:);
        
        Aeq = diag(fixed);
        beq = shape;
        
        lb = [320   370   380   420   508   550   550   550; ... center
              zeros(1,nPeak);    ... height
                0     0     0     0     5     0     0     0;    ... hwhm
              0.4 zeros(1,nPeak-1);];    % shape
        ub = [370   450   420   450   509   600   780   780; ... center
              inf*ones(1,nPeak); ... height
              300   100   100   100    15   100   100   100;    ... hwhm
              0.9 100*ones(1,nPeak-1);]; % shape
    
    elseif strcmp(dyeType,'uranine') && strcmp(sampleType,'eluant')
        clc
        nPeak = 8;
        params = [360 linspace(400, 500, 3) 514 linspace(600,740,3); ... center
                   300    10*ones(1,3)   350   10*ones(1,3); ... height
                   100 30*ones(1,3) 10 30*ones(1,3);         ... hwhm
                   0.8    2*ones(1,3)    10    30*ones(1,3)];  % shape
        params = params(:);
        
        costFnc = @(params) sqrt(mean((intensity-peakSum(wavelength,params)).^2)); 
        
        A = [];
        b = [];
        
        fixed = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 ones(1,7) ];    % shape
        fixed = fixed(:);
        
        shape = [zeros(1,nPeak); ... center
                 zeros(1,nPeak); ... height
                 zeros(1,nPeak); ... hwhm
                 0 2 2 2 10 10 5 5];    % shape
        shape = shape(:);
        
        Aeq = diag(fixed);
        beq = shape;
        
        lb = [320   370   380   420   512   514   550   550; ... center
              zeros(1,nPeak);    ... height
                0     0     0     0     5     0     0     0;    ... hwhm
              0.4 zeros(1,nPeak-1);];    % shape
        ub = [370   450   420   480   516   550   780   780; ... center
              inf*ones(1,nPeak); ... height
              300   100   100   100    15   100   100   100;    ... hwhm
              0.9 100*ones(1,nPeak-1);]; % shape
    
    end
    %% parameter fitting using fmincon
    options = optimoptions('fmincon','Display','iter');%'MaxFunctionEvaluations',1e+10,'MaxIterations',1e+5);
        
    for iter =1:10
        params = fmincon(costFnc,params,A,b,Aeq,beq,lb,ub);
    end
    
    %% visualization
    clf; 
    paramsPrint = reshape(params,4,nPeak)';
    [curve] = peakSum(wavelength,params);


    % clr = 'kkkkrbkkk';
    clr = 'kkkkkkkkk';

    subplot(100,1,1:40)
    plot(wavelength,intensity); hold on
    plot(wavelength,curve)

    for iPeak = 0:nPeak-1
        if 1%params(4*iPeak+2) > 1
            peak = Pearson7(wavelength,params(4*iPeak+1:4*iPeak+4));
            plot(wavelength,peak,[clr(iPeak+1) ':']); hold on
        end
    end

    title(flName,'Interpreter','none')

    subplot(100,1,46:55)
    plot(wavelength,intensity-curve); hold on
    plot(wavelength,zeros(size(intensity)))

    yl = ylim;
    xl = xlim;

    for iPeak = 0:nPeak-1
        if 1%params(4*iPeak+2) > 1
            plot([params(4*iPeak+1) params(4*iPeak+1)],[-15 15],clr(iPeak+1));  
        end
    end
    xlim(xl)
    ylim([-7 7])


    SSR = sum((curve - intensity).^2); % sum of squared error
    R2 = 1- SSR / sum((intensity - mean(intensity)).^2);
    weigth = 1./max(sqrt(intensity), ones(size(intensity)));
    WSSR = sum(weigth.^2.*(curve - intensity).^2); % weighted sum of squared error

    % degree of freedom = data dimensiton - #. parameter + #. fixed parameter
    DoF = length(intensity) - 4*nPeak + sum(fixed); 

    subplot(100,1,61:100)
    yl = ylim;

    txt = sprintf('1 sigma = %.2f   [sqrt(WSSR/DoF)] ', sqrt(WSSR/DoF));
    text(0,yl(2),txt,'fontsize',12)

    txt = sprintf('WSSR=%.4f   DoF=%d   WSSR/DoF=%.4f   SSR=%.4f   R2=%.4f', ...
        WSSR,DoF,WSSR/DoF,SSR,R2);
    text(0,yl(2)-diff(yl)/(nPeak+1),txt,'fontsize',12)


    txt = sprintf('%10s %12s %12s %12s %12s %12s \n',...
                    'peak', 'center','height','hwhm','shape','area');

    yl = ylim;
    text(0,yl(2)-2*diff(yl)/(nPeak+1),txt,'fontsize',12)
    area = [];
    for iPeak = 0:nPeak-1
        if 1%params(4*iPeak+2)>1
        if iPeak == 4
            peak = Pearson7(wavelength,params(4*iPeak+1:4*iPeak+4));
            area = [area trapz(wavelength,peak)];

            txt = sprintf('%10d %12.3f %12.3f %12.3f %12.3f %12.3f\n', ...
            iPeak+1, paramsPrint(iPeak+1,:), area(end));
        else
            txt = sprintf('%10d %12.3f %12.3f %12.3f %12.3f \n', ...
            iPeak+1, paramsPrint(iPeak+1,:));
        end
        text(0,yl(2)-(iPeak+3)*diff(yl)/(nPeak+1),txt,'fontsize',12)
        end
    end

    axis off

    drawnow;

    %
    fig = gcf;
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 8.5 11];

    flNameFig = sprintf([figDir flName(1:end-4) '.png']);

    print(fig,flNameFig,'-dpng')

%%
% [curve] = peakSum(wavelength,params);
%     clf
%     % subplot(100,1,1:40)
%     plot(wavelength,curve,'k','linewidth',2); hold on
%     scatter(wavelength,intensity,10,'filled','k'); hold on
% 
% 
%     for iPeak = 0:nPeak-1
%         if 1%params(4*iPeak+2) > 1
%             peak = Pearson7(wavelength,params(4*iPeak+1:4*iPeak+4));
%             plot(wavelength,peak,'k-','linewidth',1); hold on
%         end
%     end
% 
%     % title(flName,'Interpreter','none')
%     xlim([390 600])
%     ylim([-1 120])
% 
%     xlabel('Excitation Wavelength [nm]')
%     ylabel('Synchronous Scan Intensity')
% 
% 
%     set(gca,'fontsize',20)
% 
% 
% 
%     drawnow;
% 
% 
% 
% 
% 
%     fig = gcf;
%     fig.Units = 'inches';
%     fig.PaperPosition = [0 0 9 6];
% 
% 
% 
%     flNameFig = sprintf([figDir flName(1:end-4) '.png']);
% 
%     print(fig,flNameFig,'-dpng')

end


function [curve] = peakSum(x,params)    
    nPeak = length(params)/4;

    curve = Pearson7(x,params(1:4));
    for iPeak = 1:nPeak-1
        curve = curve + Pearson7(x,params(4*iPeak+1:4*iPeak+4));
    end
    
end


function [curve] = Pearson7(x,params)    
    curve = params(2)./...
        (1+((x-params(1))./params(3)).^2*(2^(1/params(4))-1)).^params(4);
end

