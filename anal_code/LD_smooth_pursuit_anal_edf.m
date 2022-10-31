%This script was developped to analyze eye movement data from the smooth
%pursuit experiment, the data was acquired in the eye tracking room 425
%using the task called LD_phantom_oscill_V3.m
%10/26/2022

folderDir = '/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/data/smooth_pursuit_10252022/';
folders = dir(folderDir);
folders = {folders(~contains({folders.name},'.') & contains({folders.name},'gamma')).name};
%compute average response curves of each subject
% meanx = nan(1321,2, length(folders));
for i =1:length(folders)
    
    fileNames = dir(strcat(folderDir, folders{i}));
    fileNames = string({fileNames(contains({fileNames.name},'smooth') | contains({fileNames.name},'demo')).name}); %
    eyeDat = load(strcat(folderDir, folders{i},'/',fileNames(contains(fileNames,'eyeDat'))));
    %edfName = char(strcat(folderDir, folders{i},'/',fileNames(contains(fileNames,'demo'))));
    %edfDat = Edf2Mat(edfName);
    exptName = strcat(folderDir, folders{i},'/',fileNames(~contains(fileNames,'eyeDat') & ~contains(fileNames,'demo')));
    exptDat = load(exptName);
    
    xbsl = exptDat.exp.resolution.width/2;
    ybsl = exptDat.exp.resolution.height/2;
    xDat = eyeDat.EyeData.mx-xbsl; %center 0 relative to middle of the screen
    yDat = eyeDat.EyeData.my-ybsl; %center 0 relative to middle of the screen
    pixperdeg = exptDat.exp.ppd;
    x = xDat/pixperdeg; %convert from pixels to degrees
    y = yDat/pixperdeg; %convert from pixels to degrees
    
    xpTr = nan(1300,10, length(exptDat.exp.conds));
    ypTr = nan(1300,10, length(exptDat.exp.conds));
    origX = nan(1300,10, length(exptDat.exp.conds));
    origY = nan(1300,10, length(exptDat.exp.conds));
    
    thresh = 1; %eye position quantile detection threshold
    for c =1:length(exptDat.exp.conds)
        blockdiff = unique([find(diff(exptDat.exp.longFormConds == c)) length(exptDat.exp.longFormConds)]); %get block onset time and offset time
        for n =1:length(blockdiff)/2
            xtr = x(blockdiff(2*(n-1)+1):blockdiff(2*n));
            ytr = y(blockdiff(2*(n-1)+1):blockdiff(2*n));
            
            origX(1:length(xtr),n,c) = xtr;
            origY(1:length(ytr),n,c) = ytr;
            
            xthr = abs(xtr)<quantile(abs(xtr(213:end-212)),thresh);
            ythr = abs(ytr)<quantile(abs(ytr(213:end-212)),thresh);
            
            xp = xtr.*xthr.*ythr - nanmean(xtr.*xthr.*ythr);
            yp = ytr.*xthr.*ythr - nanmean(ytr.*xthr.*ythr);
  
            normX = xp/max(abs(xp),[],'all');
            normY = yp/max(abs(yp),[],'all');
            xpTr(1:length(normX),n,c) = normX;
            ypTr(1:length(normY),n,c) = normY;
        end
        
    end
    
    subjDat.(folders{i}).normx = xpTr;
%     for c =1:length(exptDat.exp.conds)
%         meanx(:,c,i) = mean(subjDat.(folders{i}).normx(:,:,c),2); 
%     end
    subjDat.(folders{i}).normy = ypTr;
    
    subjDat.(folders{i}).origx = origX;
    subjDat.(folders{i}).origy = origY;
    
    %subjDat.(folders{i}).edf = edfDat;
    phases = exptDat.exp.stim.oscillation;%cos(2*pi*cycPerSec*stimDur/(stimDur*flipsPerSec)-2*pi*flipTimes./stimDur).*360;
    
    conditions = {'Minbg', 'Meanbg'};
    xax = (1:size(xpTr,1))./60;
    for n =1:size(xpTr,3)
        figure('Position',[100 100 1200 600]); %
        
        plot(xax,xpTr(:,:,n))
        hold on
        plot(xax(1:length(repmat(phases,1,6))),repmat(phases,1,6),'LineWidth', 2,'col','k')
        set(gca,'box','off')
        set(gca, 'LineWidth',2)
        sgtitle(strcat(conditions{n}, sprintf(' %s ',folders{i})),'fontweight','bold','Interpreter', 'none')
        
        xlabel('time from stimulus onset (s)','fontsize',18,'fontweight','bold')
        ylabel('x position (normalized)','fontsize',18,'fontweight','bold')
        
        yticklab = get(gca,'YTickLabel');
        set(gca,'YTickLabel',yticklab,'fontsize',12)
        xticklab = get(gca,'XTickLabel');
        set(gca,'XTickLabel',xticklab,'fontsize',12)
        
        legend('gaze', 'stimulus')%,'location','bestoutside')
    end
    
%         plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/anal_plots/');
%          saveas(gcf,strcat(plotdir, sprintf('%s_%s_100percthr_electrooculogram_10252022.png',folders{i},conditions{n})));
%          saveas(gcf,strcat(plotdir, sprintf('%s_%s_100percthr_electrooculogram_10252022.svg',folders{i},conditions{n})));
end


%% detect saccades and exclude them

%detect microsaccades
samplerate = 60;
for i =1:length(folders)
    
    eyeMovDat.(folders{i}) = eyeMovDetect(char(folderDir),char(folders{i}),subjDat.(folders{i}).origx, subjDat.(folders{i}).origy, samplerate, exptDat.exp.conds);
end


%% fit cosine model to the data

fitCoefs = struct();
amp = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));
freq = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));
phase = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));
offset = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));

for i =1:length(folders)
    for c = 1:length(exptDat.exp.conds)
        a_sol = [];
        b_sol = [];
        c_sol = [];
        d_sol = [];
        for t =1:size(subjDat.(folders{i}).normx,2)
            y = subjDat.(folders{i}).normx(1:end,t,c);
            yfit = y(213:end);
            yfit = yfit(~isnan(y(213:end)));
            x = 0:1/60:exptDat.exp.blockLength;
            x = x(213:length(y));
            xfit = x(~isnan(y(213:end)));
            data.y = yfit';
            data.x = xfit;
            x0 = [1, pi/2,0, 0];
            L = @(params)CosMse(data, params);
            %L = @(params)CosNorm(data, params);
            [xsol,fval,exitflag,output] = fminsearch(@(params) L(params), x0);
            
            a_sol = [a_sol xsol(1)];
            b_sol = [b_sol xsol(2)];
            c_sol = [c_sol xsol(3)];
            d_sol = [d_sol xsol(4)];
            amp(t,i,c) = xsol(1);
            freq(t,i,c) = xsol(2);
            phase(t,i,c) = xsol(3);
            offset(t,i,c) = xsol(4);
        end
        fitCoefs.(folders{i}).(conditions{c}).amplitude = a_sol;
        fitCoefs.(folders{i}).(conditions{c}).frequency = b_sol./(2*pi);
        fitCoefs.(folders{i}).(conditions{c}).phase = c_sol;
        fitCoefs.(folders{i}).(conditions{c}).offset = d_sol;

    end
end

%plot and compare to stimulus a and b
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmap = flip(cmaps(2).map) ;

%p = figure('Position',[100 100 900 500]);
p = figure('Position',[100 100 1200 600]); %
[h,pos]=tight_subplot(1,4,[.02 .05],[.1 .15],[.05 .02]); %use tight plot for adjusting spacing between axes, and vertical and horizontal margin size

fit = {'Amplitude (normalized)', 'Frequency (cyc/s)', 'Phase (rad)', 'Offset (normalized)'};
fitNames = {'amp', 'freq', 'phase', 'offset'};
for f =1:length(fit)
    %plot the data points
    % and the mean ±95%CI
    eval(strcat('yvar =' , sprintf('%s',fitNames{f}),';'));
    mYvar = nanmean(yvar,[1,2]);
    
    %95%CI
    ci_high = mean(yvar,[1,2]) + std(yvar,[],[1,2]);%/sqrt(size(yvar,1)*size(yvar,2));
    ci_low = mean(yvar,[1,2]) - std(yvar,[],[1,2]);%/sqrt(size(yvar,1));
    yvarlab = {sprintf('%s',fit{f})};
    
    axes(h(f));
    for c =1:length(exptDat.exp.conds)
        xplot = c*ones(1, length(yvar(:,1,c))*length(yvar(1,:,c)));
        yplot =  reshape(yvar(:,:,c), 1,size(yvar(:,:,c),1)*size(yvar(:,:,c),2));
        scatter(xplot,yplot,20,'MarkerFaceColor',cmap(4,:), 'MarkerEdgeColor',cmap(4,:),'LineWidth',0.5);
        hold on;
        scatter(c*1,mYvar(c),60,'MarkerFaceColor',cmaps(1).map(4,:), 'MarkerEdgeColor',cmaps(1).map(4,:),'LineWidth',1.5); %mean
        hold on
        line([c*1 c*1], [ci_low(c) ci_high(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI vertical
        hold on
        line([c*1-0.1 c*1+0.1], [ci_low(c) ci_low(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI whiskers
        hold on
        line([c*1-0.1 c*1+0.1], [ci_high(c) ci_high(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI whiskers
        hold on

    ylim([min(yvar(:,:,c),[],'all')-0.1, max(yvar(:,:,c),[],'all')+0.1]);
    end
    
    % Set up axes.
    xlim([0, 3]);
    curtick = get(gca, 'xTick');
    xticks(unique(round(curtick)));
    ax = gca;
    ax.XTick = [1:length(conditions)];
    ax.XTickLabels = conditions;
    ylab = yvarlab{:};
    ylabel(ylab,'fontweight','bold','fontsize',16)
    yticklab = get(gca,'YTickLabel');
    set(gca,'YTickLabel',yticklab,'fontsize',12)
  
end
legend('Trial Model fit','Mean±SD')
sgtitle('fminsearch mse results', 'FontWeight', 'Bold','fontsize',22)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/anal_plots/');
saveas(gcf,strcat(plotdir, 'model_fminsearch_params_10272022.png'));
   
%% EDF file analysis

folderDir = '/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/data/smooth_pursuit_10272022/';
folders = dir(folderDir);
folders = {folders(~contains({folders.name},'.')).name};
%compute average response curves of each subject
for i =1:length(folders)
    
    fileNames = dir(strcat(folderDir, folders{i}));
    fileNames = string({fileNames(contains({fileNames.name},'smooth') | contains({fileNames.name},'demo')).name}); %
    edfName = char(strcat(folderDir, folders{i},'/',fileNames(contains(fileNames,'demo'))));
    edfDat = Edf2Mat(edfName);
    exptName = strcat(folderDir, folders{i},'/',fileNames(~contains(fileNames,'eyeDat') & ~contains(fileNames,'demo')));
    exptDat = load(exptName);
    
    xbsl = exptDat.exp.resolution.width/2;
    ybsl = exptDat.exp.resolution.height/2;
    xDat = edfDat.RawEdf.FSAMPLE.gx(2,:)-xbsl; %center 0 relative to middle of the screen
    yDat = edfDat.RawEdf.FSAMPLE.gy(2,:)-ybsl; %center 0 relative to middle of the screen
    pixperdeg = exptDat.exp.ppd;
    x = xDat/pixperdeg; %convert from pixels to degrees
    y = yDat/pixperdeg; %convert from pixels to degrees
    
    [eyePos,blinks] = removeBlinks(x,y); %detect blinks and replace them by NaNs 
    x = eyePos(1,:);
    y = eyePos(2,:);
    x = fillmissing(x,'previous'); %replace nans by previous value in the array (needed to detect saccades)
    y = fillmissing(y,'previous'); %replace nans by previous value in the array
    
%     The commented out code detects saccades but can't really do anything
%     about it
%     figure()
%     subplot(9,1,1)
%     plot(x(354950:355500))
%     legend('xpos')
%     %detect saccades
%     gpos = sqrt(x.^2+y.^2);
%     subplot(9,1,2)
%     plot(gpos(354950:355500))
%     legend('gpos')
%     vel = diff(gpos)*1000; % gaze velocity
%     subplot(9,1,3)
%     plot(vel(354950:355500))
%     legend('vel')
%     FIRfilt = fir1(80,30/500, 'low'); % filter order = 48 lowpass filter at 30 Hz, alway normalize by nyquist frequency = 1000Hz/2
%     filtVel = filter(FIRfilt,1,vel); %filter velocity
%     subplot(9,1,4)
%     plot(filtVel(354950:355500))
%     legend('filtvel')
%     movVel = movmean(filtVel, 40); %running average of 40 ms window
%     subplot(9,1,5)    
%     plot(movVel(354950:355500))
%     legend('run avg vel')
%     acc = diff(vel); % gaze acceleration
%     subplot(9,1,6)
%     plot(acc(354950:355500))
%     legend('acc')
%     filtAcc = filter(FIRfilt,1,acc); %filter acceleration
%     subplot(9,1,7)
%     plot(filtAcc(354950:355500))
%     legend('filt acc')
%     movAcc = movmean(filtAcc, 40); %running average of 40 ms window
%     
%     subplot(9,1,8)
%     plot(movAcc(354950:355500))
%     legend('run avg acc')
%     excludeSacc= find(movAcc >= 1000);
%     
%     x(excludeSacc) = nan;
%     y(excludeSacc) = nan;
%     
%     x = fillmissing(x,'previous');
%     y = fillmissing(y,'previous');
%     
%     xq = interp1(1:length(x),x, excludeSacc,'linear');
%     yq = interp1(1:length(y),y, excludeSacc,'linear');
%     
%     x(excludeSacc) = xq;
%     y(excludeSacc) = yq;
%     subplot(9,1,9)
%     plot(x(354950:355500))
%     legend('clean x')
  
    xpTr = nan(22000,10, length(exptDat.exp.conds));
    ypTr = nan(22000,10, length(exptDat.exp.conds));
    origX = nan(22000,10, length(exptDat.exp.conds));
    origY = nan(22000,10, length(exptDat.exp.conds));
    blinkTrs = nan(22000,10, length(exptDat.exp.conds));
    %thresh = 1; %eye position quantile detection threshold
    for c =1:length(exptDat.exp.conds)
        conds = exptDat.exp.condShuffle == c;
        onsetT = edfDat.Events.Messages.time(strcmp(edfDat.Events.Messages.info, 'STIM_ONSET')); %stimulus onset times
        offsetT = edfDat.Events.Messages.time(strcmp(edfDat.Events.Messages.info, 'STIM_OFFSET')); %stimulus offset times
        onsetTcond = onsetT(conds);
        offsetTcond = offsetT(conds);
        for n =1:length(onsetTcond)
            xtr = x(find(edfDat.RawEdf.FSAMPLE.time == onsetTcond(n)):find(edfDat.RawEdf.FSAMPLE.time == offsetTcond(n)));
            ytr = y(find(edfDat.RawEdf.FSAMPLE.time == onsetTcond(n)):find(edfDat.RawEdf.FSAMPLE.time == offsetTcond(n)));
            blinktr = blinks(find(edfDat.RawEdf.FSAMPLE.time == onsetTcond(n)):find(edfDat.RawEdf.FSAMPLE.time == offsetTcond(n)));
            
            origX(1:length(xtr),n,c) = xtr;
            origY(1:length(ytr),n,c) = ytr;
    
%             xp = xtr - nanmean(xtr);
%             yp = ytr - nanmean(ytr);
%   
%             normX = xp/max(abs(xp),[],'all');
%             normY = yp/max(abs(yp),[],'all');
%             xpTr(1:length(normX),n,c) = normX;
%             ypTr(1:length(normY),n,c) = normY;
            blinkTrs(1:length(blinktr),n,c) = blinktr;
        end
        
    end
    samplerate = 1000;
    [eyeMovDat,cOrigX,cOrigY] = eyeMovDetectTr(folderDir,char(folders),origX, origY, blinkTrs, samplerate, exptDat.exp.conds); 
    
 
    subjDat.(folders{i}).origx = origX;
    subjDat.(folders{i}).origy = origY;
    
    for c =1:size(cOrigX,3)
        for tr =1:size(cOrigX,2)
            xtr = cOrigX(:,tr,c);
            ytr = cOrigY(:,tr,c);
            
            xp = xtr - nanmean(xtr);
            yp = ytr - nanmean(ytr);
            
            normX = xp/max(abs(xp),[],'all');
            normY = yp/max(abs(yp),[],'all');
            xpTr(1:length(normX),tr,c) = normX;
            ypTr(1:length(normY),tr,c) = normY;
        end
        
    end
    subjDat.(folders{i}).normx = xpTr;
    subjDat.(folders{i}).normy = ypTr;
    
    dataPts = [0:1/1000:exptDat.exp.stimDur*6];
    phases = cos(2*pi*(1/exptDat.exp.stimDur)*dataPts);%cos(2*pi*cycPerSec*stimDur/(stimDur*flipsPerSec)-2*pi*flipTimes./stimDur).*360;
    
    conditions = {'Minbg', 'Meanbg'};
    xax = (1:size(xpTr,1))./1000;
    for n =1:size(xpTr,3)
        figure('Position',[100 100 1200 600]); % 
        plot(xax(1:length(phases)), phases,'LineWidth', 2,'col','k')
        hold on
        plot(xax,xpTr(:,:,n))
        set(gca,'box','off')
        set(gca, 'LineWidth',2)
        sgtitle(strcat(conditions{n}, sprintf(' %s ',folders{i})),'fontweight','bold','Interpreter', 'none')
        
        xlabel('time from stimulus onset (s)','fontsize',18,'fontweight','bold')
        ylabel('x position (normalized)','fontsize',18,'fontweight','bold')
        
        yticklab = get(gca,'YTickLabel');
        set(gca,'YTickLabel',yticklab,'fontsize',12)
        xticklab = get(gca,'XTickLabel');
        set(gca,'XTickLabel',xticklab,'fontsize',12)
        
        legend('stimulus','gaze')%,'location','bestoutside')
        plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/anal_plots/');
        %saveas(gcf,strcat(plotdir, sprintf('%s_%s_electrooculogram_10272022.png',folders{i},conditions{n})));
        saveas(gcf,strcat(plotdir, sprintf('%s_%s_noblink_nosacc_electrooculogram_10272022.svg',folders{i},conditions{n})));

    end
    
    
%  figure(); plot(xpTr(:,8,1))
% figure(); plot(origX(:,8,1))
        plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/anal_plots/');
         %saveas(gcf,strcat(plotdir, sprintf('%s_%s_electrooculogram_10272022.png',folders{i},conditions{n})));
         saveas(gcf,strcat(plotdir, sprintf('%s_%s_noblink_nosacc_electrooculogram_10272022.svg',folders{i},conditions{n})));
end

%% Sine-Cosine fits

fitCoefs = struct();
amp = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));
freq = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));
phase = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));
offset = nan(size(subjDat.(folders{1}).normx,2),length(folders),length(exptDat.exp.conds));

for i =1:length(folders)
    for c = 1:length(exptDat.exp.conds)
        a_sol = [];
        b_sol = [];
        c_sol = [];
        d_sol = [];
        for t =1:size(subjDat.(folders{i}).normx,2)
            y = subjDat.(folders{i}).normx(1:22001,t,c);
            yfit = y(3534:end);
            yfit = yfit(~isnan(y(3534:end)));
            x = 0:1/1000:exptDat.exp.blockLength;
            x = x(3534:length(y));
            xfit = x(~isnan(y(3534:end)));
            data.y = yfit';
            data.x = xfit;
            x0 = [1, pi/2,0, 0];
            L = @(params)CosMse(data, params);
            %L = @(params)CosNorm(data, params);
            [xsol,fval,exitflag,output] = fminsearch(@(params) L(params), x0);
            
            a_sol = [a_sol xsol(1)];
            b_sol = [b_sol xsol(2)];
            c_sol = [c_sol xsol(3)];
            d_sol = [d_sol xsol(4)];
            amp(t,i,c) = xsol(1);
            freq(t,i,c) = xsol(2);
            phase(t,i,c) = xsol(3);
            offset(t,i,c) = xsol(4);
        end
        fitCoefs.(folders{i}).(conditions{c}).amplitude = a_sol;
        fitCoefs.(folders{i}).(conditions{c}).frequency = b_sol./(2*pi);
        fitCoefs.(folders{i}).(conditions{c}).phase = c_sol;
        fitCoefs.(folders{i}).(conditions{c}).offset = d_sol;

    end
end

%plot and compare to stimulus a and b
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmap = flip(cmaps(2).map) ;

%p = figure('Position',[100 100 900 500]);
p = figure('Position',[100 100 1200 600]); %
[h,pos]=tight_subplot(1,4,[.02 .05],[.1 .15],[.05 .02]); %use tight plot for adjusting spacing between axes, and vertical and horizontal margin size

fit = {'Amplitude (normalized)', 'Frequency (cyc/s)', 'Phase (rad)', 'Offset (normalized)'};
fitNames = {'amp', 'freq', 'phase', 'offset'};
for f =1:length(fit)
    %plot the data points
    % and the mean ±95%CI
    eval(strcat('yvar =' , sprintf('%s',fitNames{f}),';'));
    mYvar = nanmean(yvar,[1,2]);
    
    %95%CI
    ci_high = mean(yvar,[1,2]) + std(yvar,[],[1,2]);%/sqrt(size(yvar,1)*size(yvar,2));
    ci_low = mean(yvar,[1,2]) - std(yvar,[],[1,2]);%/sqrt(size(yvar,1));
    yvarlab = {sprintf('%s',fit{f})};
    
    axes(h(f));
    for c =1:length(exptDat.exp.conds)
        xplot = c*ones(1, length(yvar(:,1,c))*length(yvar(1,:,c)));
        yplot =  reshape(yvar(:,:,c), 1,size(yvar(:,:,c),1)*size(yvar(:,:,c),2));
        scatter(xplot,yplot,20,'MarkerFaceColor',cmap(4,:), 'MarkerEdgeColor',cmap(4,:),'LineWidth',0.5);
        hold on;
        scatter(c*1,mYvar(c),60,'MarkerFaceColor',cmaps(1).map(4,:), 'MarkerEdgeColor',cmaps(1).map(4,:),'LineWidth',1.5); %mean
        hold on
        line([c*1 c*1], [ci_low(c) ci_high(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI vertical
        hold on
        line([c*1-0.1 c*1+0.1], [ci_low(c) ci_low(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI whiskers
        hold on
        line([c*1-0.1 c*1+0.1], [ci_high(c) ci_high(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI whiskers
        hold on

    ylim([min(yvar(:,:,c),[],'all')-0.1, max(yvar(:,:,c),[],'all')+0.1]);
    end
    
    % Set up axes.
    xlim([0, 3]);
    curtick = get(gca, 'xTick');
    xticks(unique(round(curtick)));
    ax = gca;
    ax.XTick = [1:length(conditions)];
    ax.XTickLabels = conditions;
    ylab = yvarlab{:};
    ylabel(ylab,'fontweight','bold','fontsize',16)
    yticklab = get(gca,'YTickLabel');
    set(gca,'YTickLabel',yticklab,'fontsize',12)
  
end
legend('Trial Model fit','Mean±SD')
sgtitle('fminsearch mse results', 'FontWeight', 'Bold','fontsize',22)
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/smooth_pursuit/anal_plots/');
saveas(gcf,strcat(plotdir, sprintf('%s_model_fminsearch_params_10272022.png',folders{:})));
   
%% Show trial examples with fit


