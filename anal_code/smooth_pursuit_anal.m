%% load data
% 
% eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/ET_pilot_nofix.mat');
% exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/s1_smooth_pursuit_v1_sn1_rn1_date20221011T131543_nofix.mat');
% 
% xbsl = exptDat.exp.resolution.width/2;
% ybsl = exptDat.exp.resolution.height/2;
% xDat = eyeDat.EyeData.mx-xbsl;
% yDat = eyeDat.EyeData.my-ybsl;
% 
% %screen dimensions
% scrpixelwidth = exptDat.exp.resolution.width; 
% scrpixelheight = exptDat.exp.resolution.height; 
% 
% pixperdeg = exptDat.exp.ppd;
% halfscrwidth = scrpixelwidth/pixperdeg(1);
% scrheight = scrpixelheight/abs(pixperdeg(1));
% fixDiameterhorz = exptDat.exp.fixSizeDeg;
% fixDiameterVert = exptDat.exp.fixSizeDeg;
% 
% %draw screen
% xcenter =0;
% ycenter =0;
% [scrxs, scrys] = drawSquare(halfscrwidth, scrheight, xcenter,ycenter);
% [fxs,fys] = drawEllipse(fixDiameterhorz,fixDiameterVert,xcenter,ycenter);
% x = xDat/pixperdeg;
% y = yDat/pixperdeg;
% 
% figure('Position',[100 100 1800 800]);
% for c =1:length(exptDat.exp.conds)
%     subplot(1,length(exptDat.exp.conds),c)
%     plot(x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)),y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)))
%     hold on
%     plot(scrxs, scrys, 'k--', 'linewidth',2); %plot screen limits
%     hold on
%     plot(fxs, fys, 'r-','linewidth', 2); %plot fig limits
%     set(gca,'box','off');
%     set(gca,'linewidth', 2);
%     axis equal
%     xlim([-3 3])
%     ylim([-3 3])
%     title(strcat(exptDat.exp.conds{c}))
%     xlabel('x position (deg)')
%     ylabel('y position (deg)')
% 
% end
% plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots');
% saveas(gcf,strcat(plotdir, sprintf('%s_electrooculogram_minbg_10112022.png',filename)));
% saveas(gcf,strcat(plotdir, sprintf('%s_electrooculogram_minbg_10112022.svg',filename)));


%%
% Loic
fix = 0;
mini = 0;
if mini
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/ET_pilot_nofix.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/s1_smooth_pursuit_v1_sn1_rn1_date20221011T131543_nofix.mat');
    savename = 'min';
else
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_mean/ET_pilot_nofix.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_mean/s1_smooth_pursuit_v1_sn1_rn1_date20221011T121929_nofix.mat');
    savename = 'mean';
end

% 
xbsl = exptDat.exp.resolution.width/2;
ybsl = exptDat.exp.resolution.height/2;
xDat = eyeDat.EyeData.mx-xbsl; %center 0 relative to middle of the screen 
yDat = eyeDat.EyeData.my-ybsl; %center 0 relative to middle of the screen 
pixperdeg = exptDat.exp.ppd;
x = xDat/pixperdeg; %convert from pixels to degrees
y = yDat/pixperdeg; %convert from pixels to degrees

scrpixelwidth = exptDat.exp.resolution.width; 
scrpixelheight = exptDat.exp.resolution.height; 
halfscrwidth = scrpixelwidth/pixperdeg(1);
scrheight = scrpixelheight/abs(pixperdeg(1));
fixDiameterhorz = exptDat.exp.fixSizeDeg;
fixDiameterVert = exptDat.exp.fixSizeDeg;

%draw screen
xcenter =0;
ycenter =0;
[fxs,fys] = drawEllipse(fixDiameterhorz/3,fixDiameterVert/3,xcenter,ycenter);


if length(eyeDat.EyeData.Fixated) ~= length(exptDat.exp.longFormConds)
    eyeDat.EyeData.Fixated = [eyeDat.EyeData.Fixated nan(1, length(exptDat.exp.longFormConds)-length(eyeDat.EyeData.Fixated))];
end
conditions = {'Vertical Square', 'Vertical Rectangle'};
gratings.squ = reshape([squeeze(exptDat.exp.LWave(1,:,:,1))], [107 107]);
gratings.rect = reshape(squeeze(exptDat.exp.rectLWave(1,:,:,1)), [107 214]);

figure('Position',[100 100 1200 600]); %
[h,pos]=tight_subplot(1,2,[.02 .02],[.1 .15],[.03 .02]); %use tight plot for adjusting spacing between axes, and vertical and horizontal margin size

thresh = 0.99; %eye position quantile detection threshold
yc = 0.68;
voffset = 0.08; %vertical offset
hoffset = 0.33; %horizontal offset

for c =1:length(exptDat.exp.conds)
    axes(h(c));
    %     xp = x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c))- nanmean(x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)));
    %     yp = y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c))- nanmean(y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)));
    xcond = x(exptDat.exp.longFormConds == c);
    ycond =  y(exptDat.exp.longFormConds == c);
    xp = xcond(abs(ycond)<quantile(abs(ycond),thresh) & abs(xcond)<quantile(abs(xcond),thresh))- nanmean(xcond(abs(ycond)<quantile(abs(ycond),thresh) & abs(xcond)<quantile(abs(xcond),thresh)));
    yp =  ycond(abs(ycond)<quantile(abs(ycond),thresh) & abs(xcond)<quantile(abs(xcond),thresh))- nanmean(ycond(abs(ycond)<quantile(abs(ycond),thresh) & abs(xcond)<quantile(abs(xcond),thresh)));

    normX = xp/max(abs(xp),[],'all');
    normY = yp/max(abs(yp),[],'all');
    plot(normX,normY)
    hold on
    plot(fxs, fys, 'r-','linewidth', 2); %plot fix limits
    set(gca,'box','off');
    set(gca,'linewidth', 2);
    axis equal
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    title({strcat(conditions{c});''},'fontsize',20)
    xlabel('x position (normalized)','fontweight','bold','fontsize',18)
    ylabel('y position (normalized)','fontweight','bold','fontsize',18)
    xticklab = get(gca,'XTickLabel');
    set(gca,'XTickLabel',xticklab,'fontsize',14)
    yticklab = get(gca,'YTickLabel');
    set(gca,'YTickLabel',yticklab,'fontsize',14)
    
    %%%% add drawing of stimulus condition to each subplot
    ax1 = get(h,'position'); % Store handle to axes 1.
    % Create smaller axes in top right, and plot on it
    % Store handle to axes 2 in ax2.
    
    % add background
    height = 0.30;
    width = 0.30;
    ax2 = axes('Position',[pos{c}(1)+hoffset-0.1 yc-voffset-0.05 width height]);
    box on;
    hold(ax2); %
    
    if mini 
    bg = repmat(min(gratings.squ(:,:)./max(gratings.squ(:,:),[],'all'),[],'all'),size(gratings.squ,1)*4,size(gratings.rect,2)*1.5);
    else 
         bg = repmat(mean(gratings.squ(:,:)./max(gratings.squ(:,:),[],'all'),'all'),size(gratings.squ,1)*4,size(gratings.rect,2)*1.5);
    end
    if fix == 0
        imshow(bg);
%     elseif fix == 1
%         [dfxs,dfys] = drawEllipse(30,30,160,54);
%         imshow(bg);
%         hold on
%         plot(dfxs, dfys, 'w-','linewidth', 2); %plot fix limits
%         
    end
    %add stimuli
    if contains(conditions{c},'Squ')
        ax2 = axes('Position',[pos{c}(1)+hoffset yc+voffset .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings.squ(:,:)./max(gratings.squ(:,:),[],'all');
        imshow(grat);
        
        ax2 = axes('Position',[pos{c}(1)+hoffset yc-voffset .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings.squ(:,:)./max(gratings.squ(:,:),[],'all');
        imshow(grat);
    elseif contains(conditions{c},'Rect')
        ax2 = axes('Position',[pos{c}(1)+hoffset yc+voffset .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings.rect(:,:)./max(gratings.rect(:,:),[],'all');
        imshow(grat);
        
        ax2 = axes('Position',[pos{c}(1)+hoffset yc-voffset .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings.rect(:,:)./max(gratings.rect(:,:),[],'all');
        imshow(grat);
        
    end
end
sgtitle(sprintf('Smooth pursuit of visual phantom illusion'),'fontweight','bold','fontsize',22)

plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
saveas(gcf,strcat(plotdir, sprintf('%sbg_99percthr_electrooculogram_loic10112022.png',savename)));
saveas(gcf,strcat(plotdir, sprintf('%sbg_99percthr_electrooculogram_loic10112022.svg',savename)));

%% Look at eye movements over time
% Loic
fix = 0;
mini = 1;
if mini
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/ET_pilot_nofix.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/s1_smooth_pursuit_v1_sn1_rn1_date20221011T131543_nofix.mat');
    savename = 'min';
else
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_mean/ET_pilot_nofix.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_mean/s1_smooth_pursuit_v1_sn1_rn1_date20221011T121929_nofix.mat');
    savename = 'mean';
end

% 
xbsl = exptDat.exp.resolution.width/2;
ybsl = exptDat.exp.resolution.height/2;
xDat = eyeDat.EyeData.mx-xbsl; %center 0 relative to middle of the screen 
yDat = eyeDat.EyeData.my-ybsl; %center 0 relative to middle of the screen 
pixperdeg = exptDat.exp.ppd;
x = xDat/pixperdeg; %convert from pixels to degrees
y = yDat/pixperdeg; %convert from pixels to degrees

xpblock = nan(647,10, length(exptDat.exp.conds));
ypblock = nan(647,10, length(exptDat.exp.conds));


thresh = 0.99; %eye position quantile detection threshold
for c =1:length(exptDat.exp.conds)
    blockdiff = unique([find(diff(exptDat.exp.longFormConds == c)) length(exptDat.exp.longFormConds)]); %get block onset time and offset time
    for n =1:length(blockdiff)/2
        xbloc = x(blockdiff(2*(n-1)+1):blockdiff(2*n));
        ybloc = y(blockdiff(2*(n-1)+1):blockdiff(2*n));
        xp = xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
        yp = ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
        normX = xp/max(abs(xp),[],'all');
        normY = yp/max(abs(yp),[],'all');
        xpblock(1:length(normX),n,c) = normX;
        ypblock(1:length(normY),n,c) = normY;
    end

end

%% Plot x over time
% oscillatory signal function used to generate the temporally varying speed
% of the drifting stimulus
% flipsPerSec = exptDat.exp.flipsPerSec; %60;
% cycles = exptDat.exp.stim.cycles;%2;
% cycPerSec = exptDat.exp.stim.cycPerSec;%1.13;
% stimDur = (cycles/cycPerSec)*2;
% fliptimes = 1-1:size(exptDat.exp.LWave,1)-1;
phases = exptDat.exp.stim.phases;%cos(2*pi*cycPerSec*stimDur/(stimDur*flipsPerSec)-2*pi*flipTimes./stimDur).*360;

conditions = {'Square', 'Rectangle'};
xax = (1:size(xpblock,1))./60;
for n =1:size(xpblock,3)
    figure('Position',[100 100 1200 600]); %
    [h,pos]=tight_subplot(1,3,[.02 .05],[.1 .15],[.05 .02]); %use tight plot for adjusting spacing between axes, and vertical and horizontal margin size
    
    for i =1:3 %size(xpblock,2)
        %subplot(size(xpblock,3),size(xpblock,2)/2,i)
        %subplot(1,3,i)
        axes(h(i));
        plot(xax,xpblock(:,i,n))
        hold on
        plot(xax(1:length(repmat(phases,1,3)./360)),repmat(phases,1,3)./360,'LineWidth', 2,'col','k')
        set(gca,'box','off')
        set(gca, 'LineWidth',2)
        sgtitle(strcat(conditions{n}, sprintf(' %s background',savename)),'fontweight','bold')
        if i ==1
            xlabel('time from stimulus onset (s)','fontsize',18,'fontweight','bold')
            ylabel('x position (normalized)','fontsize',18,'fontweight','bold')
            
            yticklab = get(gca,'YTickLabel');
            set(gca,'YTickLabel',yticklab,'fontsize',12)
            xticklab = get(gca,'XTickLabel');
            set(gca,'XTickLabel',xticklab,'fontsize',12)
        elseif i ==3
            legend('gaze', 'stimulus')%,'location','bestoutside')
        end
        
    end
    plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
    saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_99percthr_electrooculogram_loic10112022.png',conditions{n},savename)));
    saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_99percthr_electrooculogram_loic10112022.svg',conditions{n},savename)));
end

%% analyze data of 10/18/2022

%1) let's have a first impression of the data
folderDir = '/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/smooth_pursuit_10182022/';
folders = dir(folderDir);
folders = {folders(~contains({folders.name},'.')).name};

for i =1:length(folders)
    
    fileNames = dir(strcat(folderDir, folders{i}));
    fileNames = string({fileNames(contains({fileNames.name},'smooth')).name});
    eyeDat = load(strcat(folderDir, folders{i},'/',fileNames(contains(fileNames,'eyeDat'))));
    exptDat = load(strcat(folderDir, folders{i},'/',fileNames(~contains(fileNames,'eyeDat'))));
    
    xbsl = exptDat.exp.resolution.width/2;
    ybsl = exptDat.exp.resolution.height/2;
    xDat = eyeDat.EyeData.mx-xbsl; %center 0 relative to middle of the screen
    yDat = eyeDat.EyeData.my-ybsl; %center 0 relative to middle of the screen
    pixperdeg = exptDat.exp.ppd;
    x = xDat/pixperdeg; %convert from pixels to degrees
    y = yDat/pixperdeg; %convert from pixels to degrees
    
    xpblock = nan(647,10, length(exptDat.exp.conds));
    ypblock = nan(647,10, length(exptDat.exp.conds));
    
    
    thresh = 0.99; %eye position quantile detection threshold
    for c =1:length(exptDat.exp.conds)
        blockdiff = unique([find(diff(exptDat.exp.longFormConds == c)) length(exptDat.exp.longFormConds)]); %get block onset time and offset time
        for n =1:length(blockdiff)/2
            xbloc = x(blockdiff(2*(n-1)+1):blockdiff(2*n));
            ybloc = y(blockdiff(2*(n-1)+1):blockdiff(2*n));
            xp = xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
            yp = ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
            normX = xp/max(abs(xp),[],'all');
            normY = yp/max(abs(yp),[],'all');
            xpblock(1:length(normX),n,c) = normX;
            ypblock(1:length(normY),n,c) = normY;
        end
        
    end
    
    phases = exptDat.exp.stim.oscillation;%cos(2*pi*cycPerSec*stimDur/(stimDur*flipsPerSec)-2*pi*flipTimes./stimDur).*360;
    
    conditions = {'Square'};
    xax = (1:size(xpblock,1))./60;
    for n =1:size(xpblock,3)
        figure('Position',[100 100 1200 600]); %
        [h,pos]=tight_subplot(1,3,[.02 .05],[.1 .15],[.05 .02]); %use tight plot for adjusting spacing between axes, and vertical and horizontal margin size
        
        for j =1:3 %size(xpblock,2)
            %subplot(size(xpblock,3),size(xpblock,2)/2,i)
            %subplot(1,3,i)
            axes(h(j));
            plot(xax,xpblock(:,j,n))
            hold on
            plot(xax(1:length(repmat(phases,1,3))),repmat(phases,1,3),'LineWidth', 2,'col','k')
            set(gca,'box','off')
            set(gca, 'LineWidth',2)
            %sgtitle(strcat(conditions{n}, sprintf(' %s background',savename)),'fontweight','bold')
            if j ==1
                xlabel('time from stimulus onset (s)','fontsize',18,'fontweight','bold')
                ylabel('x position (normalized)','fontsize',18,'fontweight','bold')
                
                yticklab = get(gca,'YTickLabel');
                set(gca,'YTickLabel',yticklab,'fontsize',12)
                xticklab = get(gca,'XTickLabel');
                set(gca,'XTickLabel',xticklab,'fontsize',12)
            elseif j ==3
                legend('gaze', 'stimulus')%,'location','bestoutside')
            end
            
        end
%         plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
%         saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_99percthr_electrooculogram_loic10112022.png',conditions{n},savename)));
%         saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_99percthr_electrooculogram_loic10112022.svg',conditions{n},savename)));
    end
    
end

%% Plot all trials data for each subject

%compute average response curves of each subject
meanx = nan(647, length(folders));
for i =1:length(folders)
    
    fileNames = dir(strcat(folderDir, folders{i}));
    fileNames = string({fileNames(contains({fileNames.name},'smooth')).name});
    eyeDat = load(strcat(folderDir, folders{i},'/',fileNames(contains(fileNames,'eyeDat'))));
    exptDat = load(strcat(folderDir, folders{i},'/',fileNames(~contains(fileNames,'eyeDat'))));
    
    xbsl = exptDat.exp.resolution.width/2;
    ybsl = exptDat.exp.resolution.height/2;
    xDat = eyeDat.EyeData.mx-xbsl; %center 0 relative to middle of the screen
    yDat = eyeDat.EyeData.my-ybsl; %center 0 relative to middle of the screen
    pixperdeg = exptDat.exp.ppd;
    x = xDat/pixperdeg; %convert from pixels to degrees
    y = yDat/pixperdeg; %convert from pixels to degrees
    
    xpblock = nan(647,10, length(exptDat.exp.conds));
    ypblock = nan(647,10, length(exptDat.exp.conds));
    
    
    thresh = 0.99; %eye position quantile detection threshold
    for c =1:length(exptDat.exp.conds)
        blockdiff = unique([find(diff(exptDat.exp.longFormConds == c)) length(exptDat.exp.longFormConds)]); %get block onset time and offset time
        for n =1:length(blockdiff)/2
            xbloc = x(blockdiff(2*(n-1)+1):blockdiff(2*n));
            ybloc = y(blockdiff(2*(n-1)+1):blockdiff(2*n));
            xp = xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
            yp = ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
            normX = xp/max(abs(xp),[],'all');
            normY = yp/max(abs(yp),[],'all');
            xpblock(1:length(normX),n,c) = normX;
            ypblock(1:length(normY),n,c) = normY;
        end
        
    end
    
    subjDat.(folders{i}).x = xpblock;
    meanx(:,i) = mean(subjDat.(folders{i}).x,2); 
    subjDat.(folders{i}).y = ypblock;
    
    phases = exptDat.exp.stim.oscillation;%cos(2*pi*cycPerSec*stimDur/(stimDur*flipsPerSec)-2*pi*flipTimes./stimDur).*360;
    
    conditions = {'Square'};
    xax = (1:size(xpblock,1))./60;
    for n =1:size(xpblock,3)
        figure('Position',[100 100 1200 600]); %
        
        plot(xax,xpblock(:,:,n))
        hold on
        plot(xax(1:length(repmat(phases,1,3))),repmat(phases,1,3),'LineWidth', 2,'col','k')
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
    
%         plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
%         saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_99percthr_electrooculogram_10182022.png',conditions{n},folders{i})));
%         saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_99percthr_electrooculogram_10182022.svg',conditions{n},folders{i})));
end


%% Fit a cosine function to the data using fminsearch
% use fminsearch to find optimal parameters of a cosine by minimizing the difference
%between the data and a cosine function


%%troubleshoot
y = meanx(107:end,1);
y = y(~isnan(y));
x = 0:1/60:exptDat.exp.blockLength;
x = x(1:length(meanx));
x = x(107:end);
x = x(~isnan(meanx(107:end,1)));
fun = @(a,b,x) a*cos(b*x); %cos(2*pi*(1/exptDat.exp.stimDur)*exptDat.exp.flipTimes)

%obj_fun = @(params) norm(fun(params(1), params(2), x)- y);
obj_fun = @(params) (y'-fun(params(1), params(2), x)).^2;
%obj_fun = @(params) (y-fun(params(1), params(2), x));

%x0 =rand(1,2);
x0 = [1, pi/2.5];
[xsol,fval, exitflag] = fminsearch(obj_fun,x0);

a_sol = xsol(1);
b_sol = xsol(2);

figure;
plot(x, y, '+', 'MarkerSize', 10, 'LineWidth', 2)
hold on
plot(x, fun(a_sol, b_sol, x), '-')

%% fit model using mean square error minimization

for i =1:length(folders)
    y = meanx(107:end,i);
    y = y(~isnan(y));
    x = 0:1/60:exptDat.exp.blockLength;
    x = x(1:length(meanx));
    x = x(107:end);
    x = x(~isnan(meanx(107:end,i)));
    data.y = y';
    data.x = x;
    x0 = [1, pi/2];
    L = @(params)CosFun(data, params);
    [xsol,fval,exitflag,output] = fminsearch(@(params) L(params), x0);
    
    a_sol = xsol(1);
    b_sol = xsol(2);
    
    figure('Position',[100 100 1200 600]);
    plot(x, y, 'o', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
    plot(x, a_sol*cos(b_sol*x), '-','LineWidth', 2)
    hold on
    plot(xax(1:length(repmat(phases,1,3))),repmat(phases,1,3),'LineWidth', 2,'col','k')
    hold on 
    text(8,0.5, sprintf('y = %.2f*cos(%.2f*x)',a_sol,b_sol),'fontsize',18)
    set(gca,'box','off')
    set(gca, 'LineWidth',2)
    sgtitle(strcat(conditions{:}, sprintf(' %s ',folders{i})),'fontweight','bold','Interpreter', 'none')
    
    xlabel('time from stimulus onset (s)','fontsize',18,'fontweight','bold')
    ylabel('x position (normalized)','fontsize',18,'fontweight','bold')
    
    yticklab = get(gca,'YTickLabel');
    set(gca,'YTickLabel',yticklab,'fontsize',12)
    xticklab = get(gca,'XTickLabel');
    set(gca,'XTickLabel',xticklab,'fontsize',12)
    
    legend('gaze', 'model', 'stimulus')%,'location','bestoutside')
    plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
    saveas(gcf,strcat(plotdir, sprintf('%s_%sbg_mean_model_fminsearch_10182022.png',conditions{:},folders{i})));
     
end


%% Compute stdx and stdy and compare across conditions with mean or min background luminance
% stdx = nan(10, length(exptDat.exp.conds));
% stdy = nan(10, length(exptDat.exp.conds));
% for c =1:length(exptDat.exp.conds)
% 
%     blockdiff = unique([find(diff(exptDat.exp.longFormConds == c)) length(exptDat.exp.longFormConds)]); %get block onset time and offset time
%     for n =1:length(blockdiff)/2
%         xbloc = x(blockdiff(2*(n-1)+1):blockdiff(2*n));
%         ybloc = y(blockdiff(2*(n-1)+1):blockdiff(2*n));
% 
%         xp = xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
%         yp = ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
%         normX = xp/max(abs(xp),[],'all');
%         normY = yp/max(abs(yp),[],'all');
%         stdx(n,c) = std(normX);
%         stdy(n,c) = std(normY);
%     end
% 
% end
% 
% conds = [{'Square'}, {'Rectangle'}];
% 
% % ordconds = [{'BothHor'},{'BothVert'},{'LeftHor'},{'RightHor'}];
% % ordstdx = [stdx(:,4) stdx(:,3) stdx(:,1) stdx(:,2)]; %need to reorganize in alphabetical order
% %%plot stdx
% p = spursuitJitterCiPlot(conds, stdx,  'Condition', 'Standard deviation (normalized)', 'Standard deviation of horizontal eye position', gratings, fix,mini);
% 
% 
% plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
% saveas(gcf,strcat(plotdir, sprintf('%sbg_99percthr_norm_stdx_loic10112022.png',savename)));
% saveas(gcf,strcat(plotdir, sprintf('%sbg_99percthr_norm_stdx_loic10112022.svg',savename)));
% 
% %%plot stdy
% %ordstdy = [stdy(:,4) stdy(:,3) stdy(:,1) stdy(:,2)]; %need to reorganize in alphabetical order
% plot = spursuitJitterCiPlot(conds, stdy,  'Condition', 'Standard deviation (normalized)', 'Standard deviation of vertical eye position', gratings, fix,mini);
% plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
% saveas(gcf,strcat(plotdir, sprintf('%sbg_99percthr_norm_stdy_loic10112022.png',savename)));
% saveas(gcf,strcat(plotdir, sprintf('%sbg_99percthr_norm_stdy_loic10112022.svg',savename)));
% 



%% stats
%load and preprocess data

% %load minbg condition
% eyeDat.minbg = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/ET_pilot_nofix.mat');
% exptDat.minbg = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_min/s1_smooth_pursuit_v1_sn1_rn1_date20221011T131543_nofix.mat');
% 
% xbsl = exptDat.minbg.exp.resolution.width/2;
% ybsl = exptDat.minbg.exp.resolution.height/2;
% xDat = eyeDat.minbg.EyeData.mx-xbsl; %center 0 relative to middle of the screen 
% yDat = eyeDat.minbg.EyeData.my-ybsl; %center 0 relative to middle of the screen 
% pixperdeg = exptDat.minbg.exp.ppd;
% xmin = xDat/pixperdeg; %convert from pixels to degrees
% ymin = yDat/pixperdeg; %convert from pixels to degrees
% 
% %preprocess
% thresh = 0.99; %eye position quantile detection threshold
% stdx.minbg = nan(10, length(exptDat.minbg.exp.conds));
% stdy.minbg = nan(10, length(exptDat.minbg.exp.conds));
% for c =1:length(exptDat.minbg.exp.conds)
%     
%     
%     blockdiff = unique([find(diff(exptDat.minbg.exp.longFormConds == c)) length(exptDat.minbg.exp.longFormConds)]); %get block onset time and offset time
%     for n =1:length(blockdiff)/2
%         xbloc = xmin(blockdiff(2*(n-1)+1):blockdiff(2*n));
%         ybloc = ymin(blockdiff(2*(n-1)+1):blockdiff(2*n));
% 
%         xp = xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
%         yp = ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
%  
%         normX = xp/max(abs(xp),[],'all');
%         normY = yp/max(abs(yp),[],'all');
% %         figure();
% %         plot(normX, normY)
% 
%         stdx.minbg(n,c) = std(normX);
%         stdy.minbg(n,c) = std(normY);
%     end
% 
% end
% 
% % load meanbg condition
% eyeDat.meanbg = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_mean/ET_pilot_nofix.mat');
% exptDat.meanbg = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_10112022/bg_mean/s1_smooth_pursuit_v1_sn1_rn1_date20221011T121929_nofix.mat');
% xbsl = exptDat.meanbg.exp.resolution.width/2;
% ybsl = exptDat.meanbg.exp.resolution.height/2;
% xDat = eyeDat.meanbg.EyeData.mx-xbsl; %center 0 relative to middle of the screen 
% yDat = eyeDat.meanbg.EyeData.my-ybsl; %center 0 relative to middle of the screen 
% pixperdeg = exptDat.meanbg.exp.ppd;
% xmean = xDat/pixperdeg; %convert from pixels to degrees
% ymean = yDat/pixperdeg; %convert from pixels to degrees
% %preprocess
% 
% stdx.meanbg = nan(10, length(exptDat.meanbg.exp.conds));
% stdy.meanbg = nan(10, length(exptDat.meanbg.exp.conds));
% for c =1:length(exptDat.meanbg.exp.conds)
%     
%     blockdiff = unique([find(diff(exptDat.meanbg.exp.longFormConds == c)) length(exptDat.meanbg.exp.longFormConds)]); %get block onset time and offset time
%     for n =1:length(blockdiff)/2
%         xbloc = xmean(blockdiff(2*(n-1)+1):blockdiff(2*n));
%         ybloc = ymean(blockdiff(2*(n-1)+1):blockdiff(2*n));
% 
%         xp = xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(xbloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
%         yp = ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh))- nanmean(ybloc(abs(ybloc)<quantile(abs(ybloc),thresh) & abs(xbloc)<quantile(abs(xbloc),thresh)));
%         normX = xp/max(abs(xp),[],'all');
%         normY = yp/max(abs(yp),[],'all');
% %         figure();
% %         plot(normX,normY)
% 
%         stdx.meanbg(n,c) = std(normX);
%         stdy.meanbg(n,c) = std(normY);
%     end
% 
% end
% 
% % compare square and rectangle
% % compare square minbg to square meanbg
% %= test for main effect in ANOVA2
% conds = [{'Square'}, {'Rectangle'}, {'minbg'},{'meanbg'}];
% aovx = anova2([stdx.minbg; stdx.meanbg],10, 'on');
% aovy = anova2([stdy.minbg; stdy.meanbg],10, 'on');
% 
% % compare rectangle minbg to rectangle meanbg


