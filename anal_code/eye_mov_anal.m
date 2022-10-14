
%% this script is intended to detect eye movements from the phantom illusion task and compare the effect of phantom illusion 
% based on whether the grating pair was drifting vertically or whether the
% pair was drifting horizontally.
%started on 09-15-2022


%% load data

eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/ET_pilot_09142022.mat');
exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/s1_fillingin_rsvp_v1_sn1_rn1_date20220914T130528.mat');


%% detect eye movements
dir = '/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/';
filename ='ET_pilot_09142022';
xbsl = exptDat.exp.resolution.width/2;
ybsl = exptDat.exp.resolution.height/2;
xDat = eyeDat.EyeData.mx-xbsl;
yDat = eyeDat.EyeData.my-ybsl;
%blinks = (~eyeDat.EyeData.Fixated)';
blinks = zeros(length(yDat),1)';
samplerate = 60;
eyeMovDat = eyeMovDetect(dir,filename,xDat, yDat, blinks, samplerate);

%% Plot eye movement data
%screen dimensions
scrpixelwidth = exptDat.exp.resolution.width; 
scrpixelheight = exptDat.exp.resolution.height; 

pixperdeg = exptDat.exp.ppd;
halfscrwidth = scrpixelwidth/pixperdeg(1);
scrheight = scrpixelheight/abs(pixperdeg(1));
fixDiameterhorz = exptDat.exp.fixSizeDeg;
fixDiameterVert = exptDat.exp.fixSizeDeg;

%draw screen
xcenter =0;
ycenter =0;
[scrxs, scrys] = drawSquare(halfscrwidth, scrheight, xcenter,ycenter);
[fxs,fys] = drawEllipse(fixDiameterhorz,fixDiameterVert,xcenter,ycenter);
x = xDat/pixperdeg;
y = yDat/pixperdeg;

figure('Position',[100 100 1800 800]);
for c =1:length(exptDat.exp.conds)
    subplot(1,length(exptDat.exp.conds),c)
    plot(x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)),y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)))
    hold on
    plot(scrxs, scrys, 'k--', 'linewidth',2); %plot screen limits
    hold on
    plot(fxs, fys, 'r-','linewidth', 2); %plot fig limits
    set(gca,'box','off');
    set(gca,'linewidth', 2);
    axis equal
    xlim([-3 3])
    ylim([-3 3])
    title(strcat(exptDat.exp.conds{c}))
    xlabel('x position (deg)')
    ylabel('y position (deg)')

end
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots');
saveas(gcf,strcat(plotdir, sprintf('%s_electrooculogram_fix.png',filename)));
saveas(gcf,strcat(plotdir, sprintf('%s_electrooculogram_fix.svg',filename)));



%% Compare data with fixation circle and data without
% Lasya
fix =1;
%fix
if fix == 1
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Lasya/ET_pilot_fix.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Lasya/s1_fillingin_rsvp_v1_sn1_rn1_date20220916T153317_fix.mat');
    state = 'ON';
    savename = 'fix';
else
    %no fix
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Lasya/ET_pilot_nofix.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Lasya/s1_fillingin_rsvp_v1_sn1_rn2_date20220916T152200_nofix.mat');
    state = 'OFF';
    savename = 'nofix';
end
xbsl = exptDat.exp.resolution.width/2;
ybsl = exptDat.exp.resolution.height/2;
xDat = eyeDat.EyeData.mx-xbsl;
yDat = eyeDat.EyeData.my-ybsl;


scrpixelwidth = exptDat.exp.resolution.width; 
scrpixelheight = exptDat.exp.resolution.height; 

pixperdeg = exptDat.exp.ppd;
halfscrwidth = scrpixelwidth/pixperdeg(1);
scrheight = scrpixelheight/abs(pixperdeg(1));
fixDiameterhorz = exptDat.exp.fixSizeDeg;
fixDiameterVert = exptDat.exp.fixSizeDeg;

%draw screen
xcenter =0;
ycenter =0;
[scrxs, scrys] = drawSquare(halfscrwidth, scrheight, xcenter,ycenter);
[fxs,fys] = drawEllipse(fixDiameterhorz,fixDiameterVert,xcenter,ycenter);
x = xDat/pixperdeg;
y = yDat/pixperdeg;

conditions = {'Horizontally oriented grating, left, drifting vertically', 'Horizontally oriented grating, right, drifting vertically',...
    'Vertically oriented gratings, both sides, drifting horizontally', 'Horizontally oriented gratings, both sides, drifting vertically'};
gratings = reshape([squeeze(exptDat.exp.LWave(1,:,:,2)) squeeze(exptDat.exp.LWave(1,:,:,2)) squeeze(exptDat.exp.LWave(1,:,:,1)) squeeze(exptDat.exp.LWave(1,:,:,2))], [107 107 4]);

figure('Position',[100 100 1800 800]);
for c =1:length(exptDat.exp.conds)
    %subplot(2, length(exptDat.exp.conds),c)
    %axis('on', 'image');
    h = subplot(1, length(exptDat.exp.conds),c);
    xp = x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c))- nanmean(x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)));
    yp = y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c))- nanmean(y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)));
    plot(xp,yp)
    hold on
    plot(scrxs, scrys, 'k--', 'linewidth',2); %plot screen limits
    hold on
    plot(fxs, fys, 'r-','linewidth', 2); %plot fig limits
    set(gca,'box','off');
    set(gca,'linewidth', 2);
    axis equal
    xlim([-3 3])
    ylim([-3 3])
    title(strcat(conditions{c}))
    xlabel('x position (deg)')
    ylabel('y position (deg)')
    ax1 = get(h,'position'); % Store handle to axes 1.
    % Create smaller axes in top right, and plot on it
    % Store handle to axes 2 in ax2.
    ax2 = axes('Position',[c*0.205 .58 .08 .08]);
    box on;
    hold(ax2); % 
    grat = gratings(:,:,c)./max(gratings(:,:,c),[],'all');
    imshow(grat);
    sgtitle(sprintf('Fixation %s',state))

end
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
saveas(gcf,strcat(plotdir, sprintf('electrooculogram_%s_lasya.png', savename)));
saveas(gcf,strcat(plotdir, sprintf('electrooculogram_%s_lasya.svg', savename)));

% Loic
fix = 0;
%fix
if fix == 1
    %     eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09192022/ET_pilot_fix.mat');
    %     exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09192022/s1_fillingin_rsvp_v1_sn1_rn1_date20220919T175820_fix.mat');
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09202022/ET_pilot_fix.mat'); 
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09202022/s1_fillingin_rsvp_v1_sn1_rn1_date20220920T174904_fix.mat');
    state = 'ON';
    savename = 'fix';
else
    %no fix
    %     eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09192022/ET_pilot_nofix.mat');
    %     exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09192022/s1_fillingin_rsvp_v1_sn1_rn1_date20220919T182246_nofix.mat');
    eyeDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09202022/ET_pilot_nofix2.mat');
    exptDat = load('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09202022/s1_fillingin_rsvp_v1_sn1_rn1_date20220920T180138_nofix.mat');
    
    state = 'OFF';
    savename = 'nofix';
end
% 
xbsl = exptDat.exp.resolution.width/2;
ybsl = exptDat.exp.resolution.height/2;
xDat = eyeDat.EyeData.mx-xbsl;
yDat = eyeDat.EyeData.my-ybsl;

scrpixelwidth = exptDat.exp.resolution.width; 
scrpixelheight = exptDat.exp.resolution.height; 

pixperdeg = exptDat.exp.ppd;
halfscrwidth = scrpixelwidth/pixperdeg(1);
scrheight = scrpixelheight/abs(pixperdeg(1));
fixDiameterhorz = exptDat.exp.fixSizeDeg;
fixDiameterVert = exptDat.exp.fixSizeDeg;

%draw screen
xcenter =0;
ycenter =0;
[fxs,fys] = drawEllipse(fixDiameterhorz/3,fixDiameterVert/3,xcenter,ycenter);
x = xDat/pixperdeg;
y = yDat/pixperdeg;

if length(eyeDat.EyeData.Fixated) ~= length(exptDat.exp.longFormConds)
    eyeDat.EyeData.Fixated = [eyeDat.EyeData.Fixated nan(1, length(exptDat.exp.longFormConds)-length(eyeDat.EyeData.Fixated))];
end
conditions = {'Horizontally oriented grating, left', 'Horizontally oriented grating, right',...
    'Vertically oriented gratings, both sides', 'Horizontally oriented gratings, both sides'};
gratings = reshape([squeeze(exptDat.exp.LWave(1,:,:,2)) squeeze(exptDat.exp.LWave(1,:,:,2)) squeeze(exptDat.exp.LWave(1,:,:,1)) squeeze(exptDat.exp.LWave(1,:,:,2))], [107 107 4]);

figure('Position',[100 100 1800 550]); %
[h,pos]=tight_subplot(1,4,[.03 .03],[.01 .01],[.03 .01]); %use tight plot for adjusting spacing between axes, and vertical and horizontal margin size
loffset = 0.112;
roffset = 0.164;
for c =1:length(exptDat.exp.conds)
    axes(h(c));
    xp = x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c))- nanmean(x((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)));
    yp = y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c))- nanmean(y((eyeDat.EyeData.Fixated == 1) & (exptDat.exp.longFormConds == c)));
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
    title({strcat(conditions{c});''},'fontsize',18)
    xlabel('x position (normalized)','fontweight','bold','fontsize',16)
    ylabel('y position (normalized)','fontweight','bold','fontsize',16)
    xticklab = get(gca,'XTickLabel');
    set(gca,'XTickLabel',xticklab,'fontsize',12)
    yticklab = get(gca,'YTickLabel');
    set(gca,'YTickLabel',yticklab,'fontsize',12)
    
    %%%% add drawing of stimulus condition to each subplot
    ax1 = get(h,'position'); % Store handle to axes 1.
    % Create smaller axes in top right, and plot on it
    % Store handle to axes 2 in ax2.
    
    % add background
    width = 0.24;
    height = 0.08;
    ax2 = axes('Position',[pos{c}(1)+0.0675 .755 width height]);
    box on;
    hold(ax2); %
    bg = repmat(min(gratings(:,:,c)./max(gratings(:,:,c),[],'all'),[],'all'),size(gratings,1),size(gratings,1)*3);
    if fix == 0
        imshow(bg);
    elseif fix == 1
        [dfxs,dfys] = drawEllipse(30,30,160,54);
        imshow(bg);
        hold on
        plot(dfxs, dfys, 'w-','linewidth', 2); %plot fix limits
        
    end
    %add stimuli
    if contains(conditions{c},'left')
        ax2 = axes('Position',[pos{c}(1)+loffset .755 .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings(:,:,c)./max(gratings(:,:,c),[],'all');
        imshow(grat);
    elseif contains(conditions{c},'right')
        ax2 = axes('Position',[pos{c}(1)+roffset .755 .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings(:,:,c)./max(gratings(:,:,c),[],'all');
        imshow(grat);
    elseif contains(conditions{c},'both')
        ax2 = axes('Position',[pos{c}(1)+roffset .755 .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings(:,:,c)./max(gratings(:,:,c),[],'all');
        imshow(grat);
        
        ax2 = axes('Position',[pos{c}(1)+loffset .755 .08 .08]);
        box on;
        hold(ax2); %
        grat = gratings(:,:,c)./max(gratings(:,:,c),[],'all');
        imshow(grat);
    end
end
sgtitle(sprintf('Fixation %s',state),'fontweight','bold','fontsize',20)

plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
saveas(gcf,strcat(plotdir, sprintf('electrooculogram_%s_loic09202022.png',savename)));
saveas(gcf,strcat(plotdir, sprintf('electrooculogram_%s_loic09202022.svg',savename)));


%% Compute stdx and stdy and compare across conditions with and without fixation
stdx = nan(10, length(exptDat.exp.conds));
stdy = nan(10, length(exptDat.exp.conds));
for c =1:length(exptDat.exp.conds)

    blockdiff = unique([find(diff(exptDat.exp.longFormConds == c)) length(exptDat.exp.longFormConds)]); %get block onset time and offset time
    for n =1:length(blockdiff)/2
        xbloc = x(blockdiff(2*(n-1)+1):blockdiff(2*n));
        ybloc = y(blockdiff(2*(n-1)+1):blockdiff(2*n));
        xp = xbloc((eyeDat.EyeData.Fixated(blockdiff(2*(n-1)+1):blockdiff(2*n)) == 1))- nanmean(xbloc(eyeDat.EyeData.Fixated(blockdiff(2*(n-1)+1):blockdiff(2*n)) == 1));
        yp = ybloc((eyeDat.EyeData.Fixated(blockdiff(2*(n-1)+1):blockdiff(2*n)) == 1))- nanmean(ybloc(eyeDat.EyeData.Fixated(blockdiff(2*(n-1)+1):blockdiff(2*n)) == 1));
        normX = xp/max(abs(xp),[],'all');
        normY = yp/max(abs(yp),[],'all');
        stdx(n,c) = std(normX);
        stdy(n,c) = std(normY);
    end

end

conds = [{'LeftHor'}, {'RightHor'}, {'BothVert'}, {'BothHor'}];

% ordconds = [{'BothHor'},{'BothVert'},{'LeftHor'},{'RightHor'}];
% ordstdx = [stdx(:,4) stdx(:,3) stdx(:,1) stdx(:,2)]; %need to reorganize in alphabetical order
%%plot stdx
p = jitterCiPlot(conds, stdx,  'Condition', 'Standard deviation (normalized)', 'Standard deviation of horizontal eye position', gratings, fix);


plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
saveas(gcf,strcat(plotdir, sprintf('norm_stdx_%s_loic09202022.png',savename)));
saveas(gcf,strcat(plotdir, sprintf('norm_stdx_%s_loic09202022.svg',savename)));

%%plot stdy
%ordstdy = [stdy(:,4) stdy(:,3) stdy(:,1) stdy(:,2)]; %need to reorganize in alphabetical order
plot = jitterCiPlot(conds, stdy,  'Condition', 'Standard deviation (normalized)', 'Standard deviation of vertical eye position', gratings, fix);
plotdir = strcat('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/anal_plots/');
saveas(gcf,strcat(plotdir, sprintf('norm_stdy_%s_loic09202022.png',savename)));
saveas(gcf,strcat(plotdir, sprintf('norm_stdy_%s_loic09202022.svg',savename)));


%% stats

aovx = anova1(stdx, conds);
aovy = anova1(stdy, conds);


%% Read 1000Hz data

blocks= edfImport('/Users/loicdaumail/Documents/Research_MacBook/Tong_Lab/Projects/Prototype_phantom_FEM/data/Loic_09202022/demo_fix.edf', [1 1 1], '');