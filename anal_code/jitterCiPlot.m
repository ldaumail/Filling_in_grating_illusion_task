function [p] = jitterCiPlot(xconds, yvar,  xlab, ylab, tit, gratings, fix,min);
%function to plot multiple sample groups next to each other with jitter point for
%each data point and mean +95% CI each group
%Loic Daumail 09/21/2022
draw = 1;
% xconds = conds;
% yvar = stdx;
% xlab ='Condition';
% ylab ='Standard deviation (normalized)';
% tit = 'Standard deviation of horizontal eye position';

if nargin > 5
    draw = 1;
end
nlines = 7;
cmaps = struct();
cmaps(1).map =cbrewer2('OrRd', nlines);
cmaps(2).map =cbrewer2('Blues', nlines);
cmap = flip(cmaps(2).map) ;

p = figure('Position',[100 100 900 500]);
%plot the data points
% and the mean ±95%CI
mYvar = nanmean(yvar,1);

%95% CI
ci_high = mean(yvar,1) + 1.96*std(yvar,[],1)/sqrt(size(yvar,1));
ci_low = mean(yvar,1) - 1.96*std(yvar,[],1)/sqrt(size(yvar,1));

for c =1:length(xconds)
    x = c*ones(1, length(yvar(:,c)));
    scatter(x, yvar(:,c),40,'MarkerFaceColor',cmap(4,:), 'MarkerEdgeColor',cmap(4,:),'LineWidth',1.5);
    hold on;
    scatter(c*1,mYvar(c),60,'MarkerFaceColor',cmaps(1).map(4,:), 'MarkerEdgeColor',cmaps(1).map(4,:),'LineWidth',1.5); %mean
    hold on
    line([c*1 c*1], [ci_low(c) ci_high(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI vertical
    hold on
    line([c*1-0.1 c*1+0.1], [ci_low(c) ci_low(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI whiskers
    hold on
    line([c*1-0.1 c*1+0.1], [ci_high(c) ci_high(c)], 'Color', cmaps(1).map(4,:), 'LineWidth', 2); %95%CI whiskers
    hold on
end
% Set up axes.
xlim([0, length(xconds)+1]);
ylim([0, max(yvar,[],'all')+0.1]);
ax = gca;
ax.XTick = [1:length(xconds)];
ax.XTickLabels = xconds;
xlabel(xlab,'fontweight','bold','fontsize',16)
ylabel(ylab,'fontweight','bold','fontsize',16)
yticklab = get(gca,'YTickLabel');
set(gca,'YTickLabel',yticklab,'fontsize',12)
title(tit)

loffset = 0.17;
roffset = 0.25;
bgLocs = [0.21, 0.365, 0.516, 0.671];
if draw == 1
    for c =1:length(xconds)
        %%%% add drawing of stimulus condition to each subplot % could remove
        %%%% this part for more general plots
        %ax1 = get(h,'position'); % Store handle to axes 1.
        % Create smaller axes in top right, and plot on it
        % Store handle to axes 2 in ax2.
        
        % add background
        width = 0.18;
        height = 0.06;
        ax2 = axes('Position',[bgLocs(c) .15 width height]);
        box on;
        hold(ax2); %
        if min
            bg = repmat(min(gratings.squ(:,:)./max(gratings.squ(:,:),[],'all'),[],'all'),size(gratings.squ,1)*4,size(gratings.rect,2)*1.5);
        else
            bg = repmat(mean(gratings.squ(:,:)./max(gratings.squ(:,:),[],'all'),'all'),size(gratings.squ,1)*4,size(gratings.rect,2)*1.5);
        end
        if fix == 0
            imshow(bg);
        elseif fix == 1
            [dfxs,dfys] = drawEllipse(30,30,157,54);
            imshow(bg);
            hold on
            plot(dfxs, dfys, 'w-','linewidth', 2); %plot fix limits
            hold on
        end
        %add stimuli
        if contains(xconds{c},'Squ')
            ax2 = axes('Position',[bgLocs(c)+0.013 .15 .06 .06]);
 
            %ax2 = axes('Position',[bgLocs(c)+0.082 .15 .06 .06]);
            %ax2 = axes('Position',[bgLocs(c)+0.082 .15 .06 .06]);
            box on;
            hold(ax2); %
            grat = gratings.squ(:,:)./max(gratings.squ(:,:),[],'all');
            imshow(grat);
            
            ax2 = axes('Position',[bgLocs(c)+0.013 .15 .06 .06]);
            box on;
            hold(ax2); %
            grat = gratings.squ(:,:)./max(gratings.squ(:,:),[],'all');
            imshow(grat);
            hold on
        elseif contains(xconds{c},'Rect')
            ax2 = axes('Position',[bgLocs(c)+0.082 .15 .06 .06]);
            box on;
            hold(ax2); %
            grat = gratings.rect(:,:)./max(gratings.rect(:,:),[],'all');
            imshow(grat);
            
            ax2 = axes('Position',[bgLocs(c)+0.013 .15 .06 .06]);
            box on;
            hold(ax2); %
            grat = gratings.rect(:,:)./max(gratings.rect(:,:),[],'all');
            imshow(grat);
            hold on
        end
    end
end

%The  commented code below uses the toolbox gramm and is not flexible in
%order to adjust font sizes/fontweights etc
%!!! Very important: in order to use this function, the conditions/labels
%of the data must be organized in alphabetical order and data must be
%organized the same way. Otherwise the results are not reliable

% linConds = [];
% for i =1:length(xconds) 
% linConds = [linConds repmat(xconds(i), 1,10)];
% end
% liny = reshape(yvar,1,size(yvar,1)*size(yvar,2));
% ci_high = mean(yvar,1) + 1.96*std(yvar,[],1)/sqrt(size(yvar,1));
% ci_low = mean(yvar,1) - 1.96*std(yvar,[],1)/sqrt(size(yvar,1));
% 
% nlines = 7;
% cmaps = struct();
% cmaps(1).map =cbrewer2('OrRd', nlines);
% cmaps(2).map =cbrewer2('Blues', nlines);
% cmaps(3).map =cbrewer2('Greens', nlines);
% cmap = flip(cmaps(2).map) ;
% %colormap(cmap);
% clear g
% %jitter
% plot = figure('Position',[100 100 800 450]);
% g(1,1) = gramm('x',categorical(linConds),'y',liny);
% g(1,1).geom_jitter('width',0.2,'height',0);
% g(1,1).set_names('y',ylabel,'x',xlabel);
% g(1,1).axe_property('ylim',[0 0.5]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
% g(1,1).set_color_options('map',cmap(4,:));
% g(1,1).set_title(title,'FontSize',18);
% 
% %point bar plot
% 
% g(1,1).update('x',categorical(unique(linConds)),'y',mean(yvar,1),...
%     'ymin',ci_low,'ymax',ci_high);
% g(1,1).set_color_options('map',cmaps(1).map(4,:));
% g(1,1).geom_point('dodge',0.2);
% g(1,1).geom_interval('geom','errorbar','dodge',0.2,'width',0.8);
% g(1,1).axe_property('ylim',[0 0.5]); %We have to set y scale manually, as the automatic scaling from the first plot was forgotten
% set([g.results.set_names(g.results.color==4).line_handle],'LineWidth',3);
% g.draw();