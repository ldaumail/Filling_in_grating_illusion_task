function [cleanDat, blink]= removeBlinks(xPos,yPos) %, wind)
%trigdata = groupDat; = data triggered to stimulus onset
%wind = time window to remove around the blink (former version of blink
%removal)
%developped by Loic Daumail on 10/22/2021
%removes blinks from electroooculograms. replaces portions of blinks by
%NaNs 
% 
% xPos = x;
% yPos = y;

cleanDat = nan(2,length(xPos));

xs = xPos;

%1) find the x value of the blink = common value across data. if we
%find 2 chunks with same value ==> it may be a blink value!
if nnz(xPos >= 10000)
    blinkX = unique(xPos(abs(xPos) >= 10000));
end

% 2) eliminate all blinks from all data
blink =ismember(xs,blinkX);
ys = yPos;

xs(blink) = NaN;
ys(isnan(xs)) = NaN;

%replace data
cleanDat(1,:) = xs;
cleanDat(2,:) = ys;

% xcenter =0;
% ycenter =0;
% [scrxs, scrys] = drawSquare(halfscrwidth, scrheight, xcenter,ycenter);
% figure(); plot(x1,y1); hold on; plot(scrxs, scrys);





    %hold on; plot(groupDat.(conditions{j}).xpos(:,12),groupDat.(conditions{j}).ypos(:,12))
    %y = groupDat.(conditions{j}).ypos(1:t,tr);
    %{
halfWind =round(wind/2);
cleanDat = struct();
conditions = fieldnames(trigdata);
for i =1:length(conditions)
    xs = trigdata.(conditions{i}).xpos;
    %1) find the x value of the blink = common value across trials. if we
    %find 2 trials with same value ==> it has to be the blink value!
    flag = 0;
    blinkX = [];
    clear tr1 tr2
    for tr1 =1:size(xs,2)
        xtr1 = xs(:,tr1);
        for tr2 =1:size(xs,2)
            if tr1 ~= tr2
                xtr2 = xs(:,tr2);
                if nnz(ismember(xtr1,xtr2)) >=1
                    blinkX = unique(xtr1(ismember(xtr1,xtr2)));
                    flag =1;
                    break
                    %return
                end
            end
        end
        if flag ==1
            break
        end
    end
    % 2) eliminate all blinks from all trials
    for tr = 1:size(xs,2)
        x1 = xs(:,tr);
        blink =find(diff(ismember(x1,blinkX)));
        y1 = trigdata.(conditions{i}).ypos(:,tr);
        if length(blink) >=2
            for nb =1:2:length(blink) %take all blinks within trial into account
                if mod(length(blink),2) == 0
                    endblink = blink(nb+1);
                else
                    endblink = length(x1);
                end
                if blink(nb)+round((endblink-blink(nb))/2)-halfWind >=1
                    lowB = blink(nb)+round((endblink-blink(nb))/2)-halfWind;
                else
                    lowB =1;
                end
                if blink(nb)+round((endblink-blink(nb))/2)+halfWind <= length(x1)
                    highB = blink(nb)+round((endblink-blink(nb))/2)+halfWind;
                else
                    highB = length(x1);
                end
                x1(lowB:highB) = NaN;
                y1(isnan(x1)) = NaN;
                    % figure(); plot(x1,y1); title(sprintf('%d',tr))
            end
            
        elseif length(blink) == 1
            endblink = blink;
            lowB =1;
            if lowB+round((endblink-lowB)/2)+halfWind <= length(x1)
                highB = lowB+round((endblink-lowB)/2)+halfWind;
            else
                highB = length(x1);
            end
            x1(lowB:highB) = NaN;
            y1(isnan(x1)) = NaN;

        end
    %replace data
    cleanDat.(conditions{i}).xpos(:,tr) = x1;
    cleanDat.(conditions{i}).ypos(:,tr) = y1;
    end
   
end

    %xcenter =0;
    %ycenter =0;
    %[scrxs, scrys] = drawSquare(halfscrwidth, scrheight, xcenter,ycenter);
    %figure(); plot(x1,y1); hold on; plot(scrxs, scrys);
    %hold on; plot(groupDat.(conditions{j}).xpos(:,12),groupDat.(conditions{j}).ypos(:,12))
    %y = groupDat.(conditions{j}).ypos(1:t,tr);
    
end
%}

