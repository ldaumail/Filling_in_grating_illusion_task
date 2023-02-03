function LD_localizer(subject, session, vertOffset, locSize, debug) 
% LD. 01272022
% function LGNedge_loc(locSize,subject,runNum,offset,deviceNumber) 
% flickering checkerboard full field centered on fix
% 136 TRs
% clear all;

%%%% resolution
if debug == 1
	experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop
else
    experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % scanner
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1

Screen('Preference', 'SkipSyncTests', 0);
% input('Hit enter to proceed.');

% take these out for the actual scan!
% offset = 50;
% [keyboardIndices, productNames, ~] = GetKeyboardIndices;
% deviceNumber =keyboardIndices;
% runNum = 1;
% subject = 'test';

experiment.scanNum = input('Scan number :');
experiment.runNum = input('Run number :');

%%%% input this at the beginning of the scan session for 7T
loc.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m
loc.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
loc.whichCLUT = '7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
experiment.subject = subject;
experiment.session = session;

%%%% scales all of the stimuli in DVA to the screensize
loc.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
loc.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;

%%%% set-up rand
% rand('twister', sum(100*clock));
% experiment.rand = rand;
rng(sum(100*clock));
ex.rand = rng;
%%%% files and things
experiment.root = pwd; %'/Users/Sonia/Desktop/ObjectRF/';
experiment.date = datestr(now,30);

%%%% timing
loc.blockLength = 12;            % in seconds
loc.betweenBlocks = 12;          % in seconds
loc.initialFixation = 16;        % in seconds
loc.finalFixation = 16;          % in seconds
% loc.flickerRate = .100;         % in seconds
loc.flickerRate = .050;         % in seconds
loc.numBlocks = 16;

%%%% checkerboard
loc.stim.spatialFreqDeg = 1.5;                                           % cycles per degree
loc.stim.contrast =  1;                                                  % in %, maybe??
loc.stim.stimSizeDeg = locSize;                                           % diam in degrees
loc.stim.annulusRad = 0;                                               % in degrees

%%%% conditions & layout
% loc.fixSizeDeg =  .25;            % in degrees, the size of the biggest white dot in the fixation
loc.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation
loc.littleFixDeg = loc.fixSizeDeg* .5;    % proportion of the fixSizeDeg occupied by the smaller black dot
loc.outerFixPixels = 2;          % in pixels, the black ring around fixation

loc.TRlength = 2;                % in seconds
loc.repsPerRun = 2;              % repetitions of each object type x location
loc.totalTime = loc.initialFixation + (loc.numBlocks * (loc.blockLength + loc.betweenBlocks)) - loc.betweenBlocks + loc.finalFixation;
experiment.allFlips = (0:loc.flickerRate:loc.totalTime);

%%%% screen
loc.backgroundColor = [128 128 128];  % color
loc.fontSize = 20;

%%%%%%%%%%%%%%%%%
% timing model  %
%%%%%%%%%%%%%%%%%

experiment.onSecs = [zeros(1,loc.initialFixation) repmat([ones(1,loc.blockLength) zeros(1,loc.betweenBlocks)],1,loc.numBlocks) zeros(1,loc.finalFixation-loc.betweenBlocks)];
experiment.longFormBlocks = Expand(experiment.onSecs,1/loc.flickerRate,1);
experiment.longFormFlicker = repmat([1 1],1,length(experiment.longFormBlocks)/2);
experiment.whichCheck = repmat([1 1 2 2],1,length(experiment.longFormBlocks)/4);

experiment.conds = repmat([1 2],1,loc.numBlocks/2);
experiment.longFormConds = zeros(1,loc.initialFixation);
for i = 1:loc.numBlocks-1
    experiment.longFormConds = [experiment.longFormConds, repmat(experiment.conds(i),1,loc.blockLength)]; % blocks
    experiment.longFormConds = [experiment.longFormConds, zeros(1,loc.betweenBlocks)]; % inter-block blanks
end
experiment.longFormConds = [experiment.longFormConds, repmat(experiment.conds(end),1,loc.blockLength), zeros(1,loc.finalFixation)]; % the last block
experiment.longFormConds = Expand(experiment.longFormConds, 1/loc.flickerRate,1);

%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
[win, rect]=Screen('OpenWindow',screen,loc.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat);
Screen(win, 'TextSize', loc.fontSize);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
if loc.gammaCorrect > 0
    load(loc.whichCLUT);
    Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
end

%%%% timing optimization
flipInt = Screen('GetFlipInterval',win);
slack = flipInt/2;

%%%% scale the stims for the screen
loc.ppd = pi* rect(3) / (atan(loc.screenWidth/loc.viewingDist/2)) / 360;
loc.stimSize = round(loc.stim.stimSizeDeg*loc.ppd);                 % in degrees, the size of our objects 
loc.innerAnnulus = round(loc.stim.annulusRad*loc.ppd);
loc.fixSize = round(loc.fixSizeDeg*loc.ppd);
loc.littleFix = round(loc.littleFixDeg*loc.ppd);
loc.pixPerCheck = round((loc.ppd/loc.stim.spatialFreqDeg)/2);      % half of the pix/cycle

xc = rect(3)/2; % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+loc.vertOffset;

% loc.checkerboardPix = round(loc.stim.spatialFreqDeg*loc.ppd); % pixels per 1 cycle of checkerboard (B and W square)
% loc.checkerboardSize = round(loc.stimSize/loc.checkerboardPix);
% 
% basicCheck = checkerboard(ceil(loc.checkerboardPix/4),loc.checkerboardSize,loc.checkerboardSize)>.5;

% Define a simple checkerboard with 1 pix/check
%totalCycles = ceil(loc.stim.stimSizeDeg* loc.stim.spatialFreqDeg);
totalCycles = ceil((rect(3)/loc.ppd)* loc.stim.spatialFreqDeg);
checkerboard = repmat(Expand(eye(2),loc.pixPerCheck,loc.pixPerCheck), totalCycles, totalCycles);
checkRect = CenterRectOnPoint([0 0 size(checkerboard,1) size(checkerboard,2)],xc,yc);

% Make the checkerboard into a texure (1 pix per cycle)
checkTex{1} = Screen('MakeTexture',win,255*checkerboard);
checkTex{2} = Screen('MakeTexture',win,255*abs(checkerboard-1)); % inverse

% % We will scale our texure up to x times its current size by defining a
% % larger screen destination rectangle
% [s1, s2] = size(checkerboard);
% dstRect = [0 0 s1 s2] .* loc.pixPerCheck;
% dstRect = CenterRectOnPointd(dstRect, xc, yc);

% make rectangular mask for stim
aperture2=Screen('OpenOffscreenwindow', win, 128);
Screen('FillRect',aperture2, [255 255 255 0], [xc-(1/2)*loc.stimSize yc-loc.stimSize/2 xc+loc.stimSize/2 yc+loc.stimSize/2]);

aperture=Screen('OpenOffscreenwindow', win, 128);
Screen('FillRect',aperture, [255 255 255 0], [xc-(3/2)*loc.stimSize yc-loc.stimSize/2 xc-loc.stimSize/2 yc+loc.stimSize/2]);
Screen('FillRect',aperture, [255 255 255 0], [xc+(1/2)*loc.stimSize yc-loc.stimSize/2 xc+loc.stimSize*(3/2) yc+loc.stimSize/2]);


%%%% initial window - wait for backtick
Screen(win, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(win, 'Flip', 0);

KbTriggerWait(53, deviceNumber);
KbQueueCreate(deviceNumber,responseKeys);
experiment.targetTimes=[];
experiment.responseTimes = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


flipCount = 1;

%%%%%%% START task TASK/FLIPPING
for n = 1:(length(experiment.allFlips)-1)
    KbQueueStart()
    %%%% draw check
    if experiment.longFormBlocks(n) > 0 && experiment.longFormFlicker(n) >0% zeros correspond to blanks, in which case we skip this next section
        
        Screen('DrawTexture', win, checkTex{experiment.whichCheck(n)},[],checkRect);
        if experiment.longFormConds(n) == 1
            Screen('DrawTexture',win,aperture);
        elseif experiment.longFormConds(n) == 2
            Screen('DrawTexture',win,aperture2);
        end
        %Screen('FillOval',win,loc.backgroundColor,[xc-loc.innerAnnulus yc-loc.innerAnnulus xc+loc.innerAnnulus yc+loc.innerAnnulus]);
    end
    
    t = round(rand*100);
    %%%% draw fixation big circle
    if t == 1; 
        colourA = [255 255 255]; 
        colourB = [0 0 0]; 
        experiment.targetTimes = [experiment.targetTimes GetSecs - experiment.startRun];
    else
        colourB = [255 255 255]; 
        colourA = [0 0 0]; 
    end
    Screen('FillOval', win,colourA, [xc-round(loc.fixSize/2+loc.outerFixPixels ) yc-round(loc.fixSize/2+loc.outerFixPixels ) xc+round(loc.fixSize/2+loc.outerFixPixels ) yc+round(loc.fixSize/2+loc.outerFixPixels )]); % black fixation ring
    Screen('FillOval', win,colourB, [xc-round(loc.fixSize/2) yc-round(loc.fixSize/2) xc+round(loc.fixSize/2) yc+round(loc.fixSize/2)]); % white fixation ring
    % little fixation dot
    Screen('FillOval', win,colourA, [xc-round(loc.littleFix/2) yc-round(loc.littleFix/2) xc+round(loc.littleFix/2) yc+round(loc.littleFix/2)]); % black fixation dot

    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if n == 1 [VBLT experiment.startRun FlipT missed] = Screen(win, 'Flip', 0);
        experiment.flipTime(n) = experiment.startRun;
    else
        [VBLT experiment.flipTime(n) FlipT missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(n) - slack);
    end
    
    
    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    if pressed
        experiment.responseTimes = [experiment.responseTimes GetSecs - experiment.startRun];
    end

    KbQueueFlush();
    
end
%%%% to show the very last flip screen for its 200ms
[VBLT experiment.flipTime(n+1) FlipT missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(length(experiment.allFlips)) - slack);

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

experiment.runTime = GetSecs - experiment.startRun;

% analyse task accuracy
responses = experiment.responseTimes; % make copy as we are editing this
responseWindow = 1; % responses later than this count as a miss
for t = 1:size(experiment.targetTimes,2)
    targetTime = experiment.targetTimes(t);
    experiment.hits(t) = 0; % default is a miss with RT of nan
    experiment.RTs(t) = nan;
    for r = 1:size(responses,2)
        if responses(r) > targetTime % if response happened after target
            rt = responses(r)-targetTime;
            if rt < responseWindow; % and if response happened within a second of target
                experiment.hits(t) = 1; % mark a hit
                experiment.RTs(t) = rt; % store the RT
                responses(r) = []; % delete this response so it can't be reused
                break % stop looking for responses to this target
            end
        end
    end
end
experiment.accuracy = (sum(experiment.hits)/size(experiment.targetTimes,2))*100;
experiment.meanRT = nanmean(experiment.RTs);

savedir = fullfile(experiment.root,'data',subject,session,'LGN_loc');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(subject , '_LGN_loc_sn',num2str(experiment.scanNum),'_rn',num2str(experiment.runNum),'_',experiment.date,'.mat'));    
save(savename,'experiment','loc');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;
clear all;

end


