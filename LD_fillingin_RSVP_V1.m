%function figureGround_loc_v5(subject, session, vertOffset, debug) 
subject = 1;
session = 1;
vertOffset = 0;
%%%% resolution
% %if debug == 1
%     experiment.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
%     experiment.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;
% 	experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop
% else
    experiment.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    experiment.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;
    experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % scanner
%end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2

Screen('Preference', 'SkipSyncTests', 0);

experiment.scanNum = input('Scan number :');
experiment.runNum = input('Run number :');
experiment.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m
experiment.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
experiment.whichCLUT = '7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
experiment.subject = subject;
experiment.session = session;


%%%% set-up rand
 rand('twister', sum(100*clock));
 experiment.rand = rand;

%rng(sum(100*clock));
%experiment.rand = rng;
%%%% files and things
experiment.root = pwd;
experiment.date = datestr(now,30);

%%%% timing
experiment.blockLength = 12;            % in seconds
experiment.betweenBlocks = 12;          % in seconds
experiment.initialFixation = 12;        % in seconds
experiment.finalFixation = 0;          % in seconds
experiment.flickerRate = .05;         % in seconds then actually in 1 sec the stimuli will change 10 times instead of 20 time due to the way the code was set up
experiment.numBlocks = 12;  % 6 for single and 6 for pair...
experiment.trialFreq = 1;               % frequency of fixation trials (seconds)
experiment.trialDur = .4;               % duration in seconds of each fixation trial
experiment.postTargetWindow = 1;           % duration in seconds of target-free interval after a target
flipsPerTrial = experiment.trialFreq/experiment.flickerRate;
trialOnFlips = experiment.trialDur/experiment.flickerRate;

%%%% checkerboard
experiment.stim.spatialFreqDeg = 1.5;   % cycles per degree
experiment.stim.contrast =  1;          % in %, maybe??
experiment.stim.stimSizeDeg = 14;       % keep this at 14 as it is spatially restricted later
experiment.stim.degFromFix = 1;         % distance from fixation in degrees
experiment.stim.centerRad = 2;          % radius in degrees
experiment.stim.surroundRad = 6;          % radius in degrees
experiment.stim.annulusRad = experiment.stim.centerRad;         % in degrees

%%%% conditions & layout
experiment.fixSizeDeg =  .4;            % in degrees, the size of the biggest white dot in the fixation
experiment.littleFixDeg = .25;    % proportion of the fixSizeDeg occupied by the smaller black dot
experiment.outerFixPixels = 2;          % in pixels, the black ring around fixation
experiment.TRlength = 2;                % in seconds
experiment.repsPerRun = 2;              % repetitions of each object type x experimentation
experiment.totalTime = experiment.initialFixation + (experiment.numBlocks * (experiment.blockLength + experiment.betweenBlocks)) + experiment.finalFixation;
experiment.allFlips = (0:experiment.flickerRate:experiment.totalTime);

%%%% screen
experiment.backgroundColor = [127 127 127];  % color
experiment.fontSize = 26;

%%%%%%%%%%%%%%%%%
% timing model  %
%%%%%%%%%%%%%%%%%

experiment.onSecs = [zeros(1,experiment.initialFixation)...
    repmat([ones(1,experiment.blockLength) zeros(1,experiment.betweenBlocks) 2*ones(1,experiment.blockLength) zeros(1,experiment.betweenBlocks)],1,experiment.numBlocks/2)...
    zeros(1,experiment.finalFixation)];
experiment.longFormBlocks = Expand(experiment.onSecs,1/experiment.flickerRate,1);
experiment.longFormFlicker = repmat([1 1],1,length(experiment.longFormBlocks)/2);
experiment.whichCheck = repmat([1 1 2 2],1,length(experiment.longFormBlocks)/4);
length(experiment.longFormBlocks)

%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
[win, rect]=Screen('OpenWindow',screen,experiment.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat);
Screen(win, 'TextSize', experiment.fontSize);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
% if experiment.gammaCorrect > 0
%     load(experiment.whichCLUT);
%     Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
% end

%%%% timing optimization
flipInt = Screen('GetFlipInterval',win);
slack = flipInt/2;

%%%% scale the stims for the screen
%experiment.ppd = (rect(3)/experiment.screenWidth) * 2 * experiment.viewingDist * tand(1/2);
experiment.ppd = pi* rect(3) / (atan(experiment.screenWidth/experiment.viewingDist/2)) / 360;
experiment.stimSize = round(experiment.stim.stimSizeDeg*experiment.ppd);                 % in degrees, the size of our objects
experiment.innerAnnulus = round(experiment.stim.annulusRad*experiment.ppd);
experiment.fixSize = round(experiment.fixSizeDeg*experiment.ppd);
experiment.littleFix = round(experiment.littleFixDeg*experiment.ppd);
experiment.pixPerCheck = round((experiment.ppd/experiment.stim.spatialFreqDeg)/2);      % half of the pix/cycle

xc = rect(3)/2; % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+experiment.vertOffset;

% experiment.checkerboardPix = round(experiment.stim.spatialFreqDeg*experiment.ppd); % pixels per 1 cycle of checkerboard (B and W square)
% experiment.checkerboardSize = round(experiment.stimSize/experiment.checkerboardPix);
% 
% basicCheck = checkerboard(ceil(experiment.checkerboardPix/4),experiment.checkerboardSize,experiment.checkerboardSize)>.5;

% Define a simple checkerboard with 1 pix/check
totalCycles = ceil(experiment.stim.stimSizeDeg * experiment.stim.spatialFreqDeg);
checkerboard = repmat(Expand(eye(2),experiment.pixPerCheck,experiment.pixPerCheck), totalCycles, totalCycles);
checkRect = CenterRectOnPoint([0 0 size(checkerboard,1) size(checkerboard,2)],xc,yc);
% Make the checkerboard into a texure (1 pix per cycle)
checkTex{1} = Screen('MakeTexture',win,255*checkerboard);
checkTex{2} = Screen('MakeTexture',win,255*abs(checkerboard-1)); % inverse

% make left and right circular masks for stim
apertureCenter=Screen('OpenOffscreenwindow', win, 128);
Screen('FillOval',apertureCenter,[255 255 255 0],[xc-experiment.ppd*(experiment.stim.centerRad) yc-experiment.ppd*experiment.stim.centerRad xc+experiment.ppd*(experiment.stim.centerRad) yc+experiment.ppd*experiment.stim.centerRad]);
apertureSurround=Screen('OpenOffscreenwindow', win, 128);
Screen('FillOval',apertureSurround,[255 255 255 0],[xc-experiment.ppd*(experiment.stim.surroundRad) yc-experiment.ppd*experiment.stim.surroundRad xc+experiment.ppd*(experiment.stim.surroundRad) yc+experiment.ppd*experiment.stim.surroundRad]);

%%%% initial window - wait for backtick
Screen(win, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(win, 'Flip', 0);

KbTriggerWait(53, deviceNumber);
KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment.targetTimes=[];
experiment.targets = {};
experiment.responseTimes=[];
experiment.responses = {};
experiment.allColors = {};
Colors = 'kkrg';
targetColors = ['rg'];
lastTargetCounter = 0;
lastColor = 0;
n=0;
count = 1;
%%%%%%% START task TASK/FLIPPING
% [experiment.totalTime, experiment.allFlips,n]
while n+1 < length(experiment.allFlips)
    [experiment.longFormBlocks(n+1),experiment.longFormFlicker(n+1)]

    KbQueueStart();
    
    %%%% draw check
    if experiment.longFormBlocks(n+1) == 1 && experiment.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        Screen('DrawTexture', win, checkTex{experiment.whichCheck(n+1)},[],checkRect);
        Screen('DrawTexture',win,apertureCenter);
%         Screen('FillOval',win,experiment.backgroundColor,[xc-experiment.innerAnnulus yc-experiment.innerAnnulus xc+experiment.innerAnnulus yc+experiment.innerAnnulus]);
    elseif experiment.longFormBlocks(n+1) == 2 && experiment.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        Screen('DrawTexture', win, checkTex{experiment.whichCheck(n+1)},[],checkRect);
        Screen('DrawTexture',win,apertureSurround);
        Screen('FillOval',win,experiment.backgroundColor,[xc-experiment.innerAnnulus yc-experiment.innerAnnulus xc+experiment.innerAnnulus yc+experiment.innerAnnulus]);
    end
    
    % select new character if starting new trial
    if mod(n, flipsPerTrial) == 0
        l = randperm(4);
        fixColor = Colors(l(1));
        experiment.allColors = [experiment.allColors, fixColor];
        while any(ismember(fixColor, targetColors)) && lastTargetCounter < experiment.postTargetWindow/experiment.trialDur
            experiment.allChars = [experiment.allColors, 'change'];
            l = randperm(4);
            fixColor = Colors(l(1));
            experiment.allColors = [experiment.allColors, fixColor];
        end
        lastTargetCounter = lastTargetCounter + 1;
        lastColor = fixColor;
        
        if fixColor == 'r';
            experiment.targets = [experiment.targets, 'r'];
            experiment.targetTimes = [experiment.targetTimes, GetSecs - experiment.startRun];
            lastTargetCounter = 0;
        elseif fixColor == 'g';
            experiment.targets = [experiment.targets, 'g'];
            experiment.targetTimes = [experiment.targetTimes, GetSecs - experiment.startRun];
            lastTargetCounter = 0;
        end
    end

    %%%% draw fixation circle
    
    % draw character for 3 flips of the 5 for each trials
    if mod(n, flipsPerTrial) < trialOnFlips & fixColor == 'r'
    Screen('FillOval', win,[128 128 128], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % black fixation ring
    Screen('FillOval', win,[255 0 0], [xc-round(experiment.littleFix/2) yc-round(experiment.littleFix/2) xc+round(experiment.littleFix/2) yc+round(experiment.littleFix/2)]); % black fixation ring
    elseif mod(n, flipsPerTrial) < trialOnFlips & fixColor == 'g'
    Screen('FillOval', win,[128 128 128], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % black fixation ring
    Screen('FillOval', win,[0 200 0], [xc-round(experiment.littleFix/2) yc-round(experiment.littleFix/2) xc+round(experiment.littleFix/2) yc+round(experiment.littleFix/2)]); % black fixation ring
    else
    Screen('FillOval', win,[128 128 128], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % black fixation ring
    Screen('FillOval', win,[0 0 255], [xc-round(experiment.littleFix/2) yc-round(experiment.littleFix/2) xc+round(experiment.littleFix/2) yc+round(experiment.littleFix/2)]); % black fixation ring
    end
    %Screen('DrawText', win, fixChar, -10+rect(3)/2, -14+vertOffset+rect(4)/2,[0 0 0]);
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if n == 0 [VBLT, experiment.startRun, FlipT, missed] = Screen(win, 'Flip', 0);
        experiment.flipTime(n+1) = experiment.startRun;
    else
        [VBLT, experiment.flipTime(n+1), FlipT, missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(n+1) - slack);
    end
%     n
%     fixColor
    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    if pressed
        firstPress(firstPress == 0) = nan;
        [RT,key] = min(firstPress);
        KeyName = KbName(key);
        experiment.responses = [experiment.responses, targetColors(str2num([KeyName(1)]))];
        experiment.responseTimes = [experiment.responseTimes, RT - experiment.startRun];
    end


    KbQueueFlush();
    n = n+1;
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

experiment.runTime = GetSecs - experiment.startRun;

% analyse task accuracy
responseTimes = experiment.responseTimes; % make copies as we are editing these
responses = experiment.responses;
responseWindow = 1.5; % responses later than this count as a miss
for t = 1:size(experiment.targetTimes,2)
    targetTime = experiment.targetTimes(t);
    target = experiment.targets(t);
    experiment.hits(t) = 0; % default is a miss with RT of nan
    experiment.RTs(t) = nan;
    for r = 1:size(responseTimes,2)
        if ismember(responses(r),  target) && responseTimes(r) > targetTime % if response is correct and happened after target
            rt = responseTimes(r)-targetTime;
            if rt < responseWindow; % and if response happened within a second of target
                experiment.hits(t) = 1; % mark a hit
                experiment.RTs(t) = rt; % store the RT
                responseTimes(r) = []; % delete this response so it can't be reused
                responses(r) = [];
                break % stop looking for responses to this target
            end
        end
    end
end
experiment.accuracy = (sum(experiment.hits)/size(experiment.targetTimes,2))*100;
experiment.meanRT = nanmean(experiment.RTs);

savedir = fullfile(experiment.root,'data',subject,session,'figureGround_loc_v5');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(subject , '_figureGround_loc_v5_sn',num2str(experiment.scanNum),'_rn',num2str(experiment.runNum),'_',experiment.date,'.mat'));    
save(savename,'experiment');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;
