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

%%%% sine wave grating timing

experiment.initialFixation = 6;        % in seconds
experiment.finalFixation = 0;          % in seconds
experiment.stimDur = 2.6;        % in seconds,refers to sine wave grating 2.6 = 1.3*2 = 2 phases
experiment.stimsPerBlock = 5;
experiment.blockLength = experiment.stimDur*experiment.stimsPerBlock;            % in seconds
experiment.betweenBlocks = 12;          % in seconds
experiment.flipsPerSec = 12;           % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
experiment.flipWin = 1/experiment.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%experiment.numBlocks = 12;  % 6 for single and 6 for pair...

%%%% 2D sine wave grating properties
experiment.stim.spatialFreqDeg = 1;                                           % cycles per degree of visual angle
experiment.stim.contrast =  .3;                                                 % in %, maybe??
experiment.stim.orientation = 90;                                                % in degrees
experiment.stim.degFromFix = .6;                                              % in degrees of visual angle
experiment.stim.gaborHDeg = 4;                                                  % in degrees of visual angle
experiment.stim.gaborWDeg = 4; 
experiment.stim.contrastMultiplicator = .075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
experiment.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
experiment.stim.motionRate = 1.3*360 ;                                          % 1.3 cycles per second = 360 deg of phase *1.3 per sec

%%%% conditions & layout
experiment.conds = {'singleTop','singleBottom','doubleIndir'}; %'double-oppdir'
experiment.numConds = length(experiment.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
experiment.condShuffle = Shuffle(repmat([1:experiment.numConds],1,experiment.stimsPerBlock)); %make same number of blocks with each condition, randomize order
experiment.numBlocks = length(experiment.condShuffle);
experiment.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation
experiment.repsPerRun = 2;              % repetitions of each object type x experimentation
experiment.totalTime = experiment.initialFixation + ((experiment.numBlocks-1) * (experiment.blockLength + experiment.betweenBlocks)) + experiment.blockLength + experiment.finalFixation;
experiment.allFlips = (0:experiment.flipWin:experiment.totalTime);

%%%% screen
experiment.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
experiment.fontSize = 20; %26;

%%%% RSVP task setup
experiment.targetProb = .1;              % proportion of trials where the target letters will come up
experiment.firstPossibleTarget = 20;     % no targets in the first X flips
experiment.lastPossibleTarget = 20;      % no targets in the last X flips
experiment.targetLetters = {'J' 'K'};
experiment.distractors = {'A' 'S' 'D' 'F' 'G' 'H' 'L'};
experiment.trialFreq = 0.5;               % duration of fixation trials (seconds) (time it takes to switch from one letter to blank then back to one letter)
flipsPerTrial = experiment.trialFreq/experiment.flipWin;
experiment.trialDur = experiment.trialFreq/2;               % duration in seconds of the letter presentation of each fixation trial (.25s letter ON, .25 letter OFF)
experiment.trialOnFlips = floor(experiment.trialDur/experiment.flipWin); %number of flips the letter is on over the course of 1 trial
%experiment.RSVPrate = 1;               % how fast the letters flip (second) - this likely needs to be hardcoded down below, but is saved here for bookkeeping
%params.cueColor = 0;%[50 50 255];   % letter color
experiment.totalLetters = (experiment.totalTime/experiment.trialFreq); %total number of letters that COULD BE presented during the experiment (though based on prob, much less will be presented)
experiment.responseBack = 3;    % the response is correct if the preceding N letters were the target
%experiment.letterTiming = (0:params.RSVPrate:experiment.totalTime);

%%%% character identification task
experiment.targets = [];
experiment.targetTimes = [];
experiment.responses = [];
experiment.responseTimes=[];
experiment.accuracy = 0;
experiment.meanRT = 0;

taskText = 'characters';

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

experiment.onSecs = [zeros(1,experiment.initialFixation)...
    repmat([ones(1,experiment.blockLength) zeros(1,experiment.betweenBlocks)],1,experiment.numBlocks-1)... %2*ones(1,experiment.blockLength) zeros(1,experiment.betweenBlocks)
    ones(1,experiment.blockLength) zeros(1,experiment.finalFixation)];
experiment.longFormBlocks = Expand(experiment.onSecs,experiment.flipsPerSec,1); %1 when block, 0 when between block
experiment.longFormFlicker = repmat(ones(1,1),1,length(experiment.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(experiment.longFormBlocks)

%set up the timing model for stimulus conditions (single (top, bottom), pair)
%longform condition timing, which aligns with the flicker timing
experiment.longFormConds = zeros(1,experiment.initialFixation);
for i = (1:experiment.numBlocks-1)
    experiment.longFormConds = [experiment.longFormConds, repmat(experiment.condShuffle(i),1,experiment.blockLength)]; % blocks
    experiment.longFormConds = [experiment.longFormConds, zeros(1,experiment.betweenBlocks)]; % inter-block blanks
end
experiment.longFormConds = [experiment.longFormConds, repmat(experiment.condShuffle(experiment.numBlocks),1,experiment.blockLength), zeros(1,experiment.finalFixation)]; % the last block
experiment.longFormConds = Expand(experiment.longFormConds, experiment.flipsPerSec,1);
length(experiment.longFormConds)
% %% create the timing model of stimulus conditions for this particular run
clear i
for i =1:experiment.numConds
    conditions(i).name = experiment.conds(i);
    conditions(i).startTimes = [];
end
counter = experiment.initialFixation; %will assess starting time of each block incrementally, that will be saved in the condition struct
for n=1:experiment.numBlocks
    experiment.startBlock(n) = counter; % timestamp (s) of when each block should start
    conditions(experiment.condShuffle(n)).startTimes = [conditions(experiment.condShuffle(n)).startTimes counter]; % add timestamps to the condition struct
    counter = counter + experiment.blockLength + experiment.betweenBlocks; % and progress the counter
end
%%
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

%HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
[w, rect]=Screen('OpenWindow',screen,experiment.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat);
Screen(w, 'TextSize', experiment.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
% if experiment.gammaCorrect > 0
%     load(experiment.whichCLUT);
%     Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
% end

%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate = Screen('NominalFrameRate',w);
flipTimes = [0:frameInt*frameRate/experiment.flipsPerSec:experiment.stimDur]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
experiment.stim.dphasePerFlip = experiment.stim.motionRate*frameInt * frameRate/experiment.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip

%%%% scale the stims for the screen
experiment.ppd = pi* rect(3) / (atan(experiment.screenWidth/experiment.viewingDist/2)) / 360;
experiment.fixSize = round(experiment.fixSizeDeg*experiment.ppd);
experiment.gaborHeight = round(experiment.stim.gaborHDeg*experiment.ppd);                 % in pixels, the size of our objects
experiment.gaborWidth = round(experiment.stim.gaborWDeg*experiment.ppd);                 % in pixels, the size of our objects

xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+experiment.vertOffset;

%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers
topPhase = randi(360);
bottomPhase = topPhase;
experiment.topWave = nan(length(flipTimes),experiment.gaborHeight,experiment.gaborWidth);
experiment.bottomWave = nan(length(flipTimes),experiment.gaborHeight,experiment.gaborWidth);
experiment.topWaveID = nan(length(flipTimes),1);
experiment.bottomWaveID = nan(length(flipTimes),1);

for f = 1:length(flipTimes)
    
    if f <= length(flipTimes)/2
        
    topPhase = mod(topPhase + experiment.stim.dphasePerFlip,360);
    bottomPhase = mod(bottomPhase + experiment.stim.dphasePerFlip,360);
    %ih in pixels %iw in pixels %spatial freq in cycles per dva
    %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
    %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
    %background color (unused if the grating is not an annulus)
    
    experiment.topWave(f,:,:) = makeSineGrating(experiment.gaborHeight,experiment.gaborWidth,experiment.stim.spatialFreqDeg,...
        experiment.stim.orientation,topPhase,experiment.stim.contrastOffset(1),experiment.stim.contrastMultiplicator,...
        experiment.ppd);
    experiment.bottomWave(f,:,:) = makeSineGrating(experiment.gaborHeight,experiment.gaborWidth,experiment.stim.spatialFreqDeg,...
        experiment.stim.orientation,bottomPhase,experiment.stim.contrastOffset(1),experiment.stim.contrastMultiplicator,...
        experiment.ppd);
%    figure();
 %   imshow(squeeze(topWave(f,:,:)));
%     
    elseif f > length(flipTimes)/2
        topPhase = mod(topPhase - experiment.stim.dphasePerFlip,360);
        bottomPhase = mod(bottomPhase - experiment.stim.dphasePerFlip,360);   
        experiment.topWave(f,:,:) = makeSineGrating(experiment.gaborHeight,experiment.gaborWidth,experiment.stim.spatialFreqDeg,...
            experiment.stim.orientation,topPhase,experiment.stim.contrastOffset(1),experiment.stim.contrastMultiplicator,...
            experiment.ppd);
        experiment.bottomWave(f,:,:) = makeSineGrating(experiment.gaborHeight,experiment.gaborWidth,experiment.stim.spatialFreqDeg,...
            experiment.stim.orientation,bottomPhase,experiment.stim.contrastOffset(1),experiment.stim.contrastMultiplicator,...
            experiment.ppd);
    end
    experiment.topWaveID(f) = Screen('MakeTexture', w, squeeze(experiment.topWave(f,:,:)));
    experiment.bottomWaveID(f) = Screen('MakeTexture', w, squeeze(experiment.bottomWave(f,:,:)));
            
end
%% extend stimulus matrix to include the same total number of flips as the whole experiment
experiment.topWaveID = repmat(experiment.topWaveID,1,length(experiment.longFormFlicker));
experiment.bottomWaveID = repmat(experiment.bottomWaveID,1,length(experiment.longFormFlicker));

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+experiment.vertOffset;

xtop = rect(3)/2 - (4+2)*experiment.ppd;%rect(3)/4; % = stimulus center located 4 degrees horizontal from the center
ytop = rect(4)/2 - experiment.gaborHeight;
experiment.topRect =  CenterRectOnPoint([0 0 experiment.gaborWidth experiment.gaborHeight],xtop,ytop);

xbottom = rect(3)/2 - (4+2)*experiment.ppd; %rect(3)/4; = stimulus center located 4 degrees horizontal from the center
ybottom = rect(4)/2 + experiment.gaborHeight;%rect(4)/4*3;
experiment.bottomRect =  CenterRectOnPoint([0 0 experiment.gaborWidth experiment.gaborHeight],xbottom,ybottom);


%% %%%%%%%%%%%%%%%%%%%%%%
   % Letter task set-up %
   %%%%%%%%%%%%%%%%%%%%%%

%%%% find the target positioning for this run (index of repeated letters)
experiment.numTargets = length(find(rand(1,experiment.totalLetters)<experiment.targetProb)); %draw targets from uniform distribution between 0 nd 1, only keep the density we need (0.1) of them
targetInd = zeros(1,experiment.totalLetters+1); % this has to go to +1, but it's only for the repeat check
if experiment.numTargets > 0
    while 1
        %check previous (-1,-2) and following indices (+1,+2) if they might constitute a
        %repeated letter
        maybeRep = experiment.firstPossibleTarget+Randi(experiment.totalLetters-experiment.firstPossibleTarget-experiment.lastPossibleTarget);    % a possible index for a target letter, picked randomly
        if targetInd(maybeRep-1) == 0 && targetInd(maybeRep+1) == 0 && targetInd(maybeRep-2) == 0 && targetInd(maybeRep+2) == 0 % make sure the previous  TWO and following TWO letters weren't already a repeat
            targetInd(maybeRep) = 1; %record this index if it is not a repeat
        end
        if sum(targetInd) == experiment.numTargets %if the sum of target letters reach the number of targets, stop the check
            break
        end
    end
end
targetInd = targetInd(1:experiment.totalLetters); % trim this back for sanity


%%%% fill the sequence out with actual letters
experiment.letterSequence = [];
for k = 1:length(targetInd)
    if targetInd(k) == 1 % targets - we know these don't repeat
        experiment.letterSequence{k} = experiment.targetLetters{randi(length(experiment.targetLetters))}; %pick a random target
    else % distractors - we need to make sure these don't repeat
        if k > 1 % the first letter can be whatever
            while 1
                maybeLetter = experiment.distractors{randi(length(experiment.distractors))};
                if strcmp(maybeLetter, experiment.letterSequence{k-1}) == 0
                    experiment.letterSequence{k} = maybeLetter;
                    break
                end
            end
        else %% (the first letter can be whatever)
            experiment.letterSequence{k} = experiment.distractors{randi(length(experiment.distractors))};
        end
    end
end

% %%%% scale up the task based on flips
experiment.letterSequence = Expand(experiment.letterSequence,flipsPerTrial,1);%experiment.flipsPerSec,1); 


%% %%%% initial window - wait for backtick
Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(w, 'Flip', 0);

%%%%%%%%%%%%%%%%%% Response listening %%%%%%%%%%%%%%%%%%%%%%%%
KbTriggerWait(53, deviceNumber);
KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n=0;
count = 1;
%%%%%%% START task TASK/FLIPPING
% [experiment.totalTime, experiment.allFlips,n]
gray = repmat(min(min(squeeze(experiment.topWave(1,:,:)),[],1)), [1,3]);
Screen('FillRect', w, gray);
while n+1 < length(experiment.allFlips)
    
    [experiment.longFormBlocks(n+1),experiment.longFormFlicker(n+1)]
    thisCond = experiment.longFormConds(n+1);
    tic
    KbQueueStart();
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 0 
        [VBLT, experiment.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        experiment.flipTime(n+1) = experiment.startRun;
    else
        [VBLT, experiment.flipTime(n+1), FlipT, missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(n+1) - slack);
    end
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if experiment.longFormBlocks(n+1) == 1 && experiment.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        
        % draw & increment stims
        if strcmp(conditions(thisCond).name, 'doubleIndir')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            % top stim
            Screen('DrawTexture', w, experiment.topWaveID(n+1),[],experiment.topRect);
            % bottom stim
            Screen('DrawTexture', w, experiment.bottomWaveID(n+1), [], experiment.bottomRect);
        elseif strcmp(conditions(thisCond).name, 'singleTop')
            % top stim
            Screen('DrawTexture', w, experiment.topWaveID(n+1),[],experiment.topRect);
        elseif strcmp(conditions(thisCond).name, 'singleBottom')
            % bottom stim
            Screen('DrawTexture', w, experiment.bottomWaveID(n+1), [], experiment.bottomRect);
        end
    end
    
    % select new character if starting new trial
    if mod(n, flipsPerTrial) == 0

%         % draw character
%         l = randperm(numel(letters));
%         fixChar = letters(l(1));
          fixChar = experiment.letterSequence(n+1);

       if strcmp(fixChar,'J') == 1
            experiment.targets = [experiment.targets, 1];
            experiment.targetTimes = [experiment.targetTimes, GetSecs - experiment.startRun];
        elseif strcmp(fixChar,'K') == 1
            experiment.targets = [experiment.targets, 2];
            experiment.targetTimes = [experiment.targetTimes, GetSecs - experiment.startRun];
       end


    end
    
    %%%% draw fixation letter in fixation circle
    
    if mod(n, flipsPerTrial) <= experiment.trialOnFlips
        Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]);%white fixation solid circle
        %DrawFormattedText(w, fixChar, 'center', 8+vertOffset+rect(4)/2,0); %either text function works
        %Screen('DrawText', w, fixChar, -5+rect(3)/2, -10+vertOffset+rect(4)/2,[0 0 0]);
        [width,height] = RectSize(Screen('TextBounds',w,fixChar{:}));
        Screen('DrawText', w, fixChar{:}, xc-width/2, yc-height/2,[0 0 0]);
     else
        Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]);%white fixation solid circle
   
    end

%     n

    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    
    %%%% character identification
    if (pressed == 1) && ((firstPress(KbName('1!')) > 0) || (firstPress(KbName('2@')) > 0))
        if firstPress(KbName('1!')) > 0
            experiment.responses = [experiment.responses, 1];
            experiment.responseTimes = [experiment.responseTimes, firstPress(KbName('1!')) - experiment.startRun];
        elseif firstPress(KbName('2@')) > 0
            experiment.responses = [experiment.responses, 2];
            experiment.responseTimes = [experiment.responseTimes, firstPress(KbName('2@')) - experiment.startRun];
        end
    end

    %%%% refresh queue for next character
    KbQueueFlush();
    n = n+1;
    toc
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
            if rt < responseWindow % and if response happened within a second of target
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

savedir = fullfile(experiment.root,'data',subject,session,'fillingin_rsvp_v1');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(subject , '_fillingin_rsvp_v1_sn',num2str(experiment.scanNum),'_rn',num2str(experiment.runNum),'_',experiment.date,'.mat'));
save(savename,'experiment');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;
