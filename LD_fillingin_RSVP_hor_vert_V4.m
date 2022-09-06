%function figureGround_loc_v5(subject, session, vertOffset, debug) 
subject = 1;
session = 1;
debug = 0;
vertOffset = 0;


global EyeData currPosID rect w xc yc
%%%% resolution
if debug == 1
    exp.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    exp.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;
	exp.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop
 else
    exp.screenWidth = 16;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    exp.viewingDist = 23;             % in cm; %23 in eye tracking room 425 3Tb/office=43, miniHelm=57;
    exp.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % scanner
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2

Screen('Preference', 'SkipSyncTests', 0);

exp.scanNum = input('Scan number :');
exp.runNum = input('Run number :');
exp.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m
exp.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
exp.whichCLUT = '7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
exp.subject = subject;
exp.session = session;


%%%% set-up rand
 rand('twister', sum(100*clock));
 exp.rand = rand;

%rng(sum(100*clock));
%e.rand = rng;
%%%% files and things
exp.root = pwd;
exp.date = datestr(now,30);

%%%% sine wave grating timing

exp.initialFixation = 6;        % in seconds
exp.finalFixation = 0;          % in seconds
exp.stimDur = 1.6*2;        % in seconds. 1.6 sec refers to sine wave grating 1.6 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
exp.stimsPerBlock = 5;
exp.blockLength = exp.stimDur*exp.stimsPerBlock;            % in seconds
exp.betweenBlocks = 16;          % in seconds
exp.flipsPerSec = 12;           % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
exp.flipWin = 1/exp.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%e.numBlocks = 12;  % 6 for single and 6 for pair...

%%%% 2D sine wave grating properties
exp.stim.spatialFreqDeg = 0.286;  %1                                         % cycles per degree of visual angle
exp.stim.contrast = 0.3 ;                                                 % in %, maybe??
exp.stim.orientation = [90 180];                                                % in degrees
exp.stim.degFromFix = .6;                                              % in degrees of visual angle
exp.stim.gaborHDeg = 7;                                                  % in degrees of visual angle
exp.stim.gaborWDeg = 7; 
exp.stim.contrastMultiplicator = .075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
exp.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
exp.stim.motionRate = 1.13*360 ;  %1.3                                        % 1.3 cycles per second = 360 deg of phase *1.3 per sec

%%%% conditions & layout
exp.conds = {'vertSingleTop','horSingleTop','vertSingleBottom','horSingleBottom','vertDoubleIndir','horDoubleIndir'}; %'double-oppdir'
exp.numConds = length(exp.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
exp.condShuffle = Shuffle(repmat([1:exp.numConds],1,1)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order
exp.numBlocks = length(exp.condShuffle)*2;
exp.fixSizeDeg =  .6;            % in degrees, the size of the biggest white dot in the fixation
exp.repsPerRun = 2;              % repetitions of each object type x eation
exp.totalTime = exp.initialFixation + ((exp.numBlocks-1) * (exp.blockLength + exp.betweenBlocks)) + exp.blockLength + exp.finalFixation;
exp.allFlips = (0:exp.flipWin:exp.totalTime);

%%%% screen
exp.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
exp.fontSize = 12; %26;

%%%% RSVP task setup
exp.targetProb = .1;              % proportion of trials where the target letters will come up
exp.firstPossibleTarget = 20;     % no targets in the first X flips
exp.lastPossibleTarget = 20;      % no targets in the last X flips
exp.targetLetters = {'J' 'K'};
exp.distractors = {'A' 'S' 'D' 'F' 'G' 'H' 'L'};
exp.trialFreq = 0.5;               % duration of fixation trials (seconds) (time it takes to switch from one letter to blank then back to one letter)
flipsPerTrial = exp.trialFreq/exp.flipWin;
exp.trialDur = exp.trialFreq/2;               % duration in seconds of the letter presentation of each fixation trial (.25s letter ON, .25 letter OFF)
exp.trialOnFlips = floor(exp.trialDur/exp.flipWin); %number of flips the letter is on over the course of 1 trial
%e.cueColor = 0;%[50 50 255];   % letter color
exp.totalLetters = (exp.totalTime/exp.trialFreq); %total number of letters that COULD BE presented during the e (though based on prob, much less will be presented)
exp.responseBack = 3;    % the response is correct if the preceding N letters were the target
%e.letterTiming = (0:e.RSVPrate:e.totalTime);

%%%% character identification task
exp.targets = [];
exp.targetTimes = [];
exp.responses = [];
exp.responseTimes=[];
exp.accuracy = 0;
exp.meanRT = 0;

taskText = 'characters';



%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

exp.onSecs = [zeros(1,exp.initialFixation)...
    repmat([ones(1,exp.blockLength) zeros(1,exp.betweenBlocks)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
    ones(1,exp.blockLength) zeros(1,exp.finalFixation)];
exp.longFormBlocks = Expand(exp.onSecs,exp.flipsPerSec,1); %1 when block, 0 when between block
exp.longFormFlicker = repmat(ones(1,1),1,length(exp.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(exp.longFormBlocks)


%set up the timing model for stimulus pairing conditions (single (top,
%bottom), pair) and stimulus orientation
%longform condition timing, which aligns with the flicker timing
exp.longFormConds = zeros(1,exp.initialFixation);
for n = 1:exp.numBlocks/length(exp.condShuffle)
    if n < exp.numBlocks/length(exp.condShuffle)
        len = length(exp.condShuffle);
    else
        len = length(exp.condShuffle)-1;
    end
    for i = (1:len)
        exp.longFormConds = [exp.longFormConds, repmat(exp.condShuffle(i),1,exp.blockLength)]; % blocks
        exp.longFormConds = [exp.longFormConds, zeros(1,exp.betweenBlocks)]; % inter-block blanks
    end
end
exp.longFormConds = [exp.longFormConds, repmat(exp.condShuffle(end),1,exp.blockLength), zeros(1,exp.finalFixation)]; % the last block
exp.longFormConds = Expand(exp.longFormConds, exp.flipsPerSec,1);
length(exp.longFormConds)
% %% create the timing model of stimulus conditions for this particular run
clear i
for i =1:exp.numConds
    conditions(i).name = exp.conds(i);
    conditions(i).startTimes = [];
end
% counter = e.initialFixation; %will assess starting time of each block incrementally, that will be saved in the condition struct
% for n=1:e.numBlocks
%     e.startBlock(n) = counter; % timestamp (s) of when each block should start
%     conditions(e.condShuffle(n)).startTimes = [conditions(e.condShuffle(n)).startTimes counter]; % add timestamps to the condition struct
%     counter = counter + e.blockLength + e.betweenBlocks; % and progress the counter
% end
%%
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

%HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
[w, rect]=Screen('OpenWindow',screen,exp.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
%[w, rect]=Screen('OpenWindow',screen,exp.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
Screen(w, 'TextSize', exp.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
% if e.gammaCorrect > 0
%     load(e.whichCLUT);
%     Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
% end

%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate = Screen('NominalFrameRate',w);
flipTimes = [0:frameInt*frameRate/exp.flipsPerSec:exp.stimDur]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
exp.stim.dphasePerFlip = exp.stim.motionRate*frameInt * frameRate/exp.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip

%%%% scale the stims for the screen
exp.ppd = pi* rect(3) / (atan(exp.screenWidth/exp.viewingDist/2)) / 360;
exp.fixSize = round(exp.fixSizeDeg*exp.ppd);
exp.gaborHeight = round(exp.stim.gaborHDeg*exp.ppd);                 % in pixels, the size of our objects
exp.gaborWidth = round(exp.stim.gaborWDeg*exp.ppd);                 % in pixels, the size of our objects

xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+exp.vertOffset;


%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers
topPhase = randi(360);
bottomPhase = topPhase;
exp.topWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth,length(exp.stim.orientation));
exp.bottomWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth,length(exp.stim.orientation));
exp.topWaveID = nan(length(flipTimes)*length(exp.longFormFlicker),length(exp.stim.orientation));
exp.bottomWaveID = nan(length(flipTimes)*length(exp.longFormFlicker),length(exp.stim.orientation));

for o =1:length(exp.stim.orientation)
    for f = 1:length(flipTimes)
        
        if f <= length(flipTimes)/2
            
            topPhase = mod(topPhase + exp.stim.dphasePerFlip,360);
            bottomPhase = mod(bottomPhase + exp.stim.dphasePerFlip,360);
            %ih in pixels %iw in pixels %spatial freq in cycles per dva
            %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
            %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
            %background color (unused if the grating is not an annulus)
            
            exp.topWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
                exp.stim.orientation(o),topPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
                exp.ppd);
            exp.bottomWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
                exp.stim.orientation(o),bottomPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
                exp.ppd);
            %    figure();
            %   imshow(squeeze(topWave(f,:,:)));
            %
        elseif f > length(flipTimes)/2
            topPhase = mod(topPhase - exp.stim.dphasePerFlip,360);
            bottomPhase = mod(bottomPhase - exp.stim.dphasePerFlip,360);
            exp.topWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
                exp.stim.orientation(o),topPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
                exp.ppd);
            exp.bottomWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
                exp.stim.orientation(o),bottomPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
                exp.ppd);
        end
        tmptopWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.topWave(f,:,:,o)));
        tmpbottomWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.bottomWave(f,:,:,o)));
        
    end
    
%% extend stimulus matrix to include the same total number of flips as the whole experiment
exp.topWaveID(:,o) = repmat(tmptopWaveID(:,o),length(exp.longFormFlicker),1);
exp.bottomWaveID(:,o) = repmat(tmpbottomWaveID(:,o),length(exp.longFormFlicker),1);

end

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xtop = rect(3)/2 - (8+3.5)*exp.ppd;%rect(3)/4; % = stimulus center located 8 degrees horizontal from the center
ytop = rect(4)/2 - exp.gaborHeight;
exp.topRect =  CenterRectOnPoint([0 0 exp.gaborWidth exp.gaborHeight],xtop,ytop);

xbottom = rect(3)/2 - (8+3.5)*exp.ppd; %rect(3)/4; = stimulus center located 8 degrees horizontal from the center
ybottom = rect(4)/2 + exp.gaborHeight;%rect(4)/4*3;
exp.bottomRect =  CenterRectOnPoint([0 0 exp.gaborWidth exp.gaborHeight],xbottom,ybottom);


%% %%%%%%%%%%%%%%%%%%%%%%
   % Letter task set-up %
   %%%%%%%%%%%%%%%%%%%%%%

%%%% find the target positioning for this run (index of repeated letters)
exp.numTargets = length(find(rand(1,exp.totalLetters)<exp.targetProb)); %draw targets from uniform distribution between 0 nd 1, only keep the density we need (0.1) of them
targetInd = zeros(1,exp.totalLetters+1); % this has to go to +1, but it's only for the repeat check
if exp.numTargets > 0
    while 1
        %check previous (-1,-2) and following indices (+1,+2) if they might constitute a
        %repeated letter
        maybeRep = exp.firstPossibleTarget+Randi(exp.totalLetters-exp.firstPossibleTarget-exp.lastPossibleTarget);    % a possible index for a target letter, picked randomly
        if targetInd(maybeRep-1) == 0 && targetInd(maybeRep+1) == 0 && targetInd(maybeRep-2) == 0 && targetInd(maybeRep+2) == 0 % make sure the previous  TWO and following TWO letters weren't already a repeat
            targetInd(maybeRep) = 1; %record this index if it is not a repeat
        end
        if sum(targetInd) == exp.numTargets %if the sum of target letters reach the number of targets, stop the check
            break
        end
    end
end
targetInd = targetInd(1:exp.totalLetters); % trim this back for sanity


%%%% fill the sequence out with actual letters
exp.letterSequence = [];
for k = 1:length(targetInd)
    if targetInd(k) == 1 % targets - we know these don't repeat
        exp.letterSequence{k} = exp.targetLetters{randi(length(exp.targetLetters))}; %pick a random target
    else % distractors - we need to make sure these don't repeat
        if k > 1 % the first letter can be whatever
            while 1
                maybeLetter = exp.distractors{randi(length(exp.distractors))};
                if strcmp(maybeLetter, exp.letterSequence{k-1}) == 0
                    exp.letterSequence{k} = maybeLetter;
                    break
                end
            end
        else %% (the first letter can be whatever)
            exp.letterSequence{k} = exp.distractors{randi(length(exp.distractors))};
        end
    end
end

% %%%% scale up the task based on flips
exp.letterSequence = Expand(exp.letterSequence,flipsPerTrial,1);%e.flipsPerSec,1); 

%% Eyetracking parameters
% eyetracking on (1) or off (0)
 ET = 0;
% EyelinkSetup(0);

%% %%%% initial window - wait for backtick
Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(w, 'Flip', 0);

DrawFormattedText(w,'Attend to the fixation circle: press 1 as soon as letter J appears on the screen,\n\n and press 2 as soon as letter K appears on the screen. \n\n Press Space to start'...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);

%%%%%%%%%%%%%%%%%% Response listening %%%%%%%%%%%%%%%%%%%%%%%%
%KbTriggerWait(53, deviceNumber);
KbTriggerWait(KbName('Space'), deviceNumber);
KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n=0;
count = 1;
%%%%%%% START task TASK/FLIPPING
% [e.totalTime, e.allFlips,n]
gray = repmat(min(min(squeeze(exp.topWave(1,:,:)),[],1)), [1,3]);
Screen('FillRect', w, gray);
    
if ET
    EyeData.mx{1}=[];
    EyeData.my{1}=[];
    EyeData.ma{1}=[];
    EyeData.FixDoneT{1} = [];
end
% EyeStart(1) = GetSecs(); %time we start caring about eyetracking
% currPosID = i;
% run FixCheck; %check eyetracker
while n+1 < length(exp.allFlips)

    [exp.longFormBlocks(n+1),exp.longFormFlicker(n+1)]
    thisCond = exp.longFormConds(n+1);

    KbQueueStart();
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 0
        [VBLT, exp.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        exp.flipTime(n+1) = exp.startRun;
    else
        [VBLT, exp.flipTime(n+1), FlipT, missed] = Screen(w, 'Flip', exp.startRun + exp.allFlips(n+1) - slack);
    end
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exp.longFormBlocks(n+1) == 1 && exp.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        if strfind(conditions(thisCond).name{:}, 'vert')
            ori =1;
        else
            ori = 2;
        end
        % draw & increment stims
        if strfind(conditions(thisCond).name{:}, 'DoubleIndir')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            % top stim
            Screen('DrawTexture', w, exp.topWaveID(n+1,ori),[],exp.topRect);
            % bottom stim
            Screen('DrawTexture', w, exp.bottomWaveID(n+1,ori), [], exp.bottomRect);
        elseif strfind(conditions(thisCond).name{:}, 'SingleTop')
            % top stim
            Screen('DrawTexture', w, exp.topWaveID(n+1,ori),[],exp.topRect);
        elseif strfind(conditions(thisCond).name{:}, 'SingleBottom')
            % bottom stim
            Screen('DrawTexture', w, exp.bottomWaveID(n+1,ori), [], exp.bottomRect);
        end
        
    end
    
    % select new character if starting new trial
    if mod(n, flipsPerTrial) == 0

%         % draw character
%         l = randperm(numel(letters));
%         fixChar = letters(l(1));
          fixChar = exp.letterSequence(n+1);

       if strcmp(fixChar,'J') == 1
            exp.targets = [exp.targets, 1];
            exp.targetTimes = [exp.targetTimes, GetSecs - exp.startRun];
        elseif strcmp(fixChar,'K') == 1
            exp.targets = [exp.targets, 2];
            exp.targetTimes = [exp.targetTimes, GetSecs - exp.startRun];
       end


    end
    
    %%%% draw fixation letter in fixation circle
    
    if mod(n, flipsPerTrial) <= exp.trialOnFlips
        Screen('FillOval', w,[255 255 255], [xc-round(exp.fixSize/2) yc-round(exp.fixSize/2) xc+round(exp.fixSize/2) yc+round(exp.fixSize/2)]);%white fixation solid circle
        %DrawFormattedText(w, fixChar, 'center', 8+vertOffset+rect(4)/2,0); %either text function works
        %Screen('DrawText', w, fixChar, -5+rect(3)/2, -10+vertOffset+rect(4)/2,[0 0 0]);
        [width,height] = RectSize(Screen('TextBounds',w,fixChar{:}));
        Screen('DrawText', w, fixChar{:}, xc-width/2 +0.3, yc-height/2-0.2,[0 0 0]);
     else
        Screen('FillOval', w,[255 255 255], [xc-round(exp.fixSize/2) yc-round(exp.fixSize/2) xc+round(exp.fixSize/2) yc+round(exp.fixSize/2)]);%white fixation solid circle
   
    end

%     n

    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    
    %%%% character identification
    if (pressed == 1) && ((firstPress(KbName('1!')) > 0) || (firstPress(KbName('2@')) > 0))
        if firstPress(KbName('1!')) > 0
            exp.responses = [exp.responses, 1];
            exp.responseTimes = [exp.responseTimes, firstPress(KbName('1!')) - exp.startRun];
        elseif firstPress(KbName('2@')) > 0
            exp.responses = [exp.responses, 2];
            exp.responseTimes = [exp.responseTimes, firstPress(KbName('2@')) - exp.startRun];
        end
    end

    %%%% refresh queue for next character
    KbQueueFlush();
    n = n+1;
   
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

exp.runTime = GetSecs - exp.startRun;

% analyse task accuracy
responseTimes = exp.responseTimes; % make copies as we are editing these
responses = exp.responses;
responseWindow = 1.5; % responses later than this count as a miss
for t = 1:size(exp.targetTimes,2)
    targetTime = exp.targetTimes(t);
    target = exp.targets(t);
    exp.hits(t) = 0; % default is a miss with RT of nan
    exp.RTs(t) = nan;
    for r = 1:size(responseTimes,2)
        if ismember(responses(r),  target) && responseTimes(r) > targetTime % if response is correct and happened after target
            rt = responseTimes(r)-targetTime;
            if rt < responseWindow % and if response happened within a second of target
                exp.hits(t) = 1; % mark a hit
                exp.RTs(t) = rt; % store the RT
                responseTimes(r) = []; % delete this response so it can't be reused
                responses(r) = [];
                break % stop looking for responses to this target
            end
        end
    end
end
exp.accuracy = (sum(exp.hits)/size(exp.targetTimes,2))*100;
exp.meanRT = nanmean(exp.RTs);

savedir = fullfile(exp.root,'data',subject,session,'fillingin_rsvp_v1');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, sprtinf('s%d_fillingin_rsvp_v1_sn%d_rn_%d_%d.mat',subject,num2str(exp.scanNum),num2str(exp.runNum),exp.date));
save(savename,'exp');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;


if ET
    Eyelink('StopRecording');   
    Eyelink('CloseFile');
    Eyelink('Shutdown');
    save ET_pilot.mat
end  

