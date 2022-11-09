function LD_phantom_mri_V2(subject, session, vertOffset, debug) 

% subject = 1;
% session = 1;
% debug = 0;
% vertOffset = 0;


global EyeData currPosID rect w xc yc
%%%% resolution
if debug == 1
    exp.screenWidth = 27.5;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    exp.viewingDist = 57;             % in cm; %3Tb/office=43, miniHelm=57;
    exp.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % scanner
    exp.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!

else
    exp.screenWidth = 27.5;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    exp.viewingDist = 57;             % in cm; 3Tb/office=43, miniHelm=57;
    exp.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop
    exp.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!


end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
%responseKeys(KbName('2@'))=1; % button box 2

Screen('Preference', 'SkipSyncTests', 0);

exp.scanNum = input('Scan number :');
exp.runNum = input('Run number :');
exp.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m
exp.whichCLUT = '/Users/tonglab/Desktop/Loic/CLUTs/7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
exp.subject = subject;
exp.session = session;


%%%% set-up rand
%  rand('twister', sum(100*clock));
%  exp.rand = rand;

rng(sum(100*clock));
exp.rand = rng;
%%%% files and things
exp.root = pwd;
exp.date = datestr(now,30);

%%%% 2D sine wave grating properties
exp.stim.spatialFreqDeg = 0.286;  %1                                         % cycles per degree of visual angle
exp.stim.contrast = 0.3 ;                                                 % in %, maybe??
exp.stim.orientation = [90 180];                                                % in degrees                                              % in degrees of visual angle
exp.stim.gaborHDeg = 7;                                                  % in degrees of visual angle
exp.stim.gaborWDeg = 7; 
exp.stim.contrastMultiplicator = .075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
exp.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
exp.stim.cycPerSec = 1.13; % (might need to modified to half the speed)
exp.stim.motionRate = exp.stim.cycPerSec*360;                                          % 1.13 cycles per second = 360 deg of phase *1.13 per sec
exp.stim.cycles = 2; %number of cycles shifted per lap  (modified to half the number of cycles per lap)


%%%% sine wave grating timing (within block scale)
exp.initialFixation = 6;        % in seconds
exp.finalFixation = 2;          % in seconds
exp.stimDur = (exp.stim.cycles/exp.stim.cycPerSec)*2; %1.6*2;        % in seconds. 1.6 sec refers to sine wave grating 1.6 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
exp.stimsPerTrial = 3;          % number of back-and-forth laps of the stimulus drift
exp.trialLength =  ceil(exp.stimDur*exp.stimsPerTrial);             % in seconds
exp.betweenTrial= 16;          % in seconds
exp.flipsPerSec = 12;           % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
exp.flipWin = 1/exp.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 

%%%% conditions & layout (across blocks scale)
exp.conds = {'vertSingleLeft','horSingleLeft','vertSingleRight','horSingleRight','vertDoubleIndir','horDoubleIndir'}; %'double-oppdir'
exp.numConds = length(exp.conds);
exp.repsPerScan = 2;              % repetitions of each condition per run
exp.numTrials = exp.numConds*exp.repsPerScan;
exp.condShuffle = Shuffle(repmat([1:exp.numConds],1,exp.repsPerScan)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order

exp.totalTime = exp.initialFixation + ((exp.numTrials-1) * (exp.trialLength + exp.betweenTrial)) + exp.trialLength + exp.finalFixation;
exp.allFlips = (0:exp.flipWin:exp.totalTime);

%%%% Fixation
exp.fixSizeDeg =  .6;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
exp.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
exp.fontSize = 12; %26;

%% Color discrimination task setup
exp.targetProb = .2;              % proportion of trials where the target color will come up
exp.firstPossibleTarget = 5;     % no targets in the first X flips
exp.lastPossibleTarget = 5;      % no targets in the last X flips
exp.targetColors = {'red'};
exp.distractors = {'green'};
exp.trialDur = 2;
exp.totalColors = (exp.totalTime/exp.trialDur);
exp.responseBack = 1;    % the response is correct if the preceding N colors were the target
exp.flipsPerTrial = exp.trialDur/exp.flipWin;
exp.trialOnFlips = floor(exp.trialDur/exp.flipWin);
%%%% color identification task
exp.targets = [];
exp.targetTimes = [];
exp.responses = [];
exp.responseTimes=[];
exp.accuracy = 0;
exp.meanRT = 0;

taskText = 'colors';



%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

exp.onSecs = [zeros(1,exp.initialFixation)...
    repmat([ones(1,exp.trialLength) zeros(1,exp.betweenTrial)],1,exp.numTrials-1)... %2*ones(1,e.trialLength) zeros(1,e.betweenBlocks)
    ones(1,exp.trialLength) zeros(1,exp.finalFixation)];
exp.longFormBlocks = Expand(exp.onSecs,exp.flipsPerSec,1); %1 when block, 0 when between block
exp.longFormFlicker = repmat(ones(1,1),1,length(exp.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(exp.longFormBlocks)


%set up the timing model for stimulus pairing conditions ( pair square, pair rectangle ) and stimulus orientation
%longform condition timing, which aligns with the flicker timing
exp.longFormConds = zeros(1,exp.initialFixation);
for i = 1:exp.numTrials-1
    exp.longFormConds = [exp.longFormConds, repmat(exp.condShuffle(i),1,exp.trialLength)]; % blocks
    exp.longFormConds = [exp.longFormConds, zeros(1,exp.betweenTrial)]; % inter-block blanks
    
end
exp.longFormConds = [exp.longFormConds, repmat(exp.condShuffle(end),1,exp.trialLength), zeros(1,exp.finalFixation)]; % the last block
exp.longFormConds = Expand(exp.longFormConds, exp.flipsPerSec,1);
length(exp.longFormConds)

% %% create the timing model of stimulus conditions for this particular run
% clear i
% for i =1:exp.numConds
%     conditions(i).name = exp.conds(i);
%     conditions(i).startTimes = [];
% end
%%
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

%HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
if debug == 1
    [w, rect]=Screen('OpenWindow',screen,exp.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
else
    [w, rect]=Screen('OpenWindow',screen,exp.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
end
Screen(w, 'TextSize', exp.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
if exp.gammaCorrect > 0
    load(exp.whichCLUT);
    Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
end

%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2; %used to improve flip timing accuracy
frameRate = 1/frameInt;%Screen('NominalFrameRate',w);
flipTimes = [0:frameInt*frameRate/exp.flipsPerSec:exp.stimDur]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
exp.flipTimes = flipTimes;
%(in a linear/triangle scenario) exp.stim.dphasePerFlip = exp.stim.motionRate*frameInt * frameRate/exp.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip
%%%% timing optimization
exp.stim.oscillation = cos(2*pi*(1/exp.stimDur)*flipTimes);
exp.stim.spatialPhase = 90;
exp.stim.phases = exp.stim.oscillation.*180*exp.stim.cycles+exp.stim.spatialPhase; %./exp.stimDur-2*pi*flipTimes./exp.stimDur make it oscillatory
%reasonning behind the calculation of exp.stim.phases:
%GOAL: render the back and forth of grating drifts oscillatory in time instead
%of linear to smooth the signal phase shifts at the time of drift direction
%reversal
%1) The length (x) of grating shift after each flip is the variable that will have to be oscillatory
%over time

%2)oscillatory signal formula is x(t) = A*cos(w*t+phase) (x is the oscillatory signal, w is the angular
%frequency in rad/s, t is the time in s, phase is a phase shift in rad at time t (here 0),
%A is the amplitude of the signal

%3)Here the angular frequency is calculated based on the idea of one
%oscillation per stimulus. Considering a stimulus as 1 cycle forward and 1 cycle backwards(definition might change based on how many cycles are desired per stimulus), one oscillation will represent 2 cycles =exp.stim.cycPerSec*exp.stimDur)
% thus the temporal frequency of the oscillation should be 1/exp.stimDur
% (this is w)

%%4) We multiply w by t. Here the resolution of the signal is given by the number of flips per second, so the duration is given by flipTimes.
%Finally, we multiply by 2*pi for the units to be in radian.

%5)Then we subtract the phase at each flip at time t which here is 0.

%6)The amplitude of cos() varies between -1 and +1. Here we want our
%displacement (grating shift) in degrees with a displacement of 360? (of the spatial grating) every half temporal oscillation so we multiply by 180 for the range to be [-180;+180].

%7) we multiply by the number of cycles desired to drift over one lap
%(exp.stim.cycles).


%%%% scale the stim params for the screen
exp.ppd = pi* rect(3) / (atan(exp.screenWidth/exp.viewingDist/2)) / 360;
exp.fixSize = round(exp.fixSizeDeg*exp.ppd);
exp.gaborHeight = round(exp.stim.gaborHDeg*exp.ppd);                 % in pixels, the size of our objects
exp.gaborWidth = round(exp.stim.gaborWDeg*exp.ppd);                 % in pixels, the size of our objects


%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers

exp.LWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth,length(exp.stim.orientation));
exp.RWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth,length(exp.stim.orientation));
exp.LWaveID = nan(length(exp.longFormFlicker),length(exp.stim.orientation));
exp.RWaveID = nan(length(exp.longFormFlicker),length(exp.stim.orientation));

for o =1:length(exp.stim.orientation)
    for f = 1:length(flipTimes)
            
        LPhase = exp.stim.phases(f);
        RPhase = exp.stim.phases(f);            %ih in pixels %iw in pixels %spatial freq in cycles per dva
            %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
            %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
            %background color (unused if the grating is not an annulus)
            
            exp.LWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
                exp.stim.orientation(o),LPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
                exp.ppd);
            exp.RWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
                exp.stim.orientation(o),RPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
                exp.ppd);      
            %    figure();
            %   imshow(squeeze(topWave(f,:,:)));
        tmpLWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.LWave(f,:,:,o)));
        tmpRWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.RWave(f,:,:,o)));
    end
    
%% extend stimulus matrix to include the same total number of flips as the whole experiment

exp.LWaveID(:,o) = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
    repmat([repmat(tmpLWaveID(:,o)',1,floor(exp.trialLength*exp.flipsPerSec/length(tmpLWaveID(:,o)))) tmpLWaveID(1:mod(exp.trialLength*exp.flipsPerSec,length(tmpLWaveID(:,o))),o)' zeros(1,exp.betweenTrial*exp.flipsPerSec)],1,exp.numTrials-1)... %2*ones(1,e.trialLength) zeros(1,e.betweenBlocks)
    repmat(tmpLWaveID(:,o)',1,floor(exp.trialLength*exp.flipsPerSec/length(tmpLWaveID(:,o)))) tmpLWaveID(1:mod(exp.trialLength*exp.flipsPerSec,length(tmpLWaveID(:,o))),o)' zeros(1,exp.finalFixation*exp.flipsPerSec)];

exp.RWaveID(:,o) = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
    repmat([repmat(tmpRWaveID(:,o)',1,floor(exp.trialLength*exp.flipsPerSec/length(tmpRWaveID(:,o)))) tmpRWaveID(1:mod(exp.trialLength*exp.flipsPerSec,length(tmpRWaveID(:,o))),o)' zeros(1,exp.betweenTrial*exp.flipsPerSec)],1,exp.numTrials-1)... %2*ones(1,e.trialLength) zeros(1,e.betweenBlocks)
    repmat(tmpRWaveID(:,o)',1,floor(exp.trialLength*exp.flipsPerSec/length(tmpRWaveID(:,o)))) tmpRWaveID(1:mod(exp.trialLength*exp.flipsPerSec,length(tmpRWaveID(:,o))),o)' zeros(1,exp.finalFixation*exp.flipsPerSec)];

end

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xL = rect(3)/2 - exp.gaborWidth; %rect(3)/4; % = stimulus center located 8 degrees horizontal from the center
yL = rect(4)/2; %- exp.gaborHeight;
exp.LRect =  CenterRectOnPoint([0 0 exp.gaborWidth exp.gaborHeight],xL,yL);

xR = rect(3)/2 + exp.gaborWidth; %rect(3)/4; = stimulus center located 8 degrees horizontal from the center
yR = rect(4)/2; %+ exp.gaborHeight;%rect(4)/4*3;
exp.RRect =  CenterRectOnPoint([0 0 exp.gaborWidth exp.gaborHeight],xR,yR);


%% %%%%%%%%%%%%%%%%%%%%%%
   % Letter task set-up %
   %%%%%%%%%%%%%%%%%%%%%%

%%%% find the target positioning for this run (index of repeated letters)
exp.numTargets = length(find(rand(1,exp.totalColors)<exp.targetProb)); %draw targets from uniform distribution between 0 nd 1, only keep the density we need (0.1) of them
targetInd = zeros(1,exp.totalColors+1); % this has to go to +1, but it's only for the repeat check
if exp.numTargets > 0
    while 1
        %check previous (-1,-2) and following indices (+1,+2) if they might constitute a
        %repeated color
        maybeRep = exp.firstPossibleTarget+Randi(exp.totalColors-exp.firstPossibleTarget-exp.lastPossibleTarget);    % a possible index for a target letter, picked randomly
        if targetInd(maybeRep-1) == 0 && targetInd(maybeRep+1) == 0 && targetInd(maybeRep-2) == 0 && targetInd(maybeRep+2) == 0 % make sure the previous  TWO and following TWO colors weren't already a repeat
            targetInd(maybeRep) = 1; %record this index if it is not a repeat
        end
        if sum(targetInd) == exp.numTargets %if the sum of target letters reach the number of targets, stop the check
            break
        end
    end
end
targetInd = targetInd(1:exp.totalColors); % trim this back for sanity


%%%% fill the sequence out with actual Colors
exp.colorSequence = [];
for k = 1:length(targetInd)
    if targetInd(k) == 1 % targets - we know these don't repeat
        exp.colorSequence{k} = exp.targetColors{randi(length(exp.targetColors))}; %pick a random target
    else % distractors 
        exp.colorSequence{k} = exp.distractors{randi(length(exp.distractors))};
    end
end

% %%%% scale up the task based on flips
exp.colorSequence = Expand(exp.colorSequence,exp.flipsPerTrial,1);%e.flipsPerSec,1); 

%% %%%% initial window - wait for backtick

DrawFormattedText(w,'Attend to the fixation cue in the center of the screen: press 1 as soon as color red appears on the screen. '... %\n\n (use this to change lines)
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
WaitSecs(10);
Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(w, 'Flip', 0);
KbTriggerWait(53, deviceNumber);

%%%%%%%%%%%%%%%%%% Response listening %%%%%%%%%%%%%%%%%%%%%%%%
%KbTriggerWait(53, deviceNumber);
%KbTriggerWait(KbName('Space'), deviceNumber);
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
gray = repmat(min(min(squeeze(exp.RWave(1,:,:)),[],1)), [1,3]);
Screen('FillRect', w, gray);
    
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
        if strfind(exp.conds{thisCond}, 'vert')
            ori =1;
        else
            ori = 2;
        end
        % draw & increment stims
        if strfind(exp.conds{thisCond}, 'DoubleIndir')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            % top stim
            Screen('DrawTexture', w, exp.LWaveID(n+1,ori),[],exp.LRect);
            % bottom stim
            Screen('DrawTexture', w, exp.RWaveID(n+1,ori), [], exp.RRect);
        elseif strfind(exp.conds{thisCond}, 'SingleLeft')
            % top stim
            Screen('DrawTexture', w, exp.LWaveID(n+1,ori),[],exp.LRect);
        elseif strfind(exp.conds{thisCond}, 'SingleRight')
            % bottom stim
            Screen('DrawTexture', w, exp.RWaveID(n+1,ori), [], exp.RRect);
        end
        
    end
    
    % select new color if starting new trial
    if mod(n, exp.flipsPerTrial) == 0

          fixCol = exp.colorSequence(n+1);

       if strcmp(fixCol,'red') == 1
            exp.targets = [exp.targets, 1];
            exp.targetTimes = [exp.targetTimes, GetSecs - exp.startRun];

       end


    end
    
    %%%% draw fixation color in center of the screen
    
    if mod(n, exp.flipsPerTrial) <= exp.trialOnFlips && strcmp(fixCol{:}, 'red')
        Screen('FillOval', w,[255 0 0], [xc-round(exp.fixSize/4) yc-round(exp.fixSize/4) xc+round(exp.fixSize/4) yc+round(exp.fixSize/4)]);%red fixation solid circle
        
    elseif mod(n, exp.flipsPerTrial) <= exp.trialOnFlips && strcmp(fixCol{:}, 'green')
        Screen('FillOval', w,[0 255 0], [xc-round(exp.fixSize/4) yc-round(exp.fixSize/4) xc+round(exp.fixSize/4) yc+round(exp.fixSize/4)]);%green fixation solid circle
        
    end

%     n

    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    
    %%%% color identification
    if (pressed == 1) && ((firstPress(KbName('1!')) > 0) || (firstPress(KbName('2@')) > 0))
        if firstPress(KbName('1!')) > 0
            exp.responses = [exp.responses, 1];
            exp.responseTimes = [exp.responseTimes, firstPress(KbName('1!')) - exp.startRun];
%         elseif firstPress(KbName('2@')) > 0
%             exp.responses = [exp.responses, 2];
%             exp.responseTimes = [exp.responseTimes, firstPress(KbName('2@')) - exp.startRun];
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
disp(sprintf('Accuracy: %d', exp.accuracy));
disp(sprintf('Mean reaction time: %d', exp.meanRT));



savedir = fullfile(exp.root,'data',subject,session,'phantom_v1');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir,sprintf('/s%s_phantom_v1_sn%d_rn%d_date%s.mat',subject,exp.scanNum,exp.runNum,num2str(exp.date)));
save(savename,'exp');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;



