function LD_phantom_mri_V5(subject, session, vertOffset, debug, figSizeDeg) 

%%MRI phantom task
%Loic 01312023
%In this version, we add multiple velocities
% subject = 'Dave';                                                                                                                                                                                                                                                     
% session = 1;                                                                                                                           
% debug = 1;
% vertOffset = 0;

ex.version = 'v5';
%global EyeData rect w xc yc %eye_used
%%%% resolution 
if debug == 1

    ex.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop 1920,1080
    ex.gammaCorrect = 0;       % make sure this = 1 when you're at the scanner!
else
                                                                                                                             
    ex.screenWidth = 17;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 48;             % in cm; %23 in eye tracking                                                                                                                          room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % scanner
    ex.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
    ex.scanNum = input('Scan number :');
    ex.runNum = input('Run number :');

end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('`~'))=1; % button box 2 % for backticks
Screen('Preference', 'SkipSyncTests', 0);

ex.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m
ex.whichCLUT = '/Users/tonglab/Desktop/Loic/CLUTs/7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';


%%% basic naming set-up
ex.subject = subject;
ex.session = session;


%%%% set-up rand
%  rand('twister', sum(100*clock));
%  ex.rand = rand;

rng(sum(100*clock));
ex.rand = rng;
%%%% files and things
ex.root = pwd;
ex.date = datestr(now,30);


%%%% 2D sine wave grating properties
ex.stim.spatialFreqDeg = 0.5/2; %0.286;                                          % cycles per degree of visual angle
% ex.stim.contrast = 0.15 ;                                                 % in %, maybe??
ex.stim.orientation = [90 180];                                                % in degrees
ex.stim.gaborHDeg = figSizeDeg;                                                  % in degrees of visual angle
ex.stim.gaborWDeg = figSizeDeg;
ex.stim.contrastMultiplicator = 0.075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.stim.contrastMults = linspace(0,0.075,20);                            %contrast multiplicators for the adjustable sine wave from 0 to 0.6 contrast
ex.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
ex.stim.cycPerSec = 1.20; % try multiple speeds
ex.stim.motionRate = ex.stim.cycPerSec.*360;                                          % 1.13 cycles per second = 360 deg of phase *1.13 per sec
ex.stim.cycles = 2; %number of cycles shifted per lap  (modified to half the number of cycles per lap)


%%%% sine wave grating timing (within block scale)
%ex.initialFixation = 6;        % in seconds
%ex.finalFixation = 2;          % in seconds
ex.stimDur = (ex.stim.cycles./ex.stim.cycPerSec)*2;        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
ex.stimsPerBlock = 3;%1;      % number of back-and-forth laps of the stimulus drift
ex.blockLength = ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.betweenBlocks = 16;          % in seconds
ex.flipsPerSec = 12; %60;  % 12;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 

%%%% conditions & layout (across blocks scale)

ex.conds = {'vertSingleLeft','horSingleLeft','vertSingleRight','horSingleRight','vertDoubleIndir','horDoubleIndir'};
ex.numConds = length(ex.conds);
% with line of code below we will have 6 conditions per block appear twice, randomized. 
ex.repsPerScan = 2;              % repetitions of each condition per run
ex.numBlocks = ex.numConds*ex.repsPerScan;
ex.condShuffle = Shuffle(repmat([1:ex.numConds],1,ex.repsPerScan)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order

%%%% Fixation
ex.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
ex.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 12; %26;

%% Color discrimination task setup
ex.targetProb = .4;              % proportion of trials where the target color will come up
ex.firstPossibleTarget = 5;     % no targets in the first X seconds
ex.lastPossibleTarget = 5;      % no targets in the last X seconds
ex.targetColors = {'red'};
ex.distractColors = {'green'};
ex.trialDur = 2;
ex.targetOn = 0.5;
% exp.totalColors = (exp.totalTime/exp.trialDur);
%ex.responseBack = 1;    % the response is correct if the preceding N colors were the target
ex.flipsPerTrial = ex.trialDur/ex.flipWin;
ex.trialOnFlips = floor(ex.targetOn/ex.flipWin);
%%%% color identification task
ex.targets = [];
ex.distractors = [];
ex.targetTimes = [];
ex.responses = [];
ex.responseTimes=[];
ex.accuracy = 0;
ex.meanRT = 0;

taskText = 'colors';



%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

ex.onSecs = [ones(1,ex.blockLength) zeros(1,ex.betweenBlocks)];
ex.longFormBlocks = Expand(ex.onSecs,ex.flipsPerSec,1); %1 when block, 0 when between block

% flicker = [ones(1,ex.blockLength) zeros(1,ex.betweenBlocks)]; 
% ex.longFormFlicker = repmat(ones(1,1),1,length(ex.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(ex.longFormBlocks)


%%
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%
HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
if debug
     [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
     %  [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
else   
    [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
end
Screen(w, 'TextSize', ex.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
if ex.gammaCorrect > 0
    load(ex.whichCLUT);
    Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
end

%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate =  1/frameInt;%Screen('NominalFrameRate',w);

flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.blockLength]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
stimFlipTimes = flipTimes(1:length(flipTimes)-1);
ex.stimFlipTimes = stimFlipTimes;

%%% Stimulus phases
% (in a linear/triangle scenario) ex.stim.dphasePerFlip = ex.stim.motionRate*frameInt * frameRate/ex.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip since we have 5 times less flips than frame refresh
ex.stim.tempPhase =pi; % pi/2; %this will make the oscillation start on the fast segment of the oscillation, with red dot on the center of the screen

ex.stim.oscillation = cos(2*pi*(1/ex.stimDur)*stimFlipTimes+ex.stim.tempPhase);
ex.stim.spatialPhase = 90; %this makes the grating phase so that the center of the dark stripe is on the center of the screen
ex.stim.phases = ex.stim.oscillation.*180*ex.stim.cycles+ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory

%reasonning behind the calculation of ex.stim.phases:
%GOAL: render the back and forth of grating drifts oscillatory in time instead
%of linear to smooth the signal phase shifts at the time of drift direction
%reversal
%1) The length (x) of grating shift after each flip is the variable that will have to be oscillatory
%over time

%2)oscillatory signal formula is x(t) = A*cos(w*t+phase) (x is the oscillatory signal, w is the angular
%frequency in rad/s, t is the time in s, phase is a phase shift in rad at time t (here 0),
%A is the amplitude of the signal

%3)Here the angular frequency is calculated based on the idea of one
%oscillation per stimulus. Considering a stimulus as 1 cycle forward and 1 cycle backwards(definition might change based on how many cycles are desired per stimulus), one oscillation will represent 2 cycles =ex.stim.cycPerSec*ex.stimDur)
% thus the temporal frequency of the oscillation should be 1/ex.stimDur
% (this is w)

%%4) We multiply w by t. Here the resolution of the signal is given by the number of flips per second, so the duration is given by flipTimes.
%Finally, we multiply by 2*pi for the units to be in radian.

%5)Then we subtract the phase at each flip at time t which here is 0.

%6)The amplitude of cos() varies between -1 and +1. Here we want our
%displacement (grating shift) in degrees with a displacement of 360? (of the spatial grating) every half temporal oscillation so we multiply by 180 for the range to be [-180;+180].

%7) we multiply by the number of cycles desired to drift over one lap
%(ex.stim.cycles).

%%%% scale the stim params for the screen
ex.ppd = pi* rect(3) / (atan(ex.screenWidth/ex.viewingDist/2)) / 360;
ex.fixSize = round(ex.fixSizeDeg*ex.ppd);
ex.gaborHeight = round(ex.stim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.gaborWidth = round(ex.stim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects

%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers
%rect
for o =1:length(ex.stim.orientation)
    stimFlipTimes = ex.stimFlipTimes;
    
%     LWaveIDO = nan(length(flipTimes),length(ex.stimDur));
%     RWaveIDO = nan(length(flipTimes),length(ex.stimDur));
%     
    ex.LWave = nan(length(stimFlipTimes),ex.gaborHeight,ex.gaborWidth,length(ex.stim.orientation));
    ex.RWave = nan(length(stimFlipTimes),ex.gaborHeight,ex.gaborWidth,length(ex.stim.orientation));
    
    phases = ex.stim.phases ;

    for f = 1:length(stimFlipTimes)
        
        LPhase = phases(f);
        RPhase = phases(f);
        %ih in pixels %iw in pixels %spatial freq in cycles per dva
        %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
        %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
        %background color (unused if the grating is not an annulus)
        
        ex.LWave(f,:,:,o) = makeSineGrating(ex.gaborHeight,ex.gaborWidth,ex.stim.spatialFreqDeg,...
            ex.stim.orientation(o),LPhase,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
            ex.ppd);
        ex.RWave(f,:,:,o) = makeSineGrating(ex.gaborHeight,ex.gaborWidth,ex.stim.spatialFreqDeg,...
            ex.stim.orientation(o),RPhase,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
            ex.ppd);
        %                 figure();
        %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
        %
        tmpLWaveID(f,o) = Screen('MakeTexture', w, squeeze(ex.LWave(f,:,:,o)));
        tmpRWaveID(f,o) = Screen('MakeTexture', w, squeeze(ex.RWave(f,:,:,o)));
        
    end
    
    ex.LWaveID(1:length(tmpLWaveID(:,o)),o) = tmpLWaveID(:,o);
    ex.RWaveID(1:length(tmpRWaveID(:,o)),o) = tmpRWaveID(:,o);

end
 
%% add between block flips
betweenFlipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.betweenBlocks]+flipTimes(end);
    ex.LWaveID = [ex.LWaveID; zeros(length(betweenFlipTimes)-1, length(ex.stim.orientation))];
    ex.RWaveID = [ex.RWaveID; zeros(length(betweenFlipTimes)-1, length(ex.stim.orientation))];
ex.betweenFlipTimes = betweenFlipTimes(1:length(betweenFlipTimes)-1);

ex.allFlipTimes = [ex.stimFlipTimes,ex.betweenFlipTimes];
%% Sine wave gratings locations (in the task loop since it changes)
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xL = rect(3)/2- ex.gaborWidth;  % = stimulus center located on the horizontal center of the screen
yL = rect(4)/2; % stimulus located 4 degrees above screen center

xR = rect(3)/2+ ex.gaborWidth; % = stimulus center located on the horizontal center of the screen
yR = rect(4)/2; % stimulus located 4 degrees below screen center

%% %%%%%%%%%%%%%%%%%%%%%%
   % Letter task set-up %
   %%%%%%%%%%%%%%%%%%%%%%

%%%% find the target positioning for this run (index of repeated letters)
% ex.numTargets = length(find(rand(1,ex.totalColors)<ex.targetProb)); %draw targets from uniform distribution between 0 nd 1, only keep the density we need (0.1) of them
% targetInd = zeros(1,ex.totalColors+1); % this has to go to +1, but it's only for the repeat check
% if ex.numTargets > 0
%     while 1
%         %check previous (-1,-2) and following indices (+1,+2) if they might constitute a
%         %repeated color
%         maybeRep = ex.firstPossibleTarget+Randi(ex.totalColors-ex.firstPossibleTarget-ex.lastPossibleTarget);    % a possible index for a target letter, picked randomly
%         if targetInd(maybeRep-1) == 0 && targetInd(maybeRep+1) == 0 && targetInd(maybeRep-2) == 0 && targetInd(maybeRep+2) == 0 % make sure the previous  TWO and following TWO colors weren't already a repeat
%             targetInd(maybeRep) = 1; %record this index if it is not a repeat
%         end
%         if sum(targetInd) == ex.numTargets %if the sum of target letters reach the number of targets, stop the check
%             break
%         end
%     end
% end
% targetInd = targetInd(1:ex.totalColors); % trim this back for sanity


%%%% fill the sequence out with actual Colors
% ex.colorSequence = [];
% for k = 1:length(targetInd)
%     if targetInd(k) == 1 % targets - we know these don't repeat
%         ex.colorSequence{k} = ex.targetColors{randi(length(ex.targetColors))}; %pick a random target
%     else % distractors 
%         ex.colorSequence{k} = ex.distractors{randi(length(ex.distractors))};
%     end
% end

% %%%% scale up the task based on flips
% ex.colorSequence = Expand(ex.colorSequence,ex.flipsPerTrial,1);%e.flipsPerSec,1); 

%% %%%% initial window - wait for backtick
% Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
% Screen(w, 'Flip', 0);
% KbTriggerWait(53, deviceNumber);
DrawFormattedText(w,'Attend to the fixation cue in the center of the screen: press 1 as soon as color red appears on the screen. '... %\n\n (use this to change lines)
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
WaitSecs(10);

% %%%% response listening 

KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% START task TASK/FLIPPING

gray = repmat(min(min(squeeze(ex.RWave(1,:,:)),[],1)), [1,3]);

Screen('FillRect', w, gray);
n = 1; % block flips + between blocks flips
cnt = 1; %block number
backtickBlock = 0; %record backticks emitted by scanner machine
fixCol = {'NaN'}; %initialize color
% tend = [];
probs = [1]; %first probability is 1 to allow for potential next target to occur
probsCnt = 1; %start from 1 since we already have probs = [1]
while(1)
    KbQueueStart();
    %%%%%%% specify condition to draw at the start of the next block
    if n == 1 && cnt == 1 %for first block
        start = GetSecs();
        Screen('FillRect', w, gray);
        Screen(w, 'Flip', 0);
    end
    thisCond = ex.condShuffle(cnt);
    
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ex.longFormBlocks(n) == 1
        ex.LRect =  CenterRectOnPoint([0 0 ex.gaborWidth ex.gaborHeight],xL,yL);
        ex.RRect =  CenterRectOnPoint([0 0 ex.gaborWidth ex.gaborHeight],xR,yR);
        if strfind(ex.conds{thisCond}, 'vert')
            ori =1;
        else
            ori = 2;
        end
        % draw & increment stims
        if strfind(ex.conds{thisCond}, 'DoubleIndir')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            % top stim
            Screen('DrawTexture', w, ex.LWaveID(n,ori),[],ex.LRect);
            % bottom stim
            Screen('DrawTexture', w, ex.RWaveID(n,ori), [], ex.RRect);
        elseif strfind(ex.conds{thisCond}, 'SingleLeft')
            % top stim
            Screen('DrawTexture', w, ex.LWaveID(n,ori),[],ex.LRect);
        elseif strfind(ex.conds{thisCond}, 'SingleRight')
            % bottom stim
            Screen('DrawTexture', w, ex.RWaveID(n,ori), [], ex.RRect);
        end
    end
    %%%%%%%%% Color discrimination task %%%%%
    % select new color if starting new trial
    if n < ex.flipsPerTrial
        fixCol = ex.distractColors(randi(length(ex.distractColors)));
    elseif mod(n, ex.flipsPerTrial) == 0
        prob = rand(1);
        probs = [probs, prob];
        probsCnt = probsCnt + 1;
        
        if prob < ex.targetProb %&& probs(probsCnt-1) >= ex.targetProb % make sure previous trial did not have a target
            fixCol = ex.targetColors(randi(length(ex.targetColors)));
            ex.targets = [ex.targets, 1];
            ex.targetTimes = [ex.targetTimes, GetSecs - start];
        elseif prob >= ex.targetProb %|| probs(probsCnt-1) < ex.targetProb %if previous trial already had a target, show green
            fixCol = ex.distractColors(randi(length(ex.distractColors)));
            ex.distractors = [ex.distractors, 1];
            
        end
    end
    
    %%%% draw fixation color in center of the screen
    
    if mod(n, ex.flipsPerTrial) <= ex.trialOnFlips && strcmp(fixCol{:}, 'red')
        Screen('FillOval', w,[255 0 0], [xc-round(ex.fixSize/4) yc-round(ex.fixSize/4) xc+round(ex.fixSize/4) yc+round(ex.fixSize/4)]);%red fixation solid circle
        
    elseif (mod(n, ex.flipsPerTrial) <= ex.flipsPerTrial && strcmp(fixCol{:}, 'green')) || (mod(n, ex.flipsPerTrial) > ex.trialOnFlips && strcmp(fixCol{:}, 'red')) % also accounts for after the target flashed
        Screen('FillOval', w,[0 255 0], [xc-round(ex.fixSize/4) yc-round(ex.fixSize/4) xc+round(ex.fixSize/4) yc+round(ex.fixSize/4)]);%green fixation solid circle
        
    end
    
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if n == 1
        [~,keyCode,~] = KbWait(deviceNumber,2);
    end
    
    if (n == 1 && cnt ==1) && nnz(keyCode(53)) %backTick > 0%%%%%%%
        
        [VBLT, ex.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        ex.flipTime(n,cnt) = ex.startRun;
        ex.blockOnsetTime(n,cnt) = ex.startRun;
        start = ex.startRun;
        backtickBlock = 1;
        n = n+1;
       
    elseif (n == 1 && cnt ~= 1) && nnz(keyCode(53)) %backTick > 0 %% %%%%   %use second condition if we wait for backticks from scanner
        
        [VBLT, ex.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        ex.flipTime(n,cnt) = ex.startRun;
        ex.blockOnsetTime(n,cnt) = ex.startRun;
        start = ex.startRun;
        backtickBlock = 1;
        n = n+1;

    elseif (n ~= 1 && cnt ~= 1) && backtickBlock == 1
        [VBLT, ex.flipTime(n,cnt), FlipT, missed] = Screen(w, 'Flip', start+ ex.allFlipTimes(n) - slack);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        
    elseif n ~= 1 && cnt == 1
        [VBLT, ex.flipTime(n,cnt), FlipT, missed] = Screen(w, 'Flip', start+ ex.allFlipTimes(n) - slack);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        
    end
    
    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    %     %%%% color identification
    if (pressed == 1) && (firstPress(KbName('1!')) > 0) %|| (firstPress(KbName('2@')) > 0))
        ex.responses = [ex.responses, 1];
        ex.responseTimes = [ex.responseTimes, firstPress(KbName('1!')) - start];
    end
    
    %%%% refresh queue for next character
    KbQueueFlush();

    if n ~= 1
        n = n+1;
    end
    
    if (n == (length(ex.allFlipTimes))) % for any other block , reset frame index when previous trial ends
        n = 1;
        cnt = cnt+1;
        backtickBlock = 0;
        start = GetSecs; %reset time
    end
    if cnt == ex.numBlocks+1
        break;
    end
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%


ex.runTime = GetSecs - ex.startRun;

% analyse task accuracy
responseTimes = ex.responseTimes; % make copies as we are editing these
responses = ex.responses;
responseWindow = 1.5; % responses later than this count as a miss
for t = 1:size(ex.targetTimes,2)
    targetTime = ex.targetTimes(t);
    target = ex.targets(t);
    ex.hits(t) = 0; % default is a miss with RT of nan
    ex.RTs(t) = nan;
    for r = 1:size(responseTimes,2)
        if ismember(responses(r),  target) && responseTimes(r) > targetTime % if response is correct and happened after target
            rt = responseTimes(r)-targetTime;
            if rt < responseWindow % and if response happened within a second of target
                ex.hits(t) = 1; % mark a hit
                ex.RTs(t) = rt; % store the RT
                responseTimes(r) = []; % delete this response so it can't be reused
                responses(r) = [];
                break % stop looking for responses to this target
            end
        end
    end
end
ex.accuracy = (sum(ex.hits)/size(ex.targets,2))*100;
ex.meanRT = nanmean(ex.RTs);
disp(sprintf('Accuracy: %d', ex.accuracy));
disp(sprintf('Mean reaction time: %d', ex.meanRT));

savedir = fullfile(ex.root,'data',subject,session,'phantom_v1');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir,sprintf('/s%s_phantom_v1_sn%d_rn%d_date%s.mat',subject,ex.scanNum,ex.runNum,num2str(ex.date)));
save(savename,'ex');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;

