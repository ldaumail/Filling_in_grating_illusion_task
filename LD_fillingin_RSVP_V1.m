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

% experiment.scanNum = input('Scan number :');
% experiment.runNum = input('Run number :');
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
experiment.initialFixation = 6;        % in seconds
experiment.finalFixation = 0;          % in seconds
experiment.stimDur = 2.6;        % in seconds,refers to sine wave grating
experiment.flipsPerSec = 12;           % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
experiment.flipWin = 1/experiment.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
experiment.numBlocks = 12;  % 6 for single and 6 for pair...
experiment.trialFreq = 1;               % duration of fixation trials (seconds) (time it takes to switch from red to green then back to red)
experiment.trialDur = .4;               % duration in seconds of the letter presentation of each fixation trial (.4s letter ON, .6 letter OFF)
experiment.postTargetWindow = 1;           % duration in seconds of target-free interval after a target
flipsPerTrial = experiment.trialFreq/experiment.flipWin;
trialOnFlips = experiment.trialDur/experiment.flipWin;

%%%% 2D sine wave grating properties
experiment.stim.spatialFreqDeg = 1;                                           % cycles per degree of visual angle
experiment.stim.contrast =  .3;                                                 % in %, maybe??
experiment.stim.orientation = 90;                                                % in degrees
experiment.stim.degFromFix = .6;                                              % in degrees of visual angle
experiment.stim.gaborHDeg = 4;                                                  % in degrees of visual angle
experiment.stim.gaborWDeg = 4; 
experiment.stim.contrastMultiplicator = .2;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
experiment.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
experiment.stim.motionRate = 1.3*360 ;                                          % 1.3 cycles per second = 360 deg of phase *1.3 per sec

%%%% conditions & layout
experiment.fixSizeDeg =  .4;            % in degrees, the size of the biggest white dot in the fixation
experiment.littleFixDeg = .25;    % proportion of the fixSizeDeg occupied by the smaller black dot
experiment.outerFixPixels = 2;          % in pixels, the black ring around fixation
%experiment.TRlength = 2;                % in seconds
experiment.repsPerRun = 2;              % repetitions of each object type x experimentation
experiment.totalTime = experiment.initialFixation + (experiment.numBlocks * (experiment.blockLength + experiment.betweenBlocks)) + experiment.finalFixation;
experiment.allFlips = (0:experiment.flipWin:experiment.totalTime);

%%%% screen
experiment.backgroundColor = [127 127 127];  % color
experiment.fontSize = 10; %26;

%%%%%%%%%%%%%%%%%
% timing model  %
%%%%%%%%%%%%%%%%%

experiment.onSecs = [zeros(1,experiment.initialFixation)...
    repmat([ones(1,experiment.blockLength) zeros(1,experiment.betweenBlocks) 2*ones(1,experiment.blockLength) zeros(1,experiment.betweenBlocks)],1,experiment.numBlocks/2)...
    zeros(1,experiment.finalFixation)];
experiment.longFormBlocks = Expand(experiment.onSecs,experiment.flipsPerSec,1);
experiment.longFormFlicker = repmat(ones(1,1),1,length(experiment.longFormBlocks));
experiment.waveID = repmat([1:experiment.flipsPerSec],1,length(experiment.longFormBlocks)/experiment.flipsPerSec);
length(experiment.longFormBlocks)

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
% experiment.stimSize = round(experiment.stim.stimSizeDeg*experiment.ppd);                 % in degrees, the size of our objects
% experiment.innerAnnulus = round(experiment.stim.annulusRad*experiment.ppd);
experiment.fixSize = round(experiment.fixSizeDeg*experiment.ppd);
experiment.littleFix = round(experiment.littleFixDeg*experiment.ppd);
% experiment.pixPerCheck = round((experiment.ppd/experiment.stim.spatialFreqDeg)/2);      % half of the pix/cycle

%%%% scale the stims for the screen
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
%% extend stimulus repetition to have the same total number of flips as the whole experiment
experiment.topWaveID = repmat(experiment.topWaveID,1,length(experiment.waveID));
experiment.bottomWaveID = repmat(experiment.bottomWaveID,1,length(experiment.waveID));

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+experiment.vertOffset;

xtop = rect(3)/4;
ytop = rect(4)/4;
experiment.topRect =  CenterRectOnPoint([0 0 experiment.gaborWidth experiment.gaborHeight],xtop,ytop);

xbottom = rect(3)/4;
ybottom = rect(4)/4*3;
experiment.bottomRect =  CenterRectOnPoint([0 0 experiment.gaborWidth experiment.gaborHeight],xbottom,ybottom);


%%%% initial window - wait for backtick
Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(w, 'Flip', 0);

KbTriggerWait(53, deviceNumber);
KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% experiment.targetTimes=[];
% experiment.targets = {};
% experiment.responseTimes=[];
% experiment.responses = {};
% experiment.allColors = {};
% Colors = 'kkrg';
% targetColors = ['rg'];
% lastTargetCounter = 0;
% lastColor = 0;

letters = ['ABCDEFGHIJKLMNOP'];

%%%% character identification task
targetLetters = ['O', 'K'];
experiment.target = [];
experiment.targetTimes = [];
experiment.response = [];
experiment.responseTimes=[];
experiment.accuracy = 0;
experiment.meanRT = 0;

taskText = 'characters';


n=0;
count = 1;
%%%%%%% START task TASK/FLIPPING
% [experiment.totalTime, experiment.allFlips,n]

%[VBLT experiment.startRun FlipTimestamp Missed Beampos] = Screen('Flip', w); % starts timing (experiment.startRun records the stimulus onset time)

while n+1 < length(experiment.allFlips)
    [experiment.longFormBlocks(n+1),experiment.longFormFlicker(n+1)]

    KbQueueStart();
    
    %%%% draw check
    if experiment.longFormBlocks(n+1) == 1 && experiment.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section

       % draw & increment stims

            % top stim
            Screen('DrawTexture', w, experiment.topWaveID(n+1),[],experiment.topRect);
            
            % bottom stim
            %if strcmp(conditions(thisCond).name, 'double-indir') || strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
                Screen('DrawTexture', w, experiment.bottomWaveID(n+1), [], experiment.bottomRect);
            %end
            
            % draw fixation
  %          Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % white fixation ring
            
%             % draw RSVP letter
%             if motionFlip < flipTimes(ceil(length(flipTimes)/2)) %flicker the letter every 0.5s (== 6 flip times if 12 flip per sec)
%                 
%                 [width,height] = RectSize(Screen('TextBounds',w,experiment.letterSequence{n}));
%                 Screen(w, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,experiment.cueColor);
%                     
%             end
%             %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            [VBLT experiment.flipTime(n+1) FlipT missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(n+1)+motionFlip - slack);
             %flipCount = flipCount+1;
 %       end
    
    
%    elseif experiment.longFormBlocks(n+1) == 2 && experiment.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
%         Screen('DrawTexture', w, checkTex{experiment.whichCheck(n+1)},[],checkRect);
%         Screen('DrawTexture',w,apertureSurround);

       % draw & increment stims
%         stimFlipCnt = 0;
%         for motionFlip = flipTimes
%             stimFlipCnt = stimFlipCnt+ 1;
             % top stim
%            Screen('DrawTexture', w, experiment.topWaveID(stimFlipCnt),[],experiment.topRect);
            
            % bottom stim
            %if strcmp(conditions(thisCond).name, 'double-indir') || strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
%               Screen('DrawTexture', w, experiment.bottomWaveID(stimFlipCnt), [], experiment.bottomRect);
            %end
            
            % draw fixation
%             Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % white fixation ring
            
%             % draw RSVP letter
%             if motionFlip < flipTimes(ceil(length(flipTimes)/2)) %flicker the letter every 0.5s (== 6 flip times if 12 flip per sec)
%                 
%                 [width,height] = RectSize(Screen('TextBounds',w,experiment.letterSequence{n}));
%                 Screen(w, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,experiment.cueColor);
%                     
%             end
%             %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             [VBLT experiment.flipTime(n+1) FlipT missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(n+1)+motionFlip - slack);
            % flipCount = flipCount+1;
 %       end
%        Screen('FillOval',w,experiment.backgroundColor,[xc-experiment.innerAnnulus yc-experiment.innerAnnulus xc+experiment.innerAnnulus yc+experiment.innerAnnulus]);
    end
    
    % select new character if starting new trial
    if mod(n, flipsPerTrial) == 0

        % draw character
        l = randperm(numel(letters));
        fixChar = letters(l(1));

       if fixChar == 'O'
            experiment.target = [experiment.target, 1];
            experiment.targetTimes = [experiment.targetTimes, GetSecs - experiment.startRun];
        elseif fixChar == 'K'
            experiment.target = [experiment.target, 2];
            experiment.targetTimes = [experiment.targetTimes, GetSecs - experiment.startRun];
       end


    end
    
    %%%% draw fixation letter in fixation circle
    
    if mod(n, flipsPerTrial) < trialOnFlips
        Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]);
%         DrawFormattedText(w, fixChar, 'center', 8+vertOffset+rect(4)/2,
%         0); %either text function works
        Screen('DrawText', w, fixChar, -5+rect(3)/2, -10+vertOffset+rect(4)/2,[0 0 0]);

    else
        Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]);
        
    end
 
    
       %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if n == 0 
        [VBLT, experiment.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        experiment.flipTime(n+1) = experiment.startRun;
    else
        [VBLT, experiment.flipTime(n+1), FlipT, missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(n+1) - slack);
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
