% LD 7.20.2022
% present two rectangular orientation patches on left side of fixation, with a
% uniform field surrounding them
% RSVP task at fixation throughout
% block design


%clear all;
% Clear the workspace and the screen
sca;
close all;
clear;

Screen('Preference', 'SkipSyncTests', 0);

input('Hit enter to proceed.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% take these out for the actual scan!
offset = 0;
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
deviceNumber = keyboardIndices;
runNum = 1;
subject = 'test';
Screen('Preference', 'SkipSyncTests', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% input this at the beginning of the scan session for 7T
experiment.vertOffset = offset;    % vertical offset from FindScreenSize.m
% experiment.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
% experiment.whichCLUT = '7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
experiment.subject = subject;
experiment.scanNum = runNum; % to keep both of these structs labeled

%%%% scales all of the stimuli in DVA to the screensize
experiment.screenWidth = 19;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
experiment.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;


%%%% set-up rand
rand('twister', sum(100*clock));
experiment.rand = rand;

%%%% files and things
experiment.root = pwd;
experiment.date = datestr(now,30);

%%%% timing
experiment.TRlength = 2;                % in seconds
experiment.betweenBlocks = 16;          % in seconds
experiment.initialFixation = 1;%16;     % in seconds
experiment.finalFixation = 16;          % in seconds
experiment.stimDur = 1; %.2;       % in seconds phase changes
%experiment.targetDur = 0.2;        %in seconds, duration of letter presentation
experiment.ISI = .25;                       % inter stimulus interval in seconds
experiment.IBI = 6;                      % inter block interval in seconds
experiment.blockLength = experiment.stimsPerBlock * (experiment.stimDur + experiment.ISI);
experiment.numBlocks = 8;
experiment.stimsPerBlock = 8;                      % number of sine grating repetitions per block
experiment.trialFreq = 1;               % frequency of fixation trials (seconds)
experiment.trialDur = .4;               % duration in seconds of each fixation trial
experiment.postTargetWindow = 1;           % duration in seconds of target-free interval after a target
experiment.totalTime = experiment.initialFixation + (experiment.numBlocks * (experiment.blockLength + experiment.betweenBlocks)) + experiment.finalFixation;
experiment.totalMins = experiment.totalTime/60;
flipsPerTrial = experiment.trialFreq/experiment.flickerRate;
trialOnFlips = experiment.trialDur/experiment.flickerRate;


%%%% gabor properties
experiment.stim.spatialFreqDeg = 1;                                           % cycles per degree of visual angle
experiment.stim.contrast =  .3;                                                 % in %, maybe??
experiment.stim.orientation = 90;                                                % in degrees
experiment.stim.fromFixation = .6;                                              % in degrees of visual angle
experiment.stim.gaborHDeg = 4;                                                  % in degrees of visual angle
experiment.stim.gaborWDeg = 4; 
experiment.stim.ringPix = 3;                                                    % in pixels, thickness of greyscale ring separating
experiment.stim.contrastMultiplicator = .2;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
experiment.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
experiment.stim.motionRate = 1.3*360 ;                                          % 1.3 cycles per second = 360 deg of phase *1.3 per sec



%%%% conditions & layout
experiment.numMotions = 3;
experiment.motions = {'left','right','none'}; %back and forth starts from the left side, back and forth starts from the right side, no motion if no stimulus on top or bottom
experiment.conds = {'singleTop','singleBottom','double-indir','double-oppdir'};
experiment.numConds = experiment.numMotions * length(experiment.conds);

%%%% Fixation point
experiment.fixSizeDeg =  .4;            % in degrees, the size of the biggest white dot in the fixation
experiment.littleFixDeg = .25;    % proportion of the fixSizeDeg occupied by the smaller black dot
experiment.outerFixPixels = 2;          % in pixels, the black ring around fixation

%experiment.totalTime = experiment.initialFixation+(experiment.numConds*experiment.stimsPerBlock*experiment.blockLength)+((experiment.numConds*experiment.stimsPerBlock-1)*experiment.betweenBlocks)+experiment.finalFixation;

%%%% screen
experiment.backgroundColor = [50 50 50];%[127 127 127];  % color
experiment.fontSize = 20;

%%%% task
experiment.targetProb = .1;              % proportion of trials where the target letters will come up
experiment.firstPossibleTarget = 20;     % no targets in the first X flips
experiment.lastPossibleTarget = 20;      % no targets in the last X flips
experiment.targets = {'J' 'K'};
experiment.distractors = {'A' 'S' 'D' 'F' 'G' 'H' 'L'};
experiment.RSVPrate = 1;               % how fast the letters flip (second) - this likely needs to be hardcoded down below, but is saved here for bookkeeping
experiment.cueColor = 0;%[50 50 255];   % letter color
experiment.totalLetters = (experiment.totalTime/experiment.RSVPrate); %total number of letters that COULD BE presented during the experiment (though based on prob, much less will be presented)
experiment.responseBack = 3;    % the response is correct if the preceding N letters were the target
experiment.flickerRate = .5;         % in seconds then actually in 1 sec the stimuli will change 10 times instead of 20 time due to the way the code was set up
experiment.performance = [];        % will be the performance on the one-back task in each block
experiment.letterFlips = (0:experiment.flickerRate:experiment.totalTime);


%%%% final few structs
experiment.allFlips = (0:experiment.stimDur:experiment.totalTime);

%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2


%%%%%%%%%%%%%%%%%
% timing model  %
%%%%%%%%%%%%%%%%%

%%%% set up our structs (describing all stimulation conditions tested in the
%%%% experiment)
dummyStimTop = (kron([1:experiment.numMotions]', ones(length(experiment.conds),1)))';
dummyStimBottom = repmat([1:experiment.numMotions],[1,length(experiment.conds)]);
dummyDir = repmat([1 1 1 2],[1,experiment.numMotions]);

for n = 1:length(dummyStimTop)
    conditions(n).stimTop = experiment.motions{dummyStimTop(n)};
    conditions(n).stimBottom = experiment.motions{dummyStimBottom(n)};
    conditions(n).dir = dummyDir(n);
    if strcmp(conditions(n).stimTop,conditions(n).stimBottom) ==1 && strcmp(conditions(n).stimTop, 'none') == 0  
        conditions(n).name = {'double-indir'};      
    elseif strcmp(conditions(n).stimTop,conditions(n).stimBottom) ==0 && strcmp(conditions(n).stimTop, 'none') == 0  
        conditions(n).name = {'double-oppdir'};
    elseif strcmp(conditions(n).stimTop,conditions(n).stimBottom) ==0 && strcmp(conditions(n).stimTop, 'none') == 1
        conditions(n).name = {'singleTop'};
    else
        conditions(n).name = {'singleBottom'};
    end
    conditions(n).startTimes = [];
end

 experiment.condShuffle = Shuffle(repmat([1:experiment.numConds],1,experiment.stimsPerBlock));
 %experiment.numBlocks = length(experiment.condShuffle);

%%%% longform condition timing, which aligns with the flicker timing
experiment.longFormConds = zeros(experiment.initialFixation/experiment.stimDur,1);
for n = (1:experiment.numBlocks-1)
    experiment.longFormConds = [experiment.longFormConds; repmat(experiment.condShuffle(n),experiment.blockLength/experiment.stimDur,1)]; % blocks
    experiment.longFormConds = [experiment.longFormConds; zeros(experiment.betweenBlocks/experiment.stimDur,1)]; % inter-block blanks
end
experiment.longFormConds = [experiment.longFormConds; repmat(experiment.condShuffle(experiment.numBlocks),experiment.blockLength/experiment.stimDur,1); zeros(experiment.finalFixation/experiment.stimDur,1)]; % the last block

%% create the timing model for this particular run
counter = experiment.initialFixation; %will assess starting time of each block incrementally, that will be saved in the condition struct
for n=1:experiment.numBlocks
    experiment.startBlock(n) = counter; % timestamp (s) of when each block should start
    conditions(experiment.condShuffle(n)).startTimes = [conditions(experiment.condShuffle(n)).startTimes counter]; % add timestamps to the condition struct
    counter = counter + experiment.blockLength + experiment.betweenBlocks; % and progress the counter
end

%%%%%%%%%%%%%%%%%%%%
% task set-up %
%%%%%%%%%%%%%%%%%%%%

%%%% find the target positioning for this run (index of repeated letters)
experiment.numTargets = length(find(rand(1,experiment.totalLetters)<experiment.targetProb)); %draw targets from uniform distribution between 0 nd 1, only keep the density we need (0.1) of them
targetInd = zeros(1,experiment.totalLetters+1); % this has to go to +1, but it's only for the repeat check
if experiment.numTargets > 0
    while 1
        %check previous (-1,-2) and following indices (+1,+2) if they may constitute a
        %repeated letter
        maybeRep = experiment.firstPossibleTarget+Randi(experiment.totalLetters-experiment.firstPossibleTarget-experiment.lastPossibleTarget);    % a possible index for a repeat, from 3:16
        if targetInd(maybeRep-1) == 0 && targetInd(maybeRep+1) == 0 && targetInd(maybeRep-2) == 0 && targetInd(maybeRep+2) == 0 % if the previous and following TWO images weren't already a repeat
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
        experiment.letterSequence{k} = experiment.targets{randi(length(experiment.targets))}; %pick a random target
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

%%%% set-up phase flicker conditions, and scale up the task accordingly
experiment.longFormFlicker = repmat([1 1]',round((experiment.totalTime/experiment.stimDur)/2),1);
if experiment.RSVPrate < experiment.targetDur
    experiment.letterSequence = expand(experiment.letterSequence,[experiment.RSVPrate/experiment.targetDur,1]); 
    % experiment.letterSequence = expand(experiment.letterSequence,[1, 1/experiment.RSVPrate/experiment.stimDur]);
end



%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

% HideCursor;
Priority(9);

%% open screen
screen=max(Screen('Screens'));
[w, rect]=Screen('OpenWindow',screen,experiment.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat);
Screen(w, 'TextSize', experiment.fontSize);
%Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
% if experiment.gammaCorrect > 0
%     load(experiment.whichCLUT);
%     Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
% end

%% timing optimization
flipInt = Screen('GetFlipInterval',w);
slack = flipInt/2;
frameRate = Screen('NominalFrameRate',w);
experiment.phaseRate = 12; %desired number of grating phases to store per second, which will result in adjusting the flip rate based on the nominal frame rate
flipTimes = [0:flipInt*frameRate/experiment.phaseRate:experiment.stimDur]; %multiply flipInt by 60/12 = 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
experiment.stim.motionPerFlip = experiment.stim.motionRate * flipInt*frameRate/experiment.phaseRate; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip


%%%% scale the stims for the screen
experiment.ppd = pi* rect(3) / (atan(experiment.screenWidth/experiment.viewingDist/2)) / 360; %2pi*(rect(3)/2)= pi*rect(3)
experiment.gaborHeight = round(experiment.stim.gaborHDeg*experiment.ppd);                 % in pixels, the size of our objects
experiment.gaborWidth = round(experiment.stim.gaborWDeg*experiment.ppd);                 % in pixels, the size of our objects
experiment.fromFix = round(experiment.stim.fromFixation*experiment.ppd);
experiment.fixSize = round(experiment.fixSizeDeg*experiment.ppd);

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
        
    topPhase = mod(topPhase + experiment.stim.motionPerFlip,360);
    bottomPhase = mod(bottomPhase + experiment.stim.motionPerFlip,360);
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
        topPhase = mod(topPhase - experiment.stim.motionPerFlip,360);
        bottomPhase = mod(bottomPhase - experiment.stim.motionPerFlip,360);   
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

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+experiment.vertOffset;

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
%%% for the loc
save LGNsurroundParams.mat experiment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% start timing
[vbl experiment.startRun FlipTimestamp Missed Beampos] = Screen('Flip', w); % starts timing

%%%% start recording the response


%[touch,tpress, keyCode]= KbCheck;
KbQueueStart();
experiment.numCorrect = 0;
experiment.correctTarget = [];
flipCount = 1;

%%%%%%% START task TASK/FLIPPING
for n = 1:(length(experiment.allFlips)-1)
    thisCond = experiment.longFormConds(n);
    
    %%%% draw gabors
     if thisCond > 0 && experiment.longFormFlicker(n) > 0% zeros correspond to blanks, in which case we skip this next section
        
        %         if experiment.longFormConds(n-1)==0 % if this is the beginning of a block, we need to randomize the phases
        %             topPhase = randi(360);
        %             if conditions(thisCond).phase == 2 % 1 = oppPhase, 2 = inPhase
        %                 bottomPhase = topPhase;
        %             else
        %                 bottomPhase = mod(topPhase - 180,360); % opp phase for oppPhase
        %             end
        %         end
        
        
        % draw & increment stims
        stimFlipCnt = 0;
        for motionFlip = flipTimes
            stimFlipCnt = stimFlipCnt+ 1;
            % top stim
            Screen('DrawTexture', w, experiment.topWaveID(stimFlipCnt),[],experiment.topRect);
            
            % bottom stim
            if strcmp(conditions(thisCond).name, 'double-indir') || strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
                Screen('DrawTexture', w, experiment.bottomWaveID(stimFlipCnt), [], experiment.bottomRect);
            end
            
            % draw fixation
            Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % white fixation ring
            
            % draw RSVP letter
            if motionFlip < flipTimes(ceil(length(flipTimes)/2)) %flicker the letter every 0.5s (== 6 flip times if 12 flip per sec)
                
                [width,height] = RectSize(Screen('TextBounds',w,experiment.letterSequence{n}));
                Screen(w, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,experiment.cueColor);
                    
            end
            %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [VBLT experiment.flipTime(flipCount) FlipT missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(n)+motionFlip - slack);
            flipCount = flipCount+1;
        end
     else % if phase is 0 or it's blank %need to refine this part for it to stick to ISI
         % draw fixation 
         Screen('FillOval', w,[255 255 255], [xc-round(experiment.fixSize/2) yc-round(experiment.fixSize/2) xc+round(experiment.fixSize/2) yc+round(experiment.fixSize/2)]); % white fixation ring
         % draw RSVP letter
         for motionFlip = flipTimes %need to create the motionFlip variable again, for flicker even in blank conditions
             if motionFlip < flipTimes(ceil(length(flipTimes)/2)) %~nnz(mod(floor(stimFlipCnt/(frameRate/experiment.phaseRate/2)),2))
                 
                 [width,height] = RectSize(Screen('TextBounds',w,experiment.letterSequence{n}));
                 Screen(w, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,experiment.cueColor);
                 
             end
         end
             %         %%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if n == 1 
                 [VBLT experiment.startRun FlipT missed] = Screen(w, 'Flip', 0);
                 experiment.flipTime(flipCount) = experiment.startRun;
             else
                 [VBLT experiment.flipTime(flipCount) FlipT missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(n) - slack);
             end
             
 
         flipCount = flipCount+1;
     end
    % listen for response  - correct if you respond to previous 3 letters
    if n > (experiment.firstPossibleTarget-1) % don't start response listen until targets can appear
        [pressed, firstPress]= KbQueueCheck();
        targetRange = targetInd(n - experiment.responseBack-1:n-2);
        if (pressed ==1) && (sum(targetRange)>0)
            thisTarget = find(targetRange)+n-2-experiment.responseBack; % to convert it back into targetInd's scale
            if strcmp(experiment.letterSequence{thisTarget},'J') == 1 && find(firstPress,1) == KbName('1!')
                experiment.numCorrect = experiment.numCorrect + 1; % then it's correct
                experiment.correctTarget = [experiment.correctTarget thisTarget];
            elseif strcmp(experiment.letterSequence{thisTarget},'K') == 1 && find(firstPress,1) == KbName('2@')
                experiment.numCorrect = experiment.numCorrect + 1;
                experiment.correctTarget = [experiment.correctTarget thisTarget];
            end % then it's correct
        end
        KbQueueFlush();
    end
    
end
%%%% to show the very last flip screen for its 200ms
[VBLT experiment.flipTime(n+1) FlipT missed] = Screen(w, 'Flip', experiment.startRun + experiment.allFlips(length(experiment.allFlips)) - slack);

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

experiment.runTime = GetSecs - experiment.startRun;
% experiment.performance = experiment.numCorrect/experiment.numTargets;

eval(['save data/motion_' experiment.subject '_run' num2str(experiment.scanNum) '_' experiment.date '.mat experiment conditions experiment']);


%KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
%fprintf('Hit rate this run: %.2f%%\n',100*experiment.performance)
fclose all;
clear all;

%end


