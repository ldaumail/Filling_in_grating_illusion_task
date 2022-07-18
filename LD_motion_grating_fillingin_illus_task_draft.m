
% SP 2.1.2015
% present two orientation patches at each side of fixation, with a
% full-field context orientation surrounding them
% RSVP task at fixation throughout
% block design
% VERSION2: uses create procedural sine rather than math to make the
% gratings, hopefully to fix timing issues
% Dec 2015 update: adds another condition: congruent in-phase
% Dec 2015: motion version
% Dec 15 update: instead of random phase in 3rd cond, we do opposite phase
% (old version is _rPhase.m)

%clear all;
% Clear the workspace and the screen
sca;
close all;
clear;

Screen('Preference', 'SkipSyncTests', 0);


input('Hit enter to proceed.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% take these out for the actual scan!
offset = -150;
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
deviceNumber =keyboardIndices;
runNum = 1;
subject = 'test';
Screen('Preference', 'SkipSyncTests', 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% input this at the beginning of the scan session for 7T
params.vertOffset = offset;    % vertical offset from FindScreenSize.m
% params.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
% params.whichCLUT = '7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
experiment.subject = subject;
experiment.scanNum = runNum; % to keep both of these structs labeled

%%%% scales all of the stimuli in DVA to the screensize
params.screenWidth = 19;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
params.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;


%%%% set-up rand
rand('twister', sum(100*clock));
experiment.rand = rand;

%%%% files and things
experiment.root = pwd;
experiment.date = datestr(now,30);

%%%% timing
params.blockLength = 16;            % in seconds
params.betweenBlocks = 16;          % in seconds
params.initialFixation = 16;        % in seconds
params.finalFixation = 16;          % in seconds
params.phaseFlicker = .2;           % in seconds (on for Xs, off for Xs, phase changes

%%%% gabor properties
params.stim.spatialFreqDeg = 1.5;                                           % cycles per degree
params.stim.contrast =  .3;                                                 % in %, maybe??
params.stim.orientation = 90;                                                % in degrees
params.stim.guassianSpaceConstant = .4;                                     % approx equal to the number of radians covered by one standard deviation of the radius of the gaussian mask.
params.stim.fromFixation = .6;                                              % in degrees
params.stim.gaborSizeDeg = 4;                                               % in degrees
params.stim.ringPix = 3;                                                    % in pixels, thickness of greyscale ring separating
params.stim.contrastMultiplicator = .5;                                     % for procedural gabor
params.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
params.stim.motionRate = 360*5;                                             % in degrees/second

%%%% conditions & layout
params.numMotions = 2;
params.motions = {'left','right'};
params.conds = {'single','double-inPhase','double-oppPhase'};
params.numConds = params.numMotions * length(params.conds);
params.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation
params.littleFixDeg = params.fixSizeDeg* .35;    % proportion of the fixSizeDeg occupied by the smaller black dot
params.outerFixPixels = 2;          % in pixels, the black ring around fixation
params.TRlength = 2;                % in seconds
params.repsPerRun = 1;              % repetitions of each object type x location
experiment.totalTime = params.initialFixation+(params.numConds*params.repsPerRun*params.blockLength)+((params.numConds*params.repsPerRun-1)*params.betweenBlocks)+params.finalFixation;
experiment.totalMins = experiment.totalTime/60;

%%%% screen
params.backgroundColor = [127 127 127];  % color
params.fontSize = 20;


experiment.allFlips = (0:params.phaseFlicker:experiment.totalTime);

%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2


%%%%%%%%%%%%%%%%%
% timing model  %
%%%%%%%%%%%%%%%%%

%%%% set up our structs (describing all stimulation conditions tested in the
%%%% experiment)
dummyStimTop = (kron([1:params.numMotions]', ones(length(params.conds),1)))';
dummyStimBottom = repmat([1:params.numMotions],[1,length(params.conds)]);
dummyPhase = repmat([1 1 2],[1,params.numMotions]);

for n = 1:length(dummyStimTop)
    conditions(n).stimTop = params.motions{dummyStimTop(n)};
    conditions(n).stimBottom = params.motions{dummyStimBottom(n)};
    conditions(n).phase = dummyPhase(n);
     if strcmp(conditions(n).stimTop,conditions(n).stimBottom) ==1
         if conditions(n).phase == 1
             conditions(n).name = {'double-oppPhase'};
         elseif conditions(n).phase == 2
             conditions(n).name = {'double-inPhase'}; end
     else
         conditions(n).name = {'single'}; end
    conditions(n).startTimes = [];
end

 experiment.condShuffle = Shuffle(repmat([1:params.numConds],1,params.repsPerRun));
 experiment.numBlocks = length(experiment.condShuffle);

%%%% longform condition timing, which aligns with the flicker timing
experiment.longFormConds = zeros(params.initialFixation/params.phaseFlicker,1);
for n = (1:experiment.numBlocks-1)
    experiment.longFormConds = [experiment.longFormConds; repmat(experiment.condShuffle(n),params.blockLength/params.phaseFlicker,1)]; % blocks
    experiment.longFormConds = [experiment.longFormConds; zeros(params.betweenBlocks/params.phaseFlicker,1)]; % inter-block blanks
end
experiment.longFormConds = [experiment.longFormConds; repmat(experiment.condShuffle(experiment.numBlocks),params.blockLength/params.phaseFlicker,1); zeros(params.finalFixation/params.phaseFlicker,1)]; % the last block

%%%% create the timing model for this particular run
counter = params.initialFixation;
for n=1:experiment.numBlocks
    experiment.startBlock(n) = counter; % timestamp (s) of when each block should start
    conditions(experiment.condShuffle(n)).startTimes = [conditions(experiment.condShuffle(n)).startTimes counter]; % add timestamps to the condition struct
    counter = counter + params.blockLength + params.betweenBlocks; % and progress the counter
end


%%%% set-up phase flicker conditions, and scale up the task accordingly
experiment.longFormFlicker = repmat([1 1]',round((experiment.totalTime/params.phaseFlicker)/2),1);
% if params.RSVPrate < params.phaseFlicker
%     experiment.letterSequence = expand(experiment.letterSequence,params.RSVPrate/params.phaseFlicker,1);
% end
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%

HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
[win, rect]=Screen('OpenWindow',screen,params.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat);
Screen(win, 'TextSize', params.fontSize);
%Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% gamma correction
% if params.gammaCorrect > 0
%     load(params.whichCLUT);
%     Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
% end

%%%% timing optimization
flipInt = Screen('GetFlipInterval',win);
slack = flipInt/2;
params.stim.motionPerFlip = params.stim.motionRate * flipInt;
flipTimes = [0:flipInt:params.phaseFlicker];
flipTimes = flipTimes(1:length(flipTimes)-1);


%%%% scale the stims for the screen
params.ppd = pi* rect(3) / (atan(params.screenWidth/params.viewingDist/2)) / 360;
params.freq =  (params.stim.spatialFreqDeg)*2*pi/params.ppd;
params.gaborSize = round(params.stim.gaborSizeDeg*params.ppd);                 % in degrees, the size of our objects
params.fromFix = round(params.stim.fromFixation*params.ppd);
params.fixSize = round(params.fixSizeDeg*params.ppd);
params.littleFix = round(params.littleFixDeg*params.ppd);


xc = rect(3)/2; % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+params.vertOffset;


% create sine wave gratings
[topWave, topRect] = CreateProceduralSineGrating(win, round(params.gaborSize/2), round(params.gaborSize/2), params.stim.contrastOffset, [], params.stim.contrastMultiplicator);
[bottomWave, bottomRect] = CreateProceduralSineGrating(win, round(params.gaborSize/2), round(params.gaborSize/2), params.stim.contrastOffset, [], params.stim.contrastMultiplicator);

% rects to draw the circle outline
%params.leftRect = CenterRectOnPoint(centerRect,(xc-params.fromFix-floor(params.gaborSize/2)),yc);
%params.rightRect = CenterRectOnPoint(centerRect,(xc+params.fromFix+floor(params.gaborSize/2)),yc);

% critical for the phase stuff! must account for the fact that we can
% offset the center for the scan subject
%surroundRect = CenterRectOnPoint(surroundRect,xc,yc);

%params.leftCircle = [(xc-params.fromFix-params.gaborSize) (yc-round(params.gaborSize/2)) (xc-params.fromFix) (yc+round(params.gaborSize/2))];
%params.rightCircle = [(xc+params.fromFix) (yc-round(params.gaborSize/2)) (xc+params.fromFix+params.gaborSize) (yc+round(params.gaborSize/2))];


%%%% initial window - wait for backtick
Screen(win, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(win, 'Flip', 0);

% KbTriggerWait(53, deviceNumber);
% KbQueueCreate(deviceNumber,responseKeys);

%%% for the loc
save LGNsurroundParams.mat params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% start recording the response
KbQueueStart();
experiment.numCorrect = 0;
experiment.correctTarget = [];
flipCount = 1;

%%%%%%% START task TASK/FLIPPING
for n = 1:(length(experiment.allFlips)-1)
    thisCond = experiment.longFormConds(n);
    
    %%%% draw gabors
    if thisCond > 0 && experiment.longFormFlicker(n) > 0% zeros correspond to blanks, in which case we skip this next section
        
        if experiment.longFormConds(n-1)==0 % if this is the beginning of a block, we need to randomize the phases
            topPhase = randi(360);
            if conditions(thisCond).phase == 2 % 1 = oppPhase, 2 = inPhase
                bottomPhase = topPhase; 
            else
                bottomPhase = mod(topPhase - 180,360); end % opp phase for incongruent and oppPhase
        end
        
        % surroundPhase = 0; centerPhase = 0;
        
        % draw & increment stims
        for motionFlip = flipTimes
            % top stim
            dstRect = OffsetRect(topRect, rect(3)/2-topRect(3)/2, rect(3)/4);
            Screen('DrawTexture', win, topWave, [], dstRect, params.stim.orientation, [], [], [], [], [], [topPhase, params.stim.spatialFreqDeg/params.ppd, params.stim.contrast, 0]);
  
            % bottom stim
            if conditions(thisCond).phase < 2; % draw center if it's randomized or incognruent
            Screen('DrawTexture', win, bottomWave, [], [], params.stim.orientation, [], [],...
                [], [], [], [bottomPhase, params.stim.spatialFreqDeg/params.ppd, params.stim.contrast, 0]);
            
            Screen('DrawTexture', win, bottomWave, [], [], params.stim.orientation, [], [],...
                [], [], [], [bottomPhase, params.stim.spatialFreqDeg/params.ppd, params.stim.contrast, 0]);
            end
            % Draws the grating 'gratingid' into window 'windowPtr', at position 'dstRect'
            % or in the center if dstRect is set to []. Make sure 'dstRect' has the
            % size of 'gratingrect' to avoid spatial distortions! You could do, e.g.,
            % dstRect = OffsetRect(gratingrect, xc, yc) to place the grating centered at
            % screen position (xc,yc). 'Angle' is the optional orientation angle,
            % default is zero degrees. 'modulateColor' is the base color of the grating
            % patch - it defaults to white, ie. the grating has only luminance, but no
            % color. If you'd set it to [255 0 0] you'd get a reddish grating. 'phase' is
            % the phase of the grating in degrees, 'freq' is its spatial frequency in
            % cycles per pixel, 'contrast' is the contrast of your grating.
            
            %Screen('FrameOval',win,params.backgroundColor,params.leftCircle,[],[]);
            %Screen('FrameOval',win,params.backgroundColor,params.rightCircle,[],[]);
            
            
            
            % draw fixation and RSVP letter
            Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
%             [width,height] = RectSize(Screen('TextBounds',win,experiment.letterSequence{n}));
%             Screen(win, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,params.cueColor);
%             
            %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [VBLT experiment.flipTime(flipCount) FlipT missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(n)+motionFlip - slack);
            flipCount = flipCount+1;
            
            %%%% increment phase to show motion
            if strcmp(conditions(thisCond).name,'double') ==1
                if strcmp(conditions(thisCond).stimTop,'left') == 1
                    topPhase = mod(topPhase + params.stim.motionPerFlip,360);
                    bottomPhase = mod(bottomPhase - params.stim.motionPerFlip,360);
                elseif strcmp(conditions(thisCond).stimTop,'right') == 1
                    topPhase = mod(topPhase - params.stim.motionPerFlip,360);
                    bottomPhase = mod(bottomPhase + params.stim.motionPerFlip,360);
                end
            else % single stimulus condition
                if strcmp(conditions(thisCond).stimTop,'left') == 1
                    topPhase = mod(topPhase + params.stim.motionPerFlip,360);
                   % bottomPhase = mod(bottomPhase + params.stim.motionPerFlip,360);
                elseif strcmp(conditions(thisCond).stimTop,'right') == 1
                    topPhase = mod(topPhase - params.stim.motionPerFlip,360);
                   % bottomPhase = mod(bottomPhase - params.stim.motionPerFlip,360);
                end
            end
        end
    else % if phase is 0 or it's blank
        % draw fixation and RSVP letter
        Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
%         [width,height] = RectSize(Screen('TextBounds',win,experiment.letterSequence{n}));
%         Screen(win, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,params.cueColor);
%         %%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 1 [VBLT experiment.startRun FlipT missed] = Screen(win, 'Flip', 0);
            experiment.flipTime(flipCount) = experiment.startRun;
        else [VBLT experiment.flipTime(flipCount) FlipT missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(n) - slack);end
        flipCount = flipCount+1;
     end
    

end
%%%% to show the very last flip screen for its 200ms
[VBLT experiment.flipTime(n+1) FlipT missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(length(experiment.allFlips)) - slack);

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

experiment.runTime = GetSecs - experiment.startRun;
% experiment.performance = experiment.numCorrect/experiment.numTargets;

eval(['save data/motion_' experiment.subject '_run' num2str(experiment.scanNum) '_' experiment.date '.mat params conditions experiment']);


KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fprintf('Hit rate this run: %.2f%%\n',100*experiment.performance)
fclose all;
clear all;

%end


