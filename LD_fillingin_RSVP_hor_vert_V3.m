%function figureGround_loc_v5(subject, session, vertOffset, debug) 
subject = 1;
session = 1;
vertOffset = 0;
%%%% resolution
% %if debug == 1
%     e.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
%     e.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;
% 	e.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop
% else
    e.screenWidth = 16;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    e.viewingDist = 23;             % in cm; %23 in eye tracking room 425 3Tb/office=43, miniHelm=57;
    e.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % scanner
%end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2

Screen('Preference', 'SkipSyncTests', 0);

e.scanNum = input('Scan number :');
e.runNum = input('Run number :');
e.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m
e.gammaCorrect = 1;       % make sure this = 1 when you're at the scanner!
e.whichCLUT = '7T_Sam.mat'; %'linearizedCLUT_SoniaMPB.mat';

%%% basic naming set-up
e.subject = subject;
e.session = session;


%%%% set-up rand
 rand('twister', sum(100*clock));
 e.rand = rand;

%rng(sum(100*clock));
%e.rand = rng;
%%%% files and things
e.root = pwd;
e.date = datestr(now,30);

%%%% sine wave grating timing

e.initialFixation = 6;        % in seconds
e.finalFixation = 0;          % in seconds
e.stimDur = 1.6*2;        % in seconds. 1.6 sec refers to sine wave grating 1.6 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
e.stimsPerBlock = 5;
e.blockLength = e.stimDur*e.stimsPerBlock;            % in seconds
e.betweenBlocks = 16;          % in seconds
e.flipsPerSec = 12;           % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
e.flipWin = 1/e.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%e.numBlocks = 12;  % 6 for single and 6 for pair...

%%%% 2D sine wave grating properties
e.stim.spatialFreqDeg = 0.286;  %1                                         % cycles per degree of visual angle
e.stim.contrast = 0.3 ;                                                 % in %, maybe??
e.stim.orientation = [90 180];                                                % in degrees
e.stim.degFromFix = .6;                                              % in degrees of visual angle
e.stim.gaborHDeg = 7;                                                  % in degrees of visual angle
e.stim.gaborWDeg = 7; 
e.stim.contrastMultiplicator = .075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
e.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
e.stim.motionRate = 1.13*360 ;  %1.3                                        % 1.3 cycles per second = 360 deg of phase *1.3 per sec

%%%% conditions & layout
e.conds = {'vertSingleTop','horSingleTop','vertSingleBottom','horSingleBottom','vertDoubleIndir','horDoubleIndir'}; %'double-oppdir'
e.numConds = length(e.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
e.condShuffle = Shuffle(repmat([1:e.numConds],1,1)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order
e.numBlocks = length(e.condShuffle)*2;
e.fixSizeDeg =  .6;            % in degrees, the size of the biggest white dot in the fixation
e.repsPerRun = 2;              % repetitions of each object type x eation
e.totalTime = e.initialFixation + ((e.numBlocks-1) * (e.blockLength + e.betweenBlocks)) + e.blockLength + e.finalFixation;
e.allFlips = (0:e.flipWin:e.totalTime);

%%%% screen
e.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
e.fontSize = 12; %26;

%%%% RSVP task setup
e.targetProb = .1;              % proportion of trials where the target letters will come up
e.firstPossibleTarget = 20;     % no targets in the first X flips
e.lastPossibleTarget = 20;      % no targets in the last X flips
e.targetLetters = {'J' 'K'};
e.distractors = {'A' 'S' 'D' 'F' 'G' 'H' 'L'};
e.trialFreq = 0.5;               % duration of fixation trials (seconds) (time it takes to switch from one letter to blank then back to one letter)
flipsPerTrial = e.trialFreq/e.flipWin;
e.trialDur = e.trialFreq/2;               % duration in seconds of the letter presentation of each fixation trial (.25s letter ON, .25 letter OFF)
e.trialOnFlips = floor(e.trialDur/e.flipWin); %number of flips the letter is on over the course of 1 trial
%e.cueColor = 0;%[50 50 255];   % letter color
e.totalLetters = (e.totalTime/e.trialFreq); %total number of letters that COULD BE presented during the e (though based on prob, much less will be presented)
e.responseBack = 3;    % the response is correct if the preceding N letters were the target
%e.letterTiming = (0:e.RSVPrate:e.totalTime);

%%%% character identification task
e.targets = [];
e.targetTimes = [];
e.responses = [];
e.responseTimes=[];
e.accuracy = 0;
e.meanRT = 0;

taskText = 'characters';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ET PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ET_ON = 1;
% 
% % This applies to actual trials
% e.AbortTrialWhenExceeded = [ ] ; % During actual expt
% 
% % This applies to practice trials (turned off during actual expt)
% e.ShowETbox = [ 1  ]; % during practice
% e.ETbox = [  3  ];% [ 3+.5 ]; % in VA; actually, diameter of a circle
% e.ETbox_size_pix = [ e.ETbox e.ETbox ] * e.ppd; 
% e.ETbox_color = [ 255 255 255 ]*.6;
% % e.ShowRealTimeGaze = [ 1 ]; % [] or [ something ]
% e.ShowRealTimeGaze = [  ]; % [] or [ something ]
% e.nGazetoShow = [ 60 ]; % current~past N fixations
% e.nGazeAllowed =  [ 12 ];% [ 12 ]=1st session of AP,DM; % Abort the trial if recent N trials were ALL outside the circle
% cmap_gaze = [255 0 0]; % repmat([255 0 0]', 1, s.nGazetoShow); 
% cmap_gaze_exceed = [255 0 0]; % repmat([0 0 255]', 1, s.nGazetoShow);    
% 
% e.Fixation.Color = [1 1 1]*255*0; % [1 1 1]*[ 128*.4 ]; % [ s.Background*.7 ]; % Inside Fixation(Cross) s.Black(1); %
% e.Fixation.Color2 = [ e.Fixation.Color ]; % [1 1 1]*[ 10] ; % [ 255 0 0 ]; % [0 0 0]
% % tmp_VA = [ .25 ];% .25? radius = .5? in diameter  %% [.32 ];  % [2, 6, 8]; in pixel; of Radius; was [2, 4, 6] in pilot; 
% tmp_VA = [ .1 ];% .25? radius = .5? in diameter  %% [.32 ];  % [2, 6, 8]; in pixel; of Radius; was [2, 4, 6] in pilot; 
% e.Fixation.R_pix = round(tmp_VA * e.ppd); 
% e.Fixation.R_VA = e.Fixation.R_pix / e.ppd;
% e.Fixation.R2_pix = round(e.Fixation.R_VA*[ .2 ] * e.ppd); 
% e.Fixation.R2_VA =  e.Fixation.R2_pix/e.ppd;  % [2, 6, 8]; in pixel; of Radius; was [2, 4, 6] in pilot;     
% e.Fixation.Width_pix = round(e.Fixation.R_pix*.2 ); 
% e.Fixation.Width_VA = e.Fixation.Width_pix / e.ppd; %  [ .03 ]; % Pen Width   
% 
% e.Fixation.rect  = CenterRect( [0 0 e.Fixation.R_pix*2  e.Fixation.R_pix*2], rect );
% e.Fixation.rect2 = CenterRect( [0 0 e.Fixation.R2_pix*2 e.Fixation.R2_pix*2], rect );

%%
%%%% INITIALIZE EYETRACKING
if ~isempty(ET_ON)
    if EyelinkInit()~= 1;
        error('ERROR: Eyetracking not on-line.  Make sure it is plugged in and turned on...\n')
        sca
        return;
    end;

    %%% Sets ET connection
    el = EyelinkInitDefaults(win);

    %%% Custum calibration: Adapted from EyelinkPictureCustomCalibration
    el.calibrationtargetsize=1; % from 2.5
    el.calibrationtargetwidth=.5; % from 1
    el.displayCalResults = 1; % [ avgError ?? ??? ]
    EyelinkUpdateDefaults(el)
    
    % SET UP TRACKER CONFIGURATION
    % Setting the proper recording resolution, proper calibration type,
    % as well as the data file content;
    width=e.resolution.width; height=e.resolution.height;
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
    Eyelink('command', 'calibration_type = HV13'); % HV5 HV9 HV13 
    Eyelink('command', 'generate_default_targets = NO');
    
    %%%% layout of calibration points 
    ncal_angle = 8; % 8 angles along the larger circle (one per 45?)
    cal_angles = linspace(0, 2*pi, ncal_angle+1)+pi/4; cal_angles(end)=[];
    cal_radius = [ 5 ]*e.ppd; % Eccentricity    
    ncal_angle2 = 4; % 4 angles along the smaller circle (one per 45?)
    cal_angles2 = linspace(0, 2*pi, ncal_angle2+1); cal_angles2(end)=[];
    cal_radius2 = cal_radius/2; % [ 2.5 ]*s.ppd; 
    
    cal_xx = round(  width/2 + [cos(cal_angles)*cal_radius  cos(cal_angles2)*cal_radius2 ]);
    cal_yy = round( height/2 + [sin(cal_angles)*cal_radius  sin(cal_angles2)*cal_radius2 ]);    
    cal_xx = [ width/2 cal_xx ];
    cal_yy = [ height/2 cal_yy ];         
    
    ncal_sample = length(cal_xx)  % 4+8+1(center) = 13
    cal_sequence = sprintf('%d,',[ 0 1:ncal_sample] ); cal_sequence(end)=[]
    cal_targets = [ cal_xx; cal_yy ]; cal_targets = cal_targets(:)';
    cal_targets = sprintf('%d,%d ', cal_targets)
    
    % STEP 5.1 modify calibration and validation target locations
    Eyelink('command', ['calibration_samples = ' num2str(ncal_sample+1) ]); % +1 = including initial 0(=center)
    Eyelink('command', ['calibration_sequence = ' cal_sequence  ] );% first xy = fixation (=sequence 1)
    Eyelink('command', ['calibration_targets = ' cal_targets ]);
    Eyelink('command', ['validation_samples = ' num2str(ncal_sample) ]); 
    Eyelink('command', ['validation_sequence = ' cal_sequence  ] );
    Eyelink('command', ['validation_targets = ' cal_targets ]);

    %%% Tells the ET what data to record
    Eyelink('Command', 'file_sample_data = LEFT,RIGHT,GAZE,DIAMETER');
    Eyelink('Command', 'file_event_data = GAZE,GAZEREZ,DIAMETER,HREF,VELOCITY');
    Eyelink('Command', 'file_event_filter = LEFT, RIGHT, FIXATION,SACCADE, BLINK, MESSAGE, BUTTON');
end 
   
if ~isempty(ET_ON)
    EyelinkDoTrackerSetup(el);
end

%%%% Drift correction:
if ~isempty(ET_ON)

    ETbox_rect = CenterRectOnPoint([0 0 e.ETbox_size_pix], centerX, centerY);

    % Open a file to record to
    tmpETfile = 'ET_tmp.edf';
    Eyelink('openfile', tmpETfile);

    %%%% Drift correction:
    success = 0;
    while success == 0
        
        % In fact, there's an internal while loop inside the eyelink
        % driftcorrect code. Still, we need repeat driftcorrection until it
        % succeeds. Without the loop, the trial can proceed with the ESC key
        % press despite a poor fixation

        if ~isempty(e.ShowETbox)
            Screen('FrameOval', win, e.ETbox_color, ETbox_rect) ;%  [,penWidth]);
        end
        Screen('FrameOval', win, [0 0 0], e.Fixation.rect, e.Fixation.Width_pix)
        Screen('FillOval', win, [255 0 0], e.Fixation.rect2 ); %  , s.Fixation.Width_pix
        Screen('Flip',win);
        success = EyelinkDoDriftCorrection(el, [], [], 0, 0); % You can't get out of this function unless there's a success or ESC
        
    end 
    Eyelink('StartRecording');
    eye_used = Eyelink('EyeAvailable');
    
end 
        
ETdur_per_trial = e.trialLength + [ 5 ];
FrameRate = Screen('FrameRate', win); %  hz=Screen('FrameRate', w);
max_ET_nframe = FrameRate * ETdur_per_trial;
gazesamples = cell(length(blockNum), numTrials);
for i = 1:length(blockNum)
    for j = 1:numTrials
        gazesamples{i,j}.x = nan(max_ET_nframe,1);
        gazesamples{i,j}.y = nan(max_ET_nframe,1);
        gazesamples{i,j}.pa = nan(max_ET_nframe,1);
        gazesamples{i,j}.TS = nan(max_ET_nframe,1);
        gazesamples{i,j}.D = nan(max_ET_nframe,1);
    end
end


%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

e.onSecs = [zeros(1,e.initialFixation)...
    repmat([ones(1,e.blockLength) zeros(1,e.betweenBlocks)],1,e.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
    ones(1,e.blockLength) zeros(1,e.finalFixation)];
e.longFormBlocks = Expand(e.onSecs,e.flipsPerSec,1); %1 when block, 0 when between block
e.longFormFlicker = repmat(ones(1,1),1,length(e.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(e.longFormBlocks)


%set up the timing model for stimulus pairing conditions (single (top,
%bottom), pair) and stimulus orientation
%longform condition timing, which aligns with the flicker timing
e.longFormConds = zeros(1,e.initialFixation);
for n = 1:e.numBlocks/length(e.condShuffle)
    if n < e.numBlocks/length(e.condShuffle)
        len = length(e.condShuffle);
    else
        len = length(e.condShuffle)-1;
    end
    for i = (1:len)
        e.longFormConds = [e.longFormConds, repmat(e.condShuffle(i),1,e.blockLength)]; % blocks
        e.longFormConds = [e.longFormConds, zeros(1,e.betweenBlocks)]; % inter-block blanks
    end
end
e.longFormConds = [e.longFormConds, repmat(e.condShuffle(end),1,e.blockLength), zeros(1,e.finalFixation)]; % the last block
e.longFormConds = Expand(e.longFormConds, e.flipsPerSec,1);
length(e.longFormConds)
% %% create the timing model of stimulus conditions for this particular run
clear i
for i =1:e.numConds
    conditions(i).name = e.conds(i);
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
[w, rect]=Screen('OpenWindow',screen,e.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
%[w, rect]=Screen('OpenWindow',screen,e.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
Screen(w, 'TextSize', e.fontSize);
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
flipTimes = [0:frameInt*frameRate/e.flipsPerSec:e.stimDur]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
e.stim.dphasePerFlip = e.stim.motionRate*frameInt * frameRate/e.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip

%%%% scale the stims for the screen
e.ppd = pi* rect(3) / (atan(e.screenWidth/e.viewingDist/2)) / 360;
e.fixSize = round(e.fixSizeDeg*e.ppd);
e.gaborHeight = round(e.stim.gaborHDeg*e.ppd);                 % in pixels, the size of our objects
e.gaborWidth = round(e.stim.gaborWDeg*e.ppd);                 % in pixels, the size of our objects

xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2+e.vertOffset;

%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers
topPhase = randi(360);
bottomPhase = topPhase;
e.topWave = nan(length(flipTimes),e.gaborHeight,e.gaborWidth,length(e.stim.orientation));
e.bottomWave = nan(length(flipTimes),e.gaborHeight,e.gaborWidth,length(e.stim.orientation));
e.topWaveID = nan(length(flipTimes)*length(e.longFormFlicker),length(e.stim.orientation));
e.bottomWaveID = nan(length(flipTimes)*length(e.longFormFlicker),length(e.stim.orientation));

for o =1:length(e.stim.orientation)
    for f = 1:length(flipTimes)
        
        if f <= length(flipTimes)/2
            
            topPhase = mod(topPhase + e.stim.dphasePerFlip,360);
            bottomPhase = mod(bottomPhase + e.stim.dphasePerFlip,360);
            %ih in pixels %iw in pixels %spatial freq in cycles per dva
            %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
            %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
            %background color (unused if the grating is not an annulus)
            
            e.topWave(f,:,:,o) = makeSineGrating(e.gaborHeight,e.gaborWidth,e.stim.spatialFreqDeg,...
                e.stim.orientation(o),topPhase,e.stim.contrastOffset(1),e.stim.contrastMultiplicator,...
                e.ppd);
            e.bottomWave(f,:,:,o) = makeSineGrating(e.gaborHeight,e.gaborWidth,e.stim.spatialFreqDeg,...
                e.stim.orientation(o),bottomPhase,e.stim.contrastOffset(1),e.stim.contrastMultiplicator,...
                e.ppd);
            %    figure();
            %   imshow(squeeze(topWave(f,:,:)));
            %
        elseif f > length(flipTimes)/2
            topPhase = mod(topPhase - e.stim.dphasePerFlip,360);
            bottomPhase = mod(bottomPhase - e.stim.dphasePerFlip,360);
            e.topWave(f,:,:,o) = makeSineGrating(e.gaborHeight,e.gaborWidth,e.stim.spatialFreqDeg,...
                e.stim.orientation(o),topPhase,e.stim.contrastOffset(1),e.stim.contrastMultiplicator,...
                e.ppd);
            e.bottomWave(f,:,:,o) = makeSineGrating(e.gaborHeight,e.gaborWidth,e.stim.spatialFreqDeg,...
                e.stim.orientation(o),bottomPhase,e.stim.contrastOffset(1),e.stim.contrastMultiplicator,...
                e.ppd);
        end
        tmptopWaveID(f,o) = Screen('MakeTexture', w, squeeze(e.topWave(f,:,:,o)));
        tmpbottomWaveID(f,o) = Screen('MakeTexture', w, squeeze(e.bottomWave(f,:,:,o)));
        
    end
    
%% extend stimulus matrix to include the same total number of flips as the whole experiment
e.topWaveID(:,o) = repmat(tmptopWaveID(:,o),length(e.longFormFlicker),1);
e.bottomWaveID(:,o) = repmat(tmpbottomWaveID(:,o),length(e.longFormFlicker),1);

end

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xtop = rect(3)/2 - (8+3.5)*e.ppd;%rect(3)/4; % = stimulus center located 8 degrees horizontal from the center
ytop = rect(4)/2 - e.gaborHeight;
e.topRect =  CenterRectOnPoint([0 0 e.gaborWidth e.gaborHeight],xtop,ytop);

xbottom = rect(3)/2 - (8+3.5)*e.ppd; %rect(3)/4; = stimulus center located 8 degrees horizontal from the center
ybottom = rect(4)/2 + e.gaborHeight;%rect(4)/4*3;
e.bottomRect =  CenterRectOnPoint([0 0 e.gaborWidth e.gaborHeight],xbottom,ybottom);


%% %%%%%%%%%%%%%%%%%%%%%%
   % Letter task set-up %
   %%%%%%%%%%%%%%%%%%%%%%

%%%% find the target positioning for this run (index of repeated letters)
e.numTargets = length(find(rand(1,e.totalLetters)<e.targetProb)); %draw targets from uniform distribution between 0 nd 1, only keep the density we need (0.1) of them
targetInd = zeros(1,e.totalLetters+1); % this has to go to +1, but it's only for the repeat check
if e.numTargets > 0
    while 1
        %check previous (-1,-2) and following indices (+1,+2) if they might constitute a
        %repeated letter
        maybeRep = e.firstPossibleTarget+Randi(e.totalLetters-e.firstPossibleTarget-e.lastPossibleTarget);    % a possible index for a target letter, picked randomly
        if targetInd(maybeRep-1) == 0 && targetInd(maybeRep+1) == 0 && targetInd(maybeRep-2) == 0 && targetInd(maybeRep+2) == 0 % make sure the previous  TWO and following TWO letters weren't already a repeat
            targetInd(maybeRep) = 1; %record this index if it is not a repeat
        end
        if sum(targetInd) == e.numTargets %if the sum of target letters reach the number of targets, stop the check
            break
        end
    end
end
targetInd = targetInd(1:e.totalLetters); % trim this back for sanity


%%%% fill the sequence out with actual letters
e.letterSequence = [];
for k = 1:length(targetInd)
    if targetInd(k) == 1 % targets - we know these don't repeat
        e.letterSequence{k} = e.targetLetters{randi(length(e.targetLetters))}; %pick a random target
    else % distractors - we need to make sure these don't repeat
        if k > 1 % the first letter can be whatever
            while 1
                maybeLetter = e.distractors{randi(length(e.distractors))};
                if strcmp(maybeLetter, e.letterSequence{k-1}) == 0
                    e.letterSequence{k} = maybeLetter;
                    break
                end
            end
        else %% (the first letter can be whatever)
            e.letterSequence{k} = e.distractors{randi(length(e.distractors))};
        end
    end
end

% %%%% scale up the task based on flips
e.letterSequence = Expand(e.letterSequence,flipsPerTrial,1);%e.flipsPerSec,1); 


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
% [e.totalTime, e.allFlips,n]
gray = repmat(min(min(squeeze(e.topWave(1,:,:)),[],1)), [1,3]);
Screen('FillRect', w, gray);
while n+1 < length(e.allFlips)
    
    [e.longFormBlocks(n+1),e.longFormFlicker(n+1)]
    thisCond = e.longFormConds(n+1);
    tic
    KbQueueStart();
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 0
        [VBLT, e.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        e.flipTime(n+1) = e.startRun;
    else
        [VBLT, e.flipTime(n+1), FlipT, missed] = Screen(w, 'Flip', e.startRun + e.allFlips(n+1) - slack);
    end
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if e.longFormBlocks(n+1) == 1 && e.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        if strfind(conditions(thisCond).name{:}, 'vert')
            ori =1;
        else
            ori = 2;
        end
        % draw & increment stims
        if strfind(conditions(thisCond).name{:}, 'DoubleIndir')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            % top stim
            Screen('DrawTexture', w, e.topWaveID(n+1,ori),[],e.topRect);
            % bottom stim
            Screen('DrawTexture', w, e.bottomWaveID(n+1,ori), [], e.bottomRect);
        elseif strfind(conditions(thisCond).name{:}, 'SingleTop')
            % top stim
            Screen('DrawTexture', w, e.topWaveID(n+1,ori),[],e.topRect);
        elseif strfind(conditions(thisCond).name{:}, 'SingleBottom')
            % bottom stim
            Screen('DrawTexture', w, e.bottomWaveID(n+1,ori), [], e.bottomRect);
        end
        
    end
    
    % select new character if starting new trial
    if mod(n, flipsPerTrial) == 0

%         % draw character
%         l = randperm(numel(letters));
%         fixChar = letters(l(1));
          fixChar = e.letterSequence(n+1);

       if strcmp(fixChar,'J') == 1
            e.targets = [e.targets, 1];
            e.targetTimes = [e.targetTimes, GetSecs - e.startRun];
        elseif strcmp(fixChar,'K') == 1
            e.targets = [e.targets, 2];
            e.targetTimes = [e.targetTimes, GetSecs - e.startRun];
       end


    end
    
    %%%% draw fixation letter in fixation circle
    
    if mod(n, flipsPerTrial) <= e.trialOnFlips
        Screen('FillOval', w,[255 255 255], [xc-round(e.fixSize/2) yc-round(e.fixSize/2) xc+round(e.fixSize/2) yc+round(e.fixSize/2)]);%white fixation solid circle
        %DrawFormattedText(w, fixChar, 'center', 8+vertOffset+rect(4)/2,0); %either text function works
        %Screen('DrawText', w, fixChar, -5+rect(3)/2, -10+vertOffset+rect(4)/2,[0 0 0]);
        [width,height] = RectSize(Screen('TextBounds',w,fixChar{:}));
        Screen('DrawText', w, fixChar{:}, xc-width/2 +0.3, yc-height/2-0.2,[0 0 0]);
     else
        Screen('FillOval', w,[255 255 255], [xc-round(e.fixSize/2) yc-round(e.fixSize/2) xc+round(e.fixSize/2) yc+round(e.fixSize/2)]);%white fixation solid circle
   
    end

%     n

    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    
    %%%% character identification
    if (pressed == 1) && ((firstPress(KbName('1!')) > 0) || (firstPress(KbName('2@')) > 0))
        if firstPress(KbName('1!')) > 0
            e.responses = [e.responses, 1];
            e.responseTimes = [e.responseTimes, firstPress(KbName('1!')) - e.startRun];
        elseif firstPress(KbName('2@')) > 0
            e.responses = [e.responses, 2];
            e.responseTimes = [e.responseTimes, firstPress(KbName('2@')) - e.startRun];
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

e.runTime = GetSecs - e.startRun;

% analyse task accuracy
responseTimes = e.responseTimes; % make copies as we are editing these
responses = e.responses;
responseWindow = 1.5; % responses later than this count as a miss
for t = 1:size(e.targetTimes,2)
    targetTime = e.targetTimes(t);
    target = e.targets(t);
    e.hits(t) = 0; % default is a miss with RT of nan
    e.RTs(t) = nan;
    for r = 1:size(responseTimes,2)
        if ismember(responses(r),  target) && responseTimes(r) > targetTime % if response is correct and happened after target
            rt = responseTimes(r)-targetTime;
            if rt < responseWindow % and if response happened within a second of target
                e.hits(t) = 1; % mark a hit
                e.RTs(t) = rt; % store the RT
                responseTimes(r) = []; % delete this response so it can't be reused
                responses(r) = [];
                break % stop looking for responses to this target
            end
        end
    end
end
e.accuracy = (sum(e.hits)/size(e.targetTimes,2))*100;
e.meanRT = nanmean(e.RTs);

savedir = fullfile(e.root,'data',subject,session,'fillingin_rsvp_v1');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(subject , '_fillingin_rsvp_v1_sn',num2str(e.scanNum),'_rn',num2str(e.runNum),'_',e.date,'.mat'));
save(savename,'e');

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;