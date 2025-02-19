%function phantom_v3(subject, session, vertOffset, debug) 
%In this version, we add breaks
subject = 1;
session = 1;
debug = 1;
vertOffset = 0;

global EyeData rect w xc yc %eye_used
%%%% resolution
if debug == 1
    % eyetracking on (1) or off (0)
    ET = 0;
    exp.screenWidth = 17;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    exp.viewingDist = 48;             % in cm; 3Tb/office=43, miniHelm=57;
	exp.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % laptop
    exp.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
    ET = 1;
    exp.screenWidth = 16;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    exp.viewingDist = 23;             % in cm; %23 in eye tracking room 425 3Tb/office=43, miniHelm=57;
    exp.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % scanner
    exp.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
% responseKeys = zeros(1,256);
% responseKeys(KbName('1!'))=1; % button box 1
% responseKeys(KbName('2@'))=1; % button box 2

Screen('Preference', 'SkipSyncTests', 0);

exp.scanNum = input('Scan number :');
exp.runNum = input('Run number :');
exp.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m

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


%%%% 2D sine wave grating properties
exp.stim.spatialFreqDeg = 0.5; %0.286;                                          % cycles per degree of visual angle
exp.stim.contrast = 0.3 ;                                                 % in %, maybe??
exp.stim.orientation = [90]; %[90 180];                                                % in degrees
exp.stim.degFromFix = .6;                                              % in degrees of visual angle
exp.stim.gaborHDeg = 4;                                                  % in degrees of visual angle
exp.stim.gaborWDeg = 4;
%exp.stim.rectGaborWDeg = 8;
exp.stim.contrastMultiplicator = .075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
exp.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
exp.stim.cycPerSec = 1.13*3/2; % (modified to half the speed)
exp.stim.motionRate = exp.stim.cycPerSec*360;                                          % 1.13 cycles per second = 360 deg of phase *1.13 per sec
exp.stim.cycles = 3; %number of cycles shifted per lap  (modified to half the number of cycles per lap)
%%%% sine wave grating timing (within block scale)

exp.initialFixation = 6;        % in seconds
exp.finalFixation = 2;          % in seconds
exp.trialFixation = 1;          % in seconds
exp.stimDur = (exp.stim.cycles/exp.stim.cycPerSec)*2;        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
exp.stimsPerBlock = 4.75;      % number of back-and-forth laps of the stimulus drift
exp.blockLength = exp.trialFixation+ ceil(exp.stimDur*exp.stimsPerBlock);           % in seconds
exp.betweenBlocks = 2;          % in seconds
exp.flipsPerSec = 60;  % 12;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
exp.flipWin = 1/exp.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%e.numBlocks = 12;  % 6 for single and 6 for pair...

%%%% conditions & layout (across blocks scale)
%exp.conds = %{'SquMinbgSp4', 'SquMeanbgSp4'};%...
%     'SquMinbgSp3', 'SquMeanbgSp3', ...
%     'SquMinbgSp2', 'SquMeanbgSp2' ...
exp.conds = {'RectMinbgSp4', 'RectMeanbgSp4'}; %...
%     'RectMinbgSp3', 'RectMeanbgSp3',...
%     'RectMinbgSp2', 'RectMeanbgSp2'}; %'vertDoubleRect''horSingleL','horSingleR','horDoubleIndir'
exp.numConds = length(exp.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
exp.repsPerRun = 10;              % repetitions of each condition per run
exp.numBlocks = exp.numConds*exp.repsPerRun;
exp.condShuffle = Shuffle(repmat([1:exp.numConds],1,exp.repsPerRun)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order
exp.totalTime = exp.initialFixation + ((exp.numBlocks-1) * (exp.blockLength + exp.betweenBlocks)) + exp.blockLength + exp.finalFixation;
exp.allFlips = (0:exp.flipWin:exp.totalTime);

%%%% fixation 
exp.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
exp.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
exp.fontSize = 12; %26;



%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

exp.onSecs = [zeros(1,exp.initialFixation)...
    repmat([ones(1,exp.blockLength) zeros(1,exp.betweenBlocks)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
    ones(1,exp.blockLength) zeros(1,exp.finalFixation)];
exp.longFormBlocks = Expand(exp.onSecs,exp.flipsPerSec,1); %1 when block, 0 when between block
exp.longFormFlicker = repmat(ones(1,1),1,length(exp.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(exp.longFormBlocks)


%set up the timing model for stimulus pairing conditions ( pair square, pair rectangle ) and stimulus orientation
%longform condition timing, which aligns with the flicker timing
exp.longFormConds = zeros(1,exp.initialFixation);
for i = 1:exp.numBlocks-1
    exp.longFormConds = [exp.longFormConds, repmat(exp.condShuffle(i),1,exp.blockLength)]; % blocks
    exp.longFormConds = [exp.longFormConds, zeros(1,exp.betweenBlocks)]; % inter-block blanks
    
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
HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
if debug
    [w, rect]=Screen('OpenWindow',screen,exp.backgroundColor,[100 100 900 600],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
else
    [w, rect]=Screen('OpenWindow',screen,exp.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
end
Screen(w, 'TextSize', exp.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%gamma correction, file prepared for room 425
if exp.gammaCorrection 
  %load gamma correction file
  load('/Users/tongtesting2/Desktop/MonitorCal/425dell/phase2_photometry.mat')
  Screen('LoadNormalizedGammaTable', w, inverseCLUT);
end
%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate = Screen('NominalFrameRate',w);
flipTimes = [0:frameInt*frameRate/exp.flipsPerSec:exp.stimDur]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
exp.flipTimes = flipTimes;

%%% Still red dot for 1 sec before trial starts
exp.stillDotPhase = zeros(1,exp.flipsPerSec*exp.trialFixation); 

%%% Stimulus phases
% (in a linear/triangle scenario) exp.stim.dphasePerFlip = exp.stim.motionRate*frameInt * frameRate/exp.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip since we have 5 times less flips than frame refresh
exp.stim.tempPhase = pi/2; %this will make the oscillation start on the fast segment of the oscillation, with red dot on the center of the screen
exp.stim.oscillation = cos(2*pi*(1/exp.stimDur)*flipTimes+exp.stim.tempPhase);
exp.stim.spatialPhase = 90; %this makes the grating phase so that the center of the dark stripe is on the center of the screen
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
exp.rectGaborWidth = round(exp.stim.gaborWDeg*2*exp.ppd); 
%exp.driftSpeedDeg = exp.stim.cycPerSec/exp.stim.spatialFreqDeg;
%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers

%square
exp.squLWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth,length(exp.stim.orientation));
exp.squRWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth,length(exp.stim.orientation));
exp.squLWaveID = nan(length(exp.longFormFlicker),length(exp.stim.orientation));
exp.squRWaveID = nan(length(exp.longFormFlicker),length(exp.stim.orientation));

%rectangle 
%square
exp.rectLWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth*2,length(exp.stim.orientation));
exp.rectRWave = nan(length(flipTimes),exp.gaborHeight,exp.gaborWidth*2,length(exp.stim.orientation));
exp.rectLWaveID = nan(length(exp.longFormFlicker),length(exp.stim.orientation));
exp.rectRWaveID = nan(length(exp.longFormFlicker),length(exp.stim.orientation));

for o =1:length(exp.stim.orientation)
    for f = 1:length(flipTimes)
        
        LPhase = exp.stim.phases(f);
        RPhase = exp.stim.phases(f);
        %ih in pixels %iw in pixels %spatial freq in cycles per dva
        %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
        %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
        %background color (unused if the grating is not an annulus)
        % square
        exp.squLWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
            exp.stim.orientation(o),LPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
            exp.ppd);
        exp.squRWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth,exp.stim.spatialFreqDeg,...
            exp.stim.orientation(o),RPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
            exp.ppd);
        %    figure();
        %   imshow(squeeze(topWave(f,:,:)));

        tmpsquLWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.squLWave(f,:,:,o)));
        tmpsquRWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.squRWave(f,:,:,o)));
        % rect
        exp.rectLWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth*2,exp.stim.spatialFreqDeg,...
            exp.stim.orientation(o),LPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
            exp.ppd);
        exp.rectRWave(f,:,:,o) = makeSineGrating(exp.gaborHeight,exp.gaborWidth*2,exp.stim.spatialFreqDeg,...
            exp.stim.orientation(o),RPhase,exp.stim.contrastOffset(1),exp.stim.contrastMultiplicator,...
            exp.ppd);
        %    figure();
        %   imshow(squeeze(topWave(f,:,:)));

        tmprectLWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.rectLWave(f,:,:,o)));
        tmprectRWaveID(f,o) = Screen('MakeTexture', w, squeeze(exp.rectRWave(f,:,:,o)));
  
         
    end
    
    %% extend stimulus matrix to include the same total number of flips as the whole experiment
    
    %square
    exp.squLWaveID(:,o) = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
    repmat([zeros(1,length(exp.stillDotPhase)) repmat(tmpsquLWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmpsquLWaveID(:,o)))) tmpsquLWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmpsquLWaveID(:,o))),o)' zeros(1,exp.betweenBlocks*exp.flipsPerSec)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
    zeros(1,length(exp.stillDotPhase)) repmat(tmpsquLWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmpsquLWaveID(:,o)))) tmpsquLWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmpsquLWaveID(:,o))),o)' zeros(1,exp.finalFixation*exp.flipsPerSec)]';
    
    exp.squRWaveID(:,o) = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
    repmat([zeros(1,length(exp.stillDotPhase)) repmat(tmpsquRWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmpsquRWaveID(:,o)))) tmpsquRWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmpsquRWaveID(:,o))),o)' zeros(1,exp.betweenBlocks*exp.flipsPerSec)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
    zeros(1,length(exp.stillDotPhase)) repmat(tmpsquRWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmpsquRWaveID(:,o)))) tmpsquRWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmpsquRWaveID(:,o))),o)' zeros(1,exp.finalFixation*exp.flipsPerSec)]';
    
    %rect
    exp.rectLWaveID(:,o) = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
        repmat([zeros(1,length(exp.stillDotPhase)) repmat(tmprectLWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmprectLWaveID(:,o)))) tmprectLWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmprectLWaveID(:,o))),o)' zeros(1,exp.betweenBlocks*exp.flipsPerSec)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
        zeros(1,length(exp.stillDotPhase)) repmat(tmprectLWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmprectLWaveID(:,o)))) tmprectLWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmprectLWaveID(:,o))),o)' zeros(1,exp.finalFixation*exp.flipsPerSec)];
    
    exp.rectRWaveID(:,o) = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
        repmat([zeros(1,length(exp.stillDotPhase)) repmat(tmprectRWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmprectRWaveID(:,o)))) tmprectRWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmprectRWaveID(:,o))),o)' zeros(1,exp.betweenBlocks*exp.flipsPerSec)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
        zeros(1,length(exp.stillDotPhase)) repmat(tmprectRWaveID(:,o)',1,floor((exp.blockLength-exp.trialFixation)*exp.flipsPerSec/length(tmprectRWaveID(:,o)))) tmprectRWaveID(1:mod((exp.blockLength-exp.trialFixation)*exp.flipsPerSec,length(tmprectRWaveID(:,o))),o)' zeros(1,exp.finalFixation*exp.flipsPerSec)];

end

%% Sine wave gratings locations (in the task loop since it changes)
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xL = rect(3)/2; % % = stimulus center located on the horizontal center of the screen
xR = rect(3)/2; % = stimulus center located on the horizontal center of the screen

%%% create drifting red dots position
exp.driftPosDeg = exp.stim.oscillation.*exp.stim.cycles.*1/(2*exp.stim.spatialFreqDeg);
exp.fixSpatialPhase = 0; %(1/(4*exp.stim.spatialFreqDeg))*exp.ppd;
exp.driftPos = exp.driftPosDeg.*exp.ppd +exp.fixSpatialPhase;


exp.longDriftPos = [zeros(1,exp.initialFixation*exp.flipsPerSec)...
    repmat([exp.stillDotPhase exp.driftPos exp.driftPos(1:floor(length(exp.driftPos)/2)) zeros(1,(exp.blockLength-exp.trialFixation)*exp.flipsPerSec-floor(1.5*length(exp.driftPos))) zeros(1,exp.betweenBlocks*exp.flipsPerSec)],1,exp.numBlocks-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
    exp.stillDotPhase exp.driftPos exp.driftPos(1:floor(length(exp.driftPos)/2)) zeros(1,(exp.blockLength-exp.trialFixation)*exp.flipsPerSec-floor(1.5*length(exp.driftPos))) zeros(1,exp.finalFixation*exp.flipsPerSec)];


%% Eyetracking parameters

if ET 
    EyelinkSetup(0);
    eye_used = Eyelink('EyeAvailable');
    exp.ShowRealTimeGaze = [  ]; % [] or [ something ]
    exp.nGazetoShow = [ 60 ]; % current~past N fixations
end
%% %%%% initial window - wait for backtick
Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(w, 'Flip', 0);
KbTriggerWait(53, deviceNumber);

DrawFormattedText(w,'Follow the oscillating visual phantom within the gap at the center of the screen \n\n as best as you can using the red dot as a guide, even after the red dot is gone. \n\n Press Space to start'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
KbTriggerWait(KbName('Space'), deviceNumber);

%%%%%%%%%%%%%%%%%% Response listening %%%%%%%%%%%%%%%%%%%%%%%%
%KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=0;
%%%%%%% START task TASK/FLIPPING
% [e.totalTime, e.allFlips,n]

gray = repmat(mean(squeeze(exp.squLWave(1,1,:))), [1,3]);

Screen('FillRect', w, gray);
    
if ET
    gcnt = 0; 
    EyeData.mx=nan(1,ceil(exp.totalTime*1/exp.flipsPerSec));
    EyeData.my=nan(1,ceil(exp.totalTime*1/exp.flipsPerSec));
    EyeData.ma=nan(1,ceil(exp.totalTime*1/exp.flipsPerSec));
    EyeData.FixDoneT = nan(1,ceil(exp.totalTime*1/exp.flipsPerSec));
    EyeData.gazeD = nan(1,ceil(exp.totalTime*1/exp.flipsPerSec));
    EyeData.Fixated = nan(1,ceil(exp.totalTime*1/exp.flipsPerSec));
end
%EyeStart = GetSecs(); %time we start caring about eyetracking
%currPosID = i;

onOffs = [diff(exp.longFormBlocks) 0];
flipCnt = 0;
cnt = 0;
tstartcnt = 0;
while n+1 < length(exp.allFlips)

    if ET
        run GetEyeDataLoic; %check eyetracker
    end

    [exp.longFormBlocks(n+1),exp.longFormFlicker(n+1)]

    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 0
        [VBLT, exp.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        exp.flipTime(n+1) = exp.startRun;
    else
        [VBLT, exp.flipTime(n+1), FlipT, missed] = Screen(w, 'Flip', exp.startRun + exp.allFlips(n+1) - slack);
    end
    
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exp.longFormBlocks(n+1) == 1 && exp.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        thisCond = exp.longFormConds(n+1);
        %screen background color
        if strfind(conditions(thisCond).name{:}, 'Minbg')
            gray = repmat(min(min(squeeze(exp.squLWave(1,:,:)),[],1)), [1,3]);
        elseif strfind(conditions(thisCond).name{:}, 'Meanbg')
            gray = repmat(mean(squeeze(exp.squLWave(1,1,:))), [1,3]);
        end
        Screen('FillRect', w, gray);
        % draw & increment stims
        if strfind(conditions(thisCond).name{:}, 'Sp4')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            yL = rect(4)/2 - exp.stim.gaborHDeg*exp.ppd+0*exp.ppd; % stimulus located 4 degrees above screen center
            yR = rect(4)/2+ exp.stim.gaborHDeg*exp.ppd-0*exp.ppd; % stimulus located 4 degrees below screen center
            
        elseif strfind(conditions(thisCond).name{:}, 'Sp3')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            yL = rect(4)/2 - exp.stim.gaborHDeg*exp.ppd+0.5*exp.ppd; % stimulus located 4 degrees above screen center
            yR = rect(4)/2+ exp.stim.gaborHDeg*exp.ppd-0.5*exp.ppd; % stimulus located 4 degrees below screen center

        elseif strfind(conditions(thisCond).name{:}, 'Sp2')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            yL = rect(4)/2 - exp.stim.gaborHDeg*exp.ppd+1*exp.ppd; % stimulus located 4 degrees above screen center
            yR = rect(4)/2+ exp.stim.gaborHDeg*exp.ppd-1*exp.ppd; % stimulus located 4 degrees below screen center

        end
        
        if strfind(conditions(thisCond).name{:}, 'Squ')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            
            exp.LRect =  CenterRectOnPoint([0 0 exp.gaborWidth exp.gaborHeight],xL,yL);
            exp.RRect =  CenterRectOnPoint([0 0 exp.gaborWidth exp.gaborHeight],xR,yR);
            if nnz(find(exp.squLWaveID(n+1)))
                % top stim
                Screen('DrawTexture', w, exp.squLWaveID(n+1),[],exp.LRect);
                % bottom stim
                Screen('DrawTexture', w, exp.squRWaveID(n+1), [], exp.RRect);
            end
            
         elseif strfind(conditions(thisCond).name{:}, 'Rect')  %|| strcmp(conditions(thisCond).name, 'double-oppdir')  % draw second stim if it is 'double-indir' or 'double-oppdir' %if it's randomized or incognruent
            
            exp.rectLRect =  CenterRectOnPoint([0 0 exp.rectGaborWidth exp.gaborHeight],xL,yL);
            exp.rectRRect =  CenterRectOnPoint([0 0 exp.rectGaborWidth exp.gaborHeight],xR,yR);
            if nnz(find(exp.rectLWaveID(n+1)))
                % top stim
                Screen('DrawTexture', w, exp.rectLWaveID(n+1),[],exp.rectLRect);
                % bottom stim
                Screen('DrawTexture', w, exp.rectRWaveID(n+1), [], exp.rectRRect);
            end
        end
        %draw red dot
        %if nnz(exp.longDriftPos(n+1))
        flipCnt = flipCnt+1;
        %if flipCnt < 378
        if flipCnt < length([exp.stillDotPhase exp.driftPos exp.driftPos(1:floor(length(exp.driftPos)/2))])
            xOffset = exp.longDriftPos(n+1);
            Screen('FillOval', w,[255 0 0], [xc+xOffset-round(exp.fixSize/4) yc-round(exp.fixSize/4) xc+xOffset+round(exp.fixSize/4) yc+round(exp.fixSize/4)]);%black fixation solid circle
        elseif flipCnt == exp.blockLength*exp.flipsPerSec
            flipCnt = 0;
        end
    end
    if nnz(onOffs(n+1)) == 1
        thisCond = exp.longFormConds(n+2);
        cnt = cnt +1;
        time = GetSecs;
        tstartcnt = tstartcnt +1;
        if ~isempty(find(mod(tstartcnt,2))) && ET == 1
            tr = (tstartcnt-1)/2+1;
            Eyelink('Message', char(sprintf('Cond %s', conditions(thisCond).name{:})))
            Eyelink('Message', 'TRIALID %d', tr);
            Eyelink('Message', 'STIM_ONSET');
        elseif isempty(find(mod(tstartcnt,2))) && ET == 1
            Eyelink('Message', 'STIM_OFFSET');
        end
        
        
    end
    
    if nnz(cnt) && mod(cnt,4) == 0 && GetSecs-time >= 1
        DrawFormattedText(w,'Press Space whenever you feel ready'... % : press 1 as soon as letter J appears on the screen,\n\n and press 2 as soon as letter K appears on the screen. \n\n Press Space to start'...
            ,'center', 'center',[0 0 0]);
        Screen(w, 'Flip', 0);
        KbTriggerWait(KbName('Space'), deviceNumber);
        cnt = 0;
    end
    
    n = n+1;
   
end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

exp.runTime = GetSecs - exp.startRun;


savedir = fullfile(exp.root,'data',sprintf('s%d/sess%d/',subject,session),'smooth_pursuit_v6');
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/s%d_smooth_pursuit_v6_sn%d_rn%d_date%s_fix',subject,exp.scanNum,exp.runNum,num2str(exp.date)), '.mat'));
save(savename,'exp');

%KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;


if ET
    Eyelink('StopRecording');   
    Eyelink('CloseFile');
        %download edf file
%     if downET
        Eyelink('Command', 'set_idle_mode');
    try
        fprintf('Receiving data file ''%s''\n',  EyeData.edfFile);
        status = Eyelink('ReceiveFile');
        
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(EyeData.edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n',  EyeData.edfFile, pwd );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n',  EyeData.edfFile);
    end
%     end
    
    Eyelink('Shutdown');
    savename = fullfile(savedir, strcat(sprintf('/s%d_smooth_pursuit_v6_sn%d_rn%d_date%s_fix_eyeDat',subject,exp.scanNum,exp.runNum,num2str(exp.date)), '.mat'));
    save(savename, 'EyeData')
end 