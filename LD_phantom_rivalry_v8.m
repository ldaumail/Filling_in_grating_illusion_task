function LD_phantom_rivalry_v8(subject, session, debug)

%In this version, we add multiple velocities
% subject = 'Dave';                                                                                                                                                                                                                                                     
% session = 1;                                                                                                                           
% debug = 1;


ex.version = 'v8';
%%%% resolution 
if debug == 1

    ex.screenWidth = 40;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 46;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1280,1024,85); % laptop 1920,1080/ 2880, 1800 ,0
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else                                                                                                                            
    ex.screenWidth = 40;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 46;             % in cm; %23 in eye tracking                                                                                                                          room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1280,1024,85); % ET room 1600,900,60
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1'))=1; % button box 1
responseKeys(KbName('2'))=1; % button box 2
responseKeys(KbName('3'))=1; % button box 3
% responseKeys(KbName('4'))=1; % button box 4
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2
responseKeys(KbName('3#'))=1; % button box 3
% responseKeys(KbName('4$'))=1; % button box 4

Screen('Preference', 'SkipSyncTests', 0);

% ex.scanNum = input('Scan number :');
ex.runNum = input('Run number :');

%%% basic naming set-up
ex.subject = subject;
ex.session = session;


%%%% set-up rand
ex.startTime = clock;
rng(sum(100*ex.startTime));
ex.rand = rng;

%%%% files and things
ex.root = pwd;
ex.date = datestr(now,30);


%%%% 2D sine wave grating inducers properties
ex.stim.spatialFreqDeg = 0.5/2;                                           % cycles per degree of visual angle
ex.stim.contrast = 0.3 ;                                                 % in %, maybe??
ex.stim.orientation = [180]; %[90 180];                                                % in degrees
ex.stim.degFromFix = .6;                                              % in degrees of visual angle
ex.stim.gaborHDeg = 16; %6;                                                  % in degrees of visual angle
ex.stim.gaborWDeg = 6; %16;
ex.stim.distFromFixDeg = 6; %grating edge 3 deg horizontal away from fixation (grating center 6 deg away)
ex.stim.contrast = 0.15;%linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
ex.stim.contrastMultiplicator = ex.stim.contrast/2;  % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.stim.contrastOffset = [.5 0.425 0.35]; %.5 .5 0];                                  % for procedural gabor
%ex.stim.cycPerSec = [1.13*1/2,1.13*3/2]; % try multiple speeds
%ex.stim.motionRate = ex.stim.cycPerSec.*360;                                          % 1.13 cycles per second = 360 deg of phase *1.13 per sec
%ex.stim.cycles =[1, 3]; %number of cycles shifted per lap  (modified to half the number of cycles per lap)
ex.stim.cycPerSec = [1,1.4]; %drifting speed in cycles of grating per sec
ex.stim.cycles = [1 1]; %number of cycles per lap


%%%% sine wave grating timing (within block scale)
ex.stimDur = (ex.stim.cycles./ex.stim.cycPerSec)*2;        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
ex.initialFixation = 6;        % in seconds
ex.finalFixation = 2;          % in seconds
ex.trialFixation = 1;          % in seconds
%ex.stimsPerBlock = 4.5;      % number of back-and-forth laps of the stimulus drift
ex.blockLength = ex.trialFixation + 119;%119; %ex.trialFixation+ ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.betweenBlocks = 2;          % in seconds
ex.flipsPerSec = 60;  % 60;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 


%%%% Opposite eye low contrast grating (probe)
ex.lcstim.spatialFreqDeg =1;% 2;
ex.lcstim.contrast = 0.2; %0.2; %0.03;%linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
ex.lcstim.contrastMultiplicator = ex.lcstim.contrast/2;  % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.lcstim.orientation = [45 135];
ex.lcstim.gaborHDeg = 2;                                                  % in degrees of visual angle
ex.lcstim.gaborWDeg = 2;
ex.lcstim.distFromFixDeg = 2; % in degrees of visual angle, grating center 2 deg away (edge 1 deg away)

%%%% Fixation
ex.fixSizeDeg =  .2;            % in degrees, the size of the biggest white dot in the fixation
ex.bigFixSizeDeg = 0.5;
ex.outerFixPixels = 2;          % in pixels, the black ring around fixation

%%%% conditions & layout (across blocks scale)

ex.conds = {'MinbgPairLeftVel3', 'NophLeft' ...%,'MinbgTopLeftVel3', 'MinbgBotLeftVel3', 'MeanbgPairLeftVel3', 'MeanbgTopLeftVel3', 'MeanbgBotLeftVel3', 'LightbgPairLeftVel3', 'LightbgTopLeftVel3','LightbgBotLeftVel3','MinbgRightVel3', 'MeanbgRightVel3', 'LightbgRightVel3', 'NophRight'
    ... %, Here, the Left/Right indicator in the condition name corresponds to the phantom grating pair location on the screen 
    }; 
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 3;              % repetitions of each condition per run
ex.numBlocks = ex.numConds*ex.repsPerRun;
ex.lcstim.oriVec = repmat([1 2], 1, ex.numBlocks/2);
ex.condShuffle = [];
for i =1:ex.repsPerRun
    ex.condShuffle = [ex.condShuffle, Shuffle([1:ex.numConds])];
end
% for i =1:ex.repsPerRun
%     ex.condShuffle = [ex.condShuffle, Shuffle([ex.numConds/2+1:ex.numConds])];
% end

ex.totalTime = [];
for t =1:length(ex.blockLength) %there is a different block length for every drifting speed
    if t == 1
        ex.totalTime = sum([ex.totalTime, ex.initialFixation + (ex.numBlocks/length(ex.blockLength) * (ex.blockLength(t) + ex.betweenBlocks))]);
    elseif t <length(ex.blockLength) && t > 1
             ex.totalTime = sum([ex.totalTime, (ex.numBlocks/length(ex.blockLength) * (ex.blockLength(t) + ex.betweenBlocks))]); 
    elseif t == length(ex.blockLength)
        ex.totalTime = sum([ex.totalTime, ((ex.numBlocks/length(ex.blockLength)-1) * (ex.blockLength(t) + ex.betweenBlocks)) + ex.blockLength(t) + ex.finalFixation]);
    end
end
ex.allFlips = (0:ex.flipWin:ex.totalTime);
ex.allFlips = ex.allFlips(1:end-1);
ex.trialFlips = (0:ex.flipWin:ex.blockLength(1)+ex.betweenBlocks);
ex.trialFlips = ex.trialFlips(1:end-1);

%%%% screen
ex.backgroundColor = [127.4933  127.4933  127.4933];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 26;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

ex.onSecs = [ones(1,ex.blockLength(t)) zeros(1,ex.betweenBlocks)];
ex.longFormBlocks = Expand(ex.onSecs,ex.flipsPerSec,1); %1 when block, 0 when between block
length(ex.longFormBlocks)
ex.stimOnSecs = [zeros(1,ex.trialFixation) ones(1,ex.blockLength(t)-ex.trialFixation) zeros(1,ex.betweenBlocks)];
ex.longFormStimOnSecs = Expand(ex.stimOnSecs,ex.flipsPerSec,1); %1 when stim on, 0 when fixation or between blocks

% %% create the timing model of stimulus conditions for this particular run
% clear i
for i =1:ex.numConds
    conditions(i).name = ex.conds(i);
    conditions(i).startTimes = [];
end

%%
%%%%%%%%%%%%%%%
% open screen %
%%%%%%%%%%%%%%%
HideCursor;
Priority(9);
%Priority(0);
%%%% open screen
screen=max(Screen('Screens'));
if debug
    [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat);
    %[w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
    
else
    %[w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
    
    [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
end
Screen(w, 'TextSize', ex.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%gamma correction, file prepared for room 425
if ex.gammaCorrection
    % Gamma correction (run phase2_photometry.mat in 417C computer, get gamma
    % table)
    load("phase2_photometry.mat");
    Screen('LoadNormalizedGammaTable', screen, inverseCLUT);
end
%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate =  1/frameInt;%Screen('NominalFrameRate',w);
%%% Still red dot for 1 sec before trial starts
ex.stillDotPhase = zeros(1,ex.flipsPerSec*ex.trialFixation);

flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.blockLength(1)-ex.trialFixation]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
ex.stim.flipTimes = flipTimes(1:length(flipTimes)-1);
nconds = sum(contains(ex.conds, 'bg'));
ex.stim.tempPhase1 = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun);
ex.stim.tempPhase2 = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun);
ex.stim.oscillation1 = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun,length(ex.stim.flipTimes));
ex.stim.oscillation2 = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun,length(ex.stim.flipTimes));
ex.stim.phases = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun,length(ex.stim.flipTimes));

clear c r l
for c =1:nconds
    for l =1:length(ex.stim.contrastOffset)
        for r = 1:ex.repsPerRun
            ex.stim.tempPhase1(c,l,r) = rand(1,1)*2*pi;
            ex.stim.tempPhase2(c,l,r) = rand(1,1)*2*pi;
            ex.stim.oscillation1(c,l,r,:) = cos(2*pi*(1/ex.stimDur(1))*ex.stim.flipTimes+ex.stim.tempPhase1(c,l,r));
            ex.stim.oscillation2(c,l,r,:) = cos(2*pi*(1/ex.stimDur(2))*ex.stim.flipTimes+ex.stim.tempPhase2(c,l,r));
            ex.stim.spatialPhase = 90;
            ex.stim.phases(c,l,r,:) = (ex.stim.oscillation1(c,l,r,:).*180*ex.stim.cycles(1)+ ex.stim.oscillation2(c,l,r,:).*180*ex.stim.cycles(2))/2 + ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory
        end
    end
end
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
ex.bigFixSize = round(ex.bigFixSizeDeg*ex.ppd);
ex.gaborHeight = round(ex.stim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.gaborWidth = round(ex.stim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects
ex.rawGaborHeight = ex.gaborHeight*3;
ex.rawGaborWidth = ex.gaborWidth*1.5;
ex.stim.distFromFix = round(ex.stim.distFromFixDeg*ex.ppd);

%%%% scale the prob params for the screen
ex.probeHeight = round(ex.lcstim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.probeWidth = round(ex.lcstim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects
ex.rawProbeHeight = ex.probeHeight*1;
ex.rawProbeWidth = ex.probeWidth*1;
ex.lcstim.distFromFix = round(ex.lcstim.distFromFixDeg*ex.ppd);

%% Create only one big sinewave grating image saved for each repetition and each condition

ex.rectSWave = nan(nconds, length(ex.stim.contrastOffset), ex.repsPerRun,ex.rawGaborHeight,ex.rawGaborWidth);
ex.rectSWaveID = nan(nconds, length(ex.stim.contrastOffset), ex.repsPerRun);
clear c r
for c =1:nconds
    for l = 1:length(ex.stim.contrastOffset)
        for r = 1:ex.repsPerRun %only save the first image of each trial, that we will move during the trial
            phase = ex.stim.phases(c,l,r,1);
            ex.rectSWave(c,l,r,:,:) = makeSineGrating(ex.rawGaborHeight,ex.rawGaborWidth,ex.stim.spatialFreqDeg,...
                ex.stim.orientation,phase,ex.stim.contrastOffset(l),ex.stim.contrastMultiplicator,...
                ex.ppd);
            ex.rectSWaveID(c,l,r) = Screen('MakeTexture', w, squeeze(ex.rectSWave(c,l,r,:,:)));
        end
    end
end
%check luminances ranges
% repmat(min(min(squeeze(ex.rectSWave(1,3,1,:,:)),[],1)),1)
% repmat(max(max(squeeze(ex.rectSWave(1,3,1,:,:)),[],1)),1)

%% create 1 low contrast grating image as a probe for other eye
gray2 = [108.3750  108.3750  108.3750]; %repmat(min(min(squeeze(ex.rectSWave(1,1,:,:)),[],1)), [1,3]);

phase = ex.stim.phases(c,l,r,1);
ex.lcSWave = nan(length(ex.lcstim.orientation),ex.rawProbeHeight,ex.rawProbeWidth);
for i =1:length(ex.lcstim.orientation)
    ex.lcSWave(i,:,:) = makeSineGrating(ex.rawProbeHeight,ex.rawProbeWidth,ex.lcstim.spatialFreqDeg,...
        ex.lcstim.orientation(i),phase,ex.stim.contrastOffset(2),ex.lcstim.contrastMultiplicator,...
        ex.ppd, 0, round(1*ex.ppd), gray2(1));
    ex.lcSWaveID(i) = Screen('MakeTexture', w, squeeze(ex.lcSWave(i,:,:)));
end
%check luminances ranges
% repmat(min(min(squeeze(ex.lcSWave(1,:,:)),[],1)),1)
% repmat(max(max(squeeze(ex.lcSWave(1,:,:)),[],1)),1)

% figure();
% imshow(squeeze(ex.lcSWave(1,:,:))/255)
% figure();
% imshow(squeeze(ex.lcSWave(2,:,:))/255)


%% create drifting red dots position
clear flipTimes
flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.blockLength(1)+ex.betweenBlocks]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
flipTimes = flipTimes(1:length(flipTimes)-1);
ex.stimDriftPosDeg = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun,length(ex.stim.oscillation1(1,1,1,:)));
ex.stimDriftPos = nan(nconds,length(ex.stim.contrastOffset),ex.repsPerRun,length(ex.stim.oscillation1(1,1,1,:)));
ex.stimLongDriftPos = nan(nconds,length(ex.stim.contrastOffset), ex.repsPerRun,length(flipTimes));
clear c r
for c = 1:nconds
    for l =1:length(ex.stim.contrastOffset)
        for r = 1:ex.repsPerRun
            ex.stimDriftPosDeg(c,l,r,:) = (ex.stim.oscillation1(c,l,r,:).*ex.stim.cycles(1).*1/(2*ex.stim.spatialFreqDeg)); %+ ex.stim.oscillation2(c,l,r,:).*ex.stim.cycles(2).*1/(2*ex.stim.spatialFreqDeg))/2;
            ex.stimFixSpatialPhase = 0; %-(1/(8*ex.stim.spatialFreqDeg))*ex.ppd;%(1/(4*ex.stim.spatialFreqDeg))*ex.ppd;
            ex.stimDriftPos(c,l,r,:) = ex.stimDriftPosDeg(c,l,r,:).*ex.ppd +ex.stimFixSpatialPhase;
            ex.stimStillDotPhase = ex.stimDriftPos(c,l,r,1);
            ex.stimLongDriftPos(c,l,r,:) = [repmat(ex.stimStillDotPhase,1,ex.flipsPerSec*ex.trialFixation) squeeze(ex.stimDriftPos(c,l,r,:))' ...
                zeros(1,ex.betweenBlocks*ex.flipsPerSec)];
        end
    end
end

% figure();
% subplot(4,1,1)
% plot(squeeze(ex.stimDriftPosDeg(1,1,:)))
% subplot(4,1,2)
% plot(squeeze(ex.stimDriftPosDeg(1,2,:)))
% subplot(4,1,3)
% plot(squeeze(ex.stimDriftPosDeg(1,3,:)))
% subplot(4,1,4)
% plot(squeeze(ex.stimDriftPosDeg(1,4,:)))
% ylabel('Stimulus position (dva)')
% xlabel('Flip #')

%% Sine wave gratings locations 
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

%left grating
xLl = (1/2)*xc-(1/2)*ex.gaborWidth-ex.stim.distFromFix;
xLr = (1/2)*xc+(1/2)*ex.gaborWidth-ex.stim.distFromFix;

%right grating
xRr = (1/2)*xc-(1/2)*ex.gaborWidth+ex.stim.distFromFix;
xRl = (1/2)*xc+(1/2)*ex.gaborWidth+ex.stim.distFromFix;

% grating y locations
yt = yc-(1/2)*ex.gaborHeight;
yb = yc+(1/2)*ex.gaborHeight;

 
%% Create rectangular masks for the gratings
% gray1 = [127.4933  127.4933  127.4933]; %repmat(mean(squeeze(ex.rectSWave(1,1,1,:))), [1,3]);

%phantom control condition
coLPaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',coLPaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',coLPaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',coLPaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

coLTaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',coLTaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',coLTaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

coLBaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',coLBaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',coLBaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window


coRaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',coRaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',coRaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',coRaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth+xc)-ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth-xc)/2-ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window


%phantom condition 1
ph1LPaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph1LPaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',ph1LPaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',ph1LPaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

ph1LTaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph1LTaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',ph1LTaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

ph1LBaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph1LBaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',ph1LBaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

ph1Raperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph1Raperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',ph1Raperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',ph1Raperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth + xc)-ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2-ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

ph2LPaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph2LPaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',ph2LPaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',ph2LPaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

ph2LTaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph2LTaperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',ph2LTaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

ph2LBaperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph2LBaperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',ph2LBaperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth - xc)+ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2+ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window


ph2Raperture=Screen('OpenOffscreenwindow', w, gray2);
Screen('FillRect',ph2Raperture, [255 255 255 0], [xLl yt xLr yb]); %Left grating window
Screen('FillRect',ph2Raperture, [255 255 255 0], [xRr yt xRl yb]); %Right grating window
Screen('FillRect',ph2Raperture, [255 255 255 0], [xc-(1/2)*(ex.probeWidth + xc)-ex.lcstim.distFromFix yc-ex.probeHeight/2 xc+(ex.probeWidth+ xc)/2-ex.lcstim.distFromFix yc+(1/2)*ex.probeHeight]); %opposite eye grating window

%% %%%% initial window - wait for backtick
DrawFormattedText(w,'Fixate the fixation dot as best as you can. \n\n After each perceptual change, report using the following digits \n\n 1: The central grating (probe) is fully visible \n\n 2: The probe is partially visible \n\n 3: The probe is fully faded/invisible \n\n Press Space to start'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
%WaitSecs(2);
KbTriggerWait(KbName('Space'), deviceNumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% START task TASK/FLIPPING
clear l
n = 1;
blockCnt = 1;
cnt = 0; %stim onset/ stime offset count

cntConds = zeros(nconds,1); %count number of summed gratings condition trials, as the phase is random, we need to have a different phase at each new trial

ex.responses = [];
ex.responseTimes=[];
ex.correctResp = [];


onOffs = [diff(ex.longFormBlocks) 0];
bLength = ex.blockLength(1);
ex.flipTime = nan(length(ex.trialFlips),length(ex.condShuffle));

KbQueueCreate(deviceNumber,responseKeys);
%%% initial fixation
if n == 1 && blockCnt == 1 %for first block
    ex.tasktstart = clock;
    ex.startRun = GetSecs();
    Screen('FillRect', w, ex.backgroundColor);
    Screen('DrawDots', w, [xc/2 yc], ex.fixSize, [255 255 255], [], 2);
    Screen('DrawDots', w, [xc*3/2 yc], ex.fixSize, [255 255 255], [], 2);
    Screen(w, 'Flip', 0);
    WaitSecs(ex.initialFixation);
end
%%% Launch the task
for c = 1:length(ex.condShuffle)
    cnt = cnt+1;
    condNum = ex.condShuffle(c);
    condName = conditions(condNum).name{:};
    oriNum = ex.lcstim.oriVec(c);
    %%for each condition, we specify the parameters values before we flip
    %%over the gratings phases
    %screen background color
    if condNum <= nconds
     cntConds(condNum) = cntConds(condNum)+1;
    end
    if contains(condName, 'Minbg') %contains(condName, 'Minbg')
        l =1;

    elseif contains(condName, 'Meanbg') %contains(condName, 'Minbg')
        l = 2;

    elseif contains(condName, 'Lightbg') %contains(condName, 'Minbg')
        l = 3;

    end
    
    %flip through the block and following between block time
    while n <= length(ex.trialFlips)
        KbQueueStart(); % response time
        ex.longFormBlocks(n)
        %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if contains(condName, 'Minbg')
            Screen('FillRect', w, gray2);
            if nnz(find(ex.longFormStimOnSecs(n)))
                if contains(condName, 'Left')
                    if contains(condName, 'Pair')
                        xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        % stim
                        
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture',w,ph1LPaperture);

                    elseif contains(condName, 'Top')
                        xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        % stim
                        
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture',w,ph1LTaperture);

                    elseif contains(condName, 'Bot')
                        xOffset = ex.stimLongDriftPos(condNum,l,cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l,cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        % stim
                        
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture',w,ph1LBaperture);
                    end

                elseif contains(condName, 'Right')
                    xOffset = ex.stimLongDriftPos(condNum,l,cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l,cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                    ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc+xc/2+xOffset,yc);
                    ex.lcRRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc-xc/2-ex.lcstim.distFromFix,yc);
                    % stim
                    Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                    Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcRRect);
                    Screen('DrawTexture',w,ph1Raperture);
                    
                end
                Screen('DrawDots', w, [xc/2 yc], ex.fixSize, [255 255 255], [], 2);
                Screen('DrawDots', w, [xc*3/2 yc], ex.fixSize, [255 255 255], [], 2);
            end
            
        elseif contains(condName, 'Meanbg')
            Screen('FillRect', w, gray2);
            if nnz(find(ex.longFormStimOnSecs(n)))
                
                if contains(condName, 'Left')
                    if contains(condName, 'Pair')
                        xOffset = ex.stimLongDriftPos(condNum,l,cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l,cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        
                        % stim
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture',w,coLPaperture);

                     elseif contains(condName, 'Top')
                         xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                         ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                         ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                         
                         % stim
                         Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                         Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                         Screen('DrawTexture',w,coLTaperture);

                     elseif contains(condName, 'Bot')
                         xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                         ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                         ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                         
                         % stim
                         Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                         Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                         Screen('DrawTexture',w,coLBaperture);

                    end
                         
                elseif contains(condName, 'Right')
                    xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                    ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc+xc/2+xOffset,yc);
                    ex.lcRRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc-xc/2-ex.lcstim.distFromFix,yc);
                    
                    % stim
                    Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcRRect);
                    Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                    Screen('DrawTexture',w,coRaperture);
                    
                end
                Screen('DrawDots', w, [xc/2 yc], ex.fixSize, [255 255 255], [], 2);
                Screen('DrawDots', w, [xc*3/2 yc], ex.fixSize, [255 255 255], [], 2);
            end
        elseif contains(condName, 'Lightbg')
            Screen('FillRect', w, gray2);
            if nnz(find(ex.longFormStimOnSecs(n)))
                
                if contains(condName, 'Left')
                    if contains(condName,'Pair')
                        xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        
                        % stim
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture',w,ph2LPaperture);

                    elseif contains(condName,'Top')
                        xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        
                        % stim
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture',w,ph2LTaperture);
                    elseif contains(condName,'Bot')
                        xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc-xc/2+xOffset,yc);
                        ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                        
                        % stim
                        Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                        Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                        Screen('DrawTexture',w,ph2LBaperture);
                    end
                elseif contains(condName, 'Right')
                    xOffset = ex.stimLongDriftPos(condNum,l, cntConds(condNum),n)-ex.stimLongDriftPos(condNum,l, cntConds(condNum),1); %baseline correct the position since every image already hase a spatial phase shift in the sinewave
                    ex.rectLRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xc+xc/2+xOffset,yc);
                    ex.lcRRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc-xc/2-ex.lcstim.distFromFix,yc);
                    
                    % stim
                    Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcRRect);
                    Screen('DrawTexture', w, ex.rectSWaveID(condNum,l, cntConds(condNum)),[],ex.rectLRect);
                    Screen('DrawTexture',w,ph2Raperture);
                end
                Screen('DrawDots', w, [xc/2 yc], ex.fixSize, [255 255 255], [], 2);
                Screen('DrawDots', w, [xc*3/2 yc], ex.fixSize, [255 255 255], [], 2);
            end
        elseif contains(condName, 'Noph') %
            Screen('FillRect', w, gray2);
            if nnz(find(ex.longFormStimOnSecs(n)))
                
                if contains(condName, 'Left')
                    ex.lcLRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc+xc/2+ex.lcstim.distFromFix,yc);
                    
                    % stim
                    Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcLRect);
                    Screen('DrawTexture',w,coLPaperture);

                elseif contains(condName, 'Right')
                    ex.lcRRect =  CenterRectOnPoint([0 0 ex.rawProbeWidth ex.rawProbeHeight],xc-xc/2-ex.lcstim.distFromFix,yc);
                    % stim
                    Screen('DrawTexture', w, ex.lcSWaveID(oriNum),[],ex.lcRRect);
                    Screen('DrawTexture',w,coRaperture);
                end
                Screen('DrawDots', w, [xc/2 yc], ex.fixSize, [255 255 255], [], 2);
                Screen('DrawDots', w, [xc*3/2 yc], ex.fixSize, [255 255 255], [], 2);
            end
        end
        
        %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 1
            [VBLT, ex.startTrial, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
            flipTimes = ex.startTrial;
            
        else
            [VBLT,flipTime, FlipT, missed] = Screen(w, 'Flip',ex.startTrial + ex.trialFlips(n) - slack); %,   %%% ex.flipTime(n,c)
            flipTimes = [flipTimes, flipTime];
            
        end
        
        %             DrawFormattedText(w,'Press 1: the phantom is visible \n\n 2: the low contrast grating and phantom are both visible \n\n 3: phantom totally disappears'... % : press 1 as soon as letter J appears on the screen,\n\n and press 2 as soon as letter K appears on the screen. \n\n Press Space to start'...
        %                 ,'center', 'center',[0 0 0]);
        
        KbQueueStop();
        [pressed, firstPress]= KbQueueCheck();
        %KbTriggerWait(KbName('1!'), deviceNumber);
        if  (pressed == 1) && ((firstPress(KbName('1!')) > 0 ||(firstPress(KbName('1')) > 0)) || (firstPress(KbName('2@')) > 0 || firstPress(KbName('2')) > 0) || (firstPress(KbName('3#')) > 0 || firstPress(KbName('3')) > 0)) %|| (firstPress(KbName('4$')) > 0 || firstPress(KbName('4')) > 0)
            ex.responses = [ex.responses, 1];
            if (firstPress(KbName('1')) > 0)
                ex.correctResp = [ex.correctResp, 1];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('1')) - ex.startRun];
            elseif (firstPress(KbName('2')) > 0)
                ex.correctResp = [ex.correctResp, 2];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('2')) - ex.startRun];
            elseif (firstPress(KbName('3')) > 0)
                ex.correctResp = [ex.correctResp, 3];
                ex.responseTimes = [ex.responseTimes, firstPress(KbName('3')) - ex.startRun];
%             elseif (firstPress(KbName('4')) > 0)
%                 ex.correctResp = [ex.correctResp, 3];
%                 ex.responseTimes = [ex.responseTimes, firstPress(KbName('4')) - ex.startRun];
            end
            pressed = 0;
        end

        if nnz(onOffs(n)) == 1
            time = GetSecs;
            cnt = cnt+1;
        end
        if (cnt/2 == 1 && GetSecs-time >= 1) && c ~= length(ex.condShuffle)
            
            %             [ex.respT(cnt),~,~] =KbWait(deviceNumber,2);
            DrawFormattedText(w,'Press Space whenever you feel ready'... % : press 1 as soon as letter J appears on the screen,\n\n and press 2 as soon as letter K appears on the screen. \n\n Press Space to start'...
                ,'center', 'center',[0 0 0]);
            Screen(w, 'Flip', 0);
            [~,~,~] =KbWait(deviceNumber,2);
            %             KbTriggerWait(KbName('Space'), deviceNumber);
            cnt = 0;
        end
        %%%% refresh queue for next character
        KbQueueFlush();
        n = n+1;
    end
    
    ex.flipTime(:,c) = flipTimes;
    n = 1;
end

%     if n == 1382%360%421%420%359%1382%1561%
%         break;
%     end

%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

ex.runTime = GetSecs - ex.startRun;

savedir = fullfile(ex.root,'data',sprintf('%s_%s/',subject,ex.version));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/%s_binocular_rivalry_%s_date%s_fix',subject,ex.version,num2str(ex.date)), '.mat'));
%save(savename,'ex');
save(savename,'ex','-v7.3')

ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;                  