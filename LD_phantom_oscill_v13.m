%In this version, we add multiple velocities
subject = 'Lasya'                                                                                                                              ';
session = 1;
debug = 0;
vertOffset = 0;

ex.version = 'v13';
global EyeData rect w xc yc %eye_used
%%%% resolution 
if debug == 1
    % eyetracking on (1) or off (0)
    ET = 0;
    ex.screenWidth = 53.1;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1600,900,60); % laptop 1920,1080
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
    ET = 1;
    ex.screenWidth = 53.1;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; %23 in eye tracking room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1600,900,60); % scanner
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
% responseKeys = zeros(1,256);
% responseKeys(KbName('1!'))=1; % button box 1
% responseKeys(KbName('2@'))=1; % button box 2

Screen('Preference', 'SkipSyncTests', 0);

% ex.scanNum = input('Scan number :');
% ex.runNum = input('Run number :');
ex.vertOffset = vertOffset;    % vertical offset from FindScreenSize.m


%%% basic naming set-up
ex.subject = subject;
ex.session = session;


%%%% set-up rand
 rand('twister', sum(100*clock));
 ex.rand = rand;

% rng(sum(100*clock));
% ex.rand = rng;
%%%% files and things
ex.root = pwd;
ex.date = datestr(now,30);


%%%% 2D sine wave grating properties
ex.stim.spatialFreqDeg = 0.5/2; %0.286;                                          % cycles per degree of visual angle
ex.stim.contrast = 0.3 ;                                                 % in %, maybe??
ex.stim.orientation = [90]; %[90 180];                                                % in degrees
ex.stim.degFromFix = .6;                                              % in degrees of visual angle
ex.stim.gaborHDeg = 8;                                                  % in degrees of visual angle
ex.stim.gaborWDeg = 8;
%ex.stim.rectGaborWDeg = 8;
ex.stim.contrastMultiplicator = .075;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
ex.stim.cycPerSec = [1.13*1/2,1.13*3/2]; % try multiple speeds
ex.stim.motionRate = ex.stim.cycPerSec.*360;                                          % 1.13 cycles per second = 360 deg of phase *1.13 per sec
ex.stim.cycles =[1, 3]; %number of cycles shifted per lap  (modified to half the number of cycles per lap)
%%%% sine wave grating timing (within block scale)

ex.initialFixation = 6;        % in seconds
ex.finalFixation = 2;          % in seconds
ex.trialFixation = 1;          % in seconds
ex.stimDur = (ex.stim.cycles./ex.stim.cycPerSec)*2;        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
ex.stimsPerBlock = 4.5;      % number of back-and-forth laps of the stimulus drift
ex.blockLength = ex.trialFixation+ ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.betweenBlocks = 2;          % in seconds
ex.flipsPerSec = 60;  % 12;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%e.numBlocks = 12;  % 6 for single and 6 for pair...

%%%% conditions & layout (across blocks scale)

ex.conds = {'MeanbgRedDotVel1',...
    'RectMinbgSp8Vel1', 'RectMeanbgSp8Vel1',...   
    'MeanbgRedDotVel2','RectMinbgSp8Vel2', 'RectMeanbgSp8Vel2',...
    }; %...
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 20;              % repetitions of each condition per run
ex.numBlocks = ex.numConds*ex.repsPerRun;
condShuffle1 = Shuffle(repmat([1:ex.numConds/2],1,ex.repsPerRun)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order
condShuffle2 = Shuffle(repmat([ex.numConds/2+1:ex.numConds],1,ex.repsPerRun));
%condShuffle3 = Shuffle(repmat([ex.numConds*2/(length(ex.conds)/2)+1:ex.numConds],1,ex.repsPerRun));

ex.condShuffle = [condShuffle1 condShuffle2];

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

%%%% fixation 
ex.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
ex.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 12; %26;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%
   clear t
for t =1:length(ex.blockLength) %there is a different block length for every drifting speed
    if t == 1
        ex.onSecs = [zeros(1,ex.initialFixation)...
            repmat([ones(1,ex.blockLength(t)) zeros(1,ex.betweenBlocks)],1,ex.numBlocks/length(ex.blockLength))];
    elseif t <length(ex.blockLength) && t > 1
        ex.onSecs = [ex.onSecs...
            repmat([ones(1,ex.blockLength(t)) zeros(1,ex.betweenBlocks)],1,ex.numBlocks/length(ex.blockLength))];

    elseif t == length(ex.blockLength)
                ex.onSecs = [ex.onSecs...
                    repmat([ones(1,ex.blockLength(t)) zeros(1,ex.betweenBlocks)],1,ex.numBlocks/length(ex.blockLength)-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
                    ones(1,ex.blockLength(t)) zeros(1,ex.finalFixation)];
    end
end
ex.longFormBlocks = Expand(ex.onSecs,ex.flipsPerSec,1); %1 when block, 0 when between block
ex.longFormFlicker = repmat(ones(1,1),1,length(ex.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(ex.longFormBlocks)


%set up the timing model for stimulus pairing conditions ( pair square, pair rectangle ) and stimulus orientation
%longform condition timing, which aligns with the flicker timing
clear t
ex.longFormConds = zeros(1,ex.initialFixation);
for t =1:length(ex.blockLength) %length(blockLength) indicates the number of different velocities
    if t <length(ex.blockLength)
        for i = 1+(ex.numBlocks/length(ex.blockLength))*(t-1):ex.numBlocks/length(ex.blockLength)+(ex.numBlocks/length(ex.blockLength))*(t-1)
            ex.longFormConds = [ex.longFormConds, repmat(ex.condShuffle(i),1,ex.blockLength(t))]; % blocks
            ex.longFormConds = [ex.longFormConds, zeros(1,ex.betweenBlocks)]; % inter-block blanks
        end
    elseif t == length(ex.blockLength)
%         take care of the last speed trials
        for i =  1+(ex.numBlocks/length(ex.blockLength))*(t-1):ex.numBlocks/length(ex.blockLength)+(ex.numBlocks/length(ex.blockLength))*(t-1)-1
            ex.longFormConds = [ex.longFormConds, repmat(ex.condShuffle(i),1,ex.blockLength(end))]; % blocks
            ex.longFormConds = [ex.longFormConds, zeros(1,ex.betweenBlocks)]; % inter-block blanks
        end
    end
end
ex.longFormConds = [ex.longFormConds, repmat(ex.condShuffle(end),1,ex.blockLength(end)), zeros(1,ex.finalFixation)]; % the last block
ex.longFormConds = Expand(ex.longFormConds, ex.flipsPerSec,1);
length(ex.longFormConds)
% %% create the timing model of stimulus conditions for this particular run
clear i
for i =1:ex.numConds
    conditions(i).name = ex.conds(i);
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
%Priority(0);
%%%% open screen
screen=max(Screen('Screens'));
if debug
     [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
%      [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425

else
    [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
end
Screen(w, 'TextSize', ex.fontSize);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%gamma correction, file prepared for room 425
if ex.gammaCorrection 
  %load gamma correction file
  load('/Users/tonglab/Desktop/monitor_calibration/425dell_22-12-09/phase2_photometry_22-12-09.mat')
  Screen('LoadNormalizedGammaTable', w, inverseCLUT);
end
%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate =  1/frameInt;%Screen('NominalFrameRate',w);
clear t
for t =1:length(ex.stimDur)
    flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.stimDur(t)]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
    flipTimes = flipTimes(1:length(flipTimes)-1);
    flipt =sprintf('flipTimes%d',t);
    ex.(flipt) = flipTimes;
    
    %%% Still red dot for 1 sec before trial starts
    ex.stillDotPhase = zeros(1,ex.flipsPerSec*ex.trialFixation);
    
    %%% Stimulus phases
    % (in a linear/triangle scenario) ex.stim.dphasePerFlip = ex.stim.motionRate*frameInt * frameRate/ex.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip since we have 5 times less flips than frame refresh
    ex.stim.tempPhase =pi; % pi/2; %this will make the oscillation start on the fast segment of the oscillation, with red dot on the center of the screen
    
    oscillation = sprintf('oscillation%d',t);
    ex.stim.(oscillation) = cos(2*pi*(1/ex.stimDur(t))*flipTimes+ex.stim.tempPhase);
    ex.stim.spatialPhase = 90; %this makes the grating phase so that the center of the dark stripe is on the center of the screen
    phases = sprintf('phases%d',t);
    ex.stim.(phases) = ex.stim.(oscillation).*180*ex.stim.cycles(t)+ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory
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
ex.gaborHeight = round(ex.stim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.gaborWidth = round(ex.stim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects
ex.rectGaborWidth = round(ex.stim.gaborWDeg*2*ex.ppd); 
%ex.driftSpeedDeg = ex.stim.cycPerSec/ex.stim.spatialFreqDeg;
%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers
%rect
rectLWaveIDO = [zeros(ex.initialFixation*ex.flipsPerSec,1)];
rectRWaveIDO = [zeros(ex.initialFixation*ex.flipsPerSec,1)];
clear t
for o =1:length(ex.stim.orientation)
    for t =1:length(ex.stimDur)
        flipt =sprintf('flipTimes%d',t);
        flipTimes = ex.(flipt);
        
        %rect
        rectLWave = sprintf('rectLWave%d',t);
        rectRWave = sprintf('rectRWave%d',t);
        
        ex.(rectLWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2,length(ex.stim.orientation));
        ex.(rectRWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2,length(ex.stim.orientation));
        %     ex.rectLWaveID = nan(length(ex.longFormFlicker),length(ex.stim.orientation));
        %     ex.rectRWaveID = nan(length(ex.longFormFlicker),length(ex.stim.orientation));
        
        phaseNum = sprintf('phases%d',t);
        phases = ex.stim.(phaseNum) ;
 %       clear tmprectLWaveID tmprectRWaveID
        for f = 1:length(flipTimes)

            LPhase = phases(f);
            RPhase = phases(f);
            %ih in pixels %iw in pixels %spatial freq in cycles per dva
            %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
            %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
            %background color (unused if the grating is not an annulus)
            
            % rect
            ex.(rectLWave)(f,:,:,o) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                ex.stim.orientation(o),LPhase,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
                ex.ppd);
            ex.(rectRWave)(f,:,:,o) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                ex.stim.orientation(o),RPhase,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
                ex.ppd);
%                 figure();
%                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
%             
            tmprectLWaveID(f,o) = Screen('MakeTexture', w, squeeze(ex.(rectLWave)(f,:,:,o)));
            tmprectRWaveID(f,o) = Screen('MakeTexture', w, squeeze(ex.(rectRWave)(f,:,:,o)));
            
        end
        
        %% extend stimulus matrix to include the same total number of flips as the whole experiment
        
        if t<length(ex.stimDur)
            %rect
            rectLWaveIDO = [rectLWaveIDO;...
                repmat([zeros(length(ex.stillDotPhase),1); repmat(tmprectLWaveID(:,o),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectLWaveID(:,o))),1); tmprectLWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectLWaveID(:,o))),o); zeros(ex.betweenBlocks*ex.flipsPerSec,1)],ex.numBlocks/length(ex.blockLength),1)];
            
            rectRWaveIDO = [rectRWaveIDO;...
                repmat([zeros(length(ex.stillDotPhase),1); repmat(tmprectRWaveID(:,o),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectRWaveID(:,o))),1); tmprectRWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectRWaveID(:,o))),o); zeros(ex.betweenBlocks*ex.flipsPerSec,1)],ex.numBlocks/length(ex.blockLength),1)];
            
        elseif t == length(ex.stimDur)
            
            %rect
            rectLWaveIDO = [rectLWaveIDO;...
                repmat([zeros(length(ex.stillDotPhase),1); repmat(tmprectLWaveID(:,o),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectLWaveID(:,o))),1); tmprectLWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectLWaveID(:,o))),o); zeros(ex.betweenBlocks*ex.flipsPerSec,1)],ex.numBlocks/length(ex.blockLength)-1,1);... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
                zeros(length(ex.stillDotPhase),1); repmat(tmprectLWaveID(:,o),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectLWaveID(:,o))),1); tmprectLWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectLWaveID(:,o))),o); zeros(ex.finalFixation*ex.flipsPerSec,1)];
            ex.rectLWaveID(1:length(rectLWaveIDO),o) = rectLWaveIDO;
            
            rectRWaveIDO = [rectRWaveIDO;...
                repmat([zeros(length(ex.stillDotPhase),1); repmat(tmprectRWaveID(:,o),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectRWaveID(:,o))),1); tmprectRWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectRWaveID(:,o))),o); zeros(ex.betweenBlocks*ex.flipsPerSec,1)],ex.numBlocks/length(ex.blockLength)-1,1);... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
                zeros(length(ex.stillDotPhase),1); repmat(tmprectRWaveID(:,o),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectRWaveID(:,o))),1); tmprectRWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectRWaveID(:,o))),o); zeros(ex.finalFixation*ex.flipsPerSec,1)];
            ex.rectRWaveID(1:length(rectRWaveIDO),o) =rectRWaveIDO ;
        end
    end
end
%% Sine wave gratings locations (in the task loop since it changes)
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xL = rect(3)/2; % % = stimulus center located on the horizontal center of the screen
xR = rect(3)/2; % = stimulus center located on the horizontal center of the screen

yL = rect(4)/2 - ex.stim.gaborHDeg*ex.ppd+0*ex.ppd; % stimulus located 4 degrees above screen center
yR = rect(4)/2+ ex.stim.gaborHDeg*ex.ppd-0*ex.ppd; % stimulus located 4 degrees below screen cente

%% create drifting red dots position
for t =1:length(ex.stimDur)
    oscillation = sprintf('oscillation%d',t);
    driftPosDeg = sprintf('driftPosDeg%d',t);
    ex.(driftPosDeg) = ex.stim.(oscillation).*ex.stim.cycles(t).*1/(2*ex.stim.spatialFreqDeg);
    ex.fixSpatialPhase = 0; %(1/(4*ex.stim.spatialFreqDeg))*ex.ppd;
    driftPos = sprintf('driftPos%d',t);
    ex.(driftPos) = ex.(driftPosDeg).*ex.ppd +ex.fixSpatialPhase;
    
    stillDotPhase = sprintf('stillDotPhase%d',t);
    ex.(stillDotPhase)= ex.(driftPos)(1);
    
    if t == 1
        ex.longDriftPos = [zeros(1,ex.initialFixation*ex.flipsPerSec)...
            repmat([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:floor(length(ex.(driftPos))/4)) ...
            zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(1.25*length(ex.(driftPos))))...
            zeros(1,ex.betweenBlocks*ex.flipsPerSec)],1,ex.numBlocks/length(ex.blockLength))];
      
        % make lonDriftPos for the red dot only condition
        ex.longDriftPos2 = [zeros(1,ex.initialFixation*ex.flipsPerSec)...
            repmat([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) repmat(ex.(driftPos),1,4) ... 
            ex.(driftPos)(1:floor(length(ex.(driftPos))/2)) zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(4.5*length(ex.(driftPos))))...
            zeros(1,ex.betweenBlocks*ex.flipsPerSec)],1,ex.numBlocks/length(ex.blockLength))];
        
    elseif t <length(ex.blockLength) && t > 1
        ex.longDriftPos = [ex.longDriftPos...
            repmat([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:floor(length(ex.(driftPos))/4))...
            zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(1.25*length(ex.(driftPos)))) ...
            zeros(1,ex.betweenBlocks*ex.flipsPerSec)],1,ex.numBlocks/length(ex.blockLength))];
        
         % make lonDriftPos for the red dot only condition
        ex.longDriftPos2 = [ex.longDriftPos2...
            repmat([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) repmat(ex.(driftPos),1,4) ... 
            ex.(driftPos)(1:floor(length(ex.(driftPos))/2)) zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(4.5*length(ex.(driftPos))))...
            zeros(1,ex.betweenBlocks*ex.flipsPerSec)],1,ex.numBlocks/length(ex.blockLength))];
        
    elseif t == length(ex.blockLength)
        
        ex.longDriftPos = [ex.longDriftPos...
            repmat([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:floor(length(ex.(driftPos))/4))...
            zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(1.25*length(ex.(driftPos))))...
            zeros(1,ex.betweenBlocks*ex.flipsPerSec)],1,ex.numBlocks/length(ex.blockLength)-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
            repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:floor(length(ex.(driftPos))/4))...
            zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(1.25*length(ex.(driftPos))))...
            zeros(1,ex.finalFixation*ex.flipsPerSec)];
        
                 % make lonDriftPos for the red dot only condition
        ex.longDriftPos2 = [ex.longDriftPos2...
            repmat([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) repmat(ex.(driftPos),1,4) ... 
            ex.(driftPos)(1:floor(length(ex.(driftPos))/2)) zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(4.5*length(ex.(driftPos))))...
            zeros(1,ex.betweenBlocks*ex.flipsPerSec)],1,ex.numBlocks/length(ex.blockLength)-1)... %2*ones(1,e.blockLength) zeros(1,e.betweenBlocks)
            repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation)...
            repmat(ex.(driftPos),1,4) ex.(driftPos)(1:floor(length(ex.(driftPos))/2))...
            zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(4.5*length(ex.(driftPos))))...
            zeros(1,ex.finalFixation*ex.flipsPerSec)];
    end

end

%% Eyetracking parameters

if ET 
    EyelinkSetup(0);
    eye_used = Eyelink('EyeAvailable');
    ex.ShowRealTimeGaze = [  ]; % [] or [ something ]
    ex.nGazetoShow = [ 60 ]; % current~past N fixations
end
%% %%%% initial window - wait for backtick
% Screen(w, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
% Screen(w, 'Flip', 0);
% KbTriggerWait(53, deviceNumber);

DrawFormattedText(w,'Follow the oscillating visual phantom within the gap at the center of the screen \n\n as best as you can using the red dot as a guide, even after the red dot is gone. \n\n Press Space to start'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
%WaitSecs(10);  
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

gray = repmat(mean(squeeze(ex.rectLWave1(1,1,:))), [1,3]);

ex.fixCol1Grad = linspace(255,gray(1),90);% make red dot disapear in 90 flips = 1.5 sec ; logspace(log10(255),log10(gray(1)));
ex.fixCol2Grad = linspace(0,gray(1),90);

Screen('FillRect', w, gray);
    
if ET
    gcnt = 0; 
    EyeData.mx=nan(1,ceil(ex.totalTime*1/ex.flipsPerSec));
    EyeData.my=nan(1,ceil(ex.totalTime*1/ex.flipsPerSec));
    EyeData.ma=nan(1,ceil(ex.totalTime*1/ex.flipsPerSec));
    EyeData.FixDoneT = nan(1,ceil(ex.totalTime*1/ex.flipsPerSec));
    EyeData.gazeD = nan(1,ceil(ex.totalTime*1/ex.flipsPerSec));
    EyeData.Fixated = nan(1,ceil(ex.totalTime*1/ex.flipsPerSec));
end
%EyeStart = GetSecs(); %time we start caring about eyetracking
%currPosID = i;

onOffs = [diff(ex.longFormBlocks) 0];
flipCnt = 0;
cnt = 0;
tstartcnt = 0;
while n+1 < length(ex.allFlips)

    if ET
        run GetEyeDataLoic; %check eyetracker
    end

    [ex.longFormBlocks(n+1),ex.longFormFlicker(n+1)]

    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 0
        [VBLT, ex.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        ex.flipTime(n+1) = ex.startRun;
    else
        [VBLT, ex.flipTime(n+1), FlipT, missed] = Screen(w, 'Flip', ex.startRun + ex.allFlips(n+1) - slack);
    end
    
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ex.longFormBlocks(n+1) == 1 && ex.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
        thisCond = ex.longFormConds(n+1);
        
        %screen background color
        if strfind(conditions(thisCond).name{:}, 'Minbg')
            gray = repmat(min(min(squeeze(ex.rectLWave1(1,:,:)),[],1)), [1,3]);
        elseif strfind(conditions(thisCond).name{:}, 'Meanbg')
            gray = repmat(mean(squeeze(ex.rectLWave1(1,1,:))), [1,3]);
        end
        Screen('FillRect', w, gray);

        %draw red dot
        if strfind(conditions(thisCond).name{:}, 'Vel1')
            stillDotPhase = 'stillDotPhase1';
            driftPos = 'driftPos1';
            %longDriftPos = 'longDriftPos1';
            bLength = ex.blockLength(1);
        elseif strfind(conditions(thisCond).name{:}, 'Vel2')
            stillDotPhase = 'stillDotPhase2';
            driftPos = 'driftPos2';
            %longDriftPos = 'longDriftPos2';
            bLength = ex.blockLength(2);
        end
        if strfind(conditions(thisCond).name{:}, 'Rect')  
            
            ex.rectLRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xL,yL);
            ex.rectRRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xR,yR);
            
            if nnz(find(ex.rectLWaveID(n+1)))
                % top stim
                Screen('DrawTexture', w, ex.rectLWaveID(n+1),[],ex.rectLRect);
                % bottom stim
                Screen('DrawTexture', w, ex.rectRWaveID(n+1), [], ex.rectRRect);
            end
            %if nnz(ex.longDriftPos(n+1))
            flipCnt = flipCnt+1;
            %if flipCnt < 378
            timeOn = length([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:floor(length(ex.(driftPos))/4))]);
            if flipCnt <= timeOn - length(ex.fixCol1Grad)
                xOffset = ex.longDriftPos(n+1);
                Screen('FillOval', w,[255 0 0], [xc+xOffset-round(ex.fixSize/4) yc-round(ex.fixSize/4) xc+xOffset+round(ex.fixSize/4) yc+round(ex.fixSize/4)]);%black fixation solid circle
            elseif flipCnt > timeOn - length(ex.fixCol1Grad) && flipCnt < timeOn
                xOffset = ex.longDriftPos(n+1);
                Screen('FillOval', w,[ex.fixCol1Grad(flipCnt-timeOn+length(ex.fixCol1Grad)) ex.fixCol2Grad(flipCnt-timeOn+length(ex.fixCol1Grad)) ex.fixCol2Grad(flipCnt-timeOn+length(ex.fixCol1Grad))], [xc+xOffset-round(ex.fixSize/4) yc-round(ex.fixSize/4) xc+xOffset+round(ex.fixSize/4) yc+round(ex.fixSize/4)]);%black fixation solid circle
            elseif flipCnt == bLength*ex.flipsPerSec
                flipCnt = 0;
            end
            
        elseif strfind(conditions(thisCond).name{:}, 'RedDot')
            flipCnt = flipCnt+1;
            timeOn = length([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation)  repmat(ex.(driftPos),1,4) ex.(driftPos)(1:floor(length(ex.(driftPos))/2))]);
            if flipCnt < timeOn
                xOffset = ex.longDriftPos2(n+1);
                Screen('FillOval', w,[255 0 0], [xc+xOffset-round(ex.fixSize/4) yc-round(ex.fixSize/4) xc+xOffset+round(ex.fixSize/4) yc+round(ex.fixSize/4)]);%black fixation solid circle
            elseif flipCnt >= timeOn && flipCnt < (bLength)*ex.flipsPerSec
                Screen('FillOval', w,[gray], [xc+xOffset-round(ex.fixSize/4) yc-round(ex.fixSize/4) xc+xOffset+round(ex.fixSize/4) yc+round(ex.fixSize/4)]);%black fixation solid circle
            elseif flipCnt == (bLength)*ex.flipsPerSec
                flipCnt =0;
                clear timeOn
            end
        end
    end
    if nnz(onOffs(n+1)) == 1
        thisCond = ex.longFormConds(n+2);  
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
    
    if nnz(cnt) && mod(cnt,2) == 0 && GetSecs-time >= 1
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

ex.runTime = GetSecs - ex.startRun;


savedir = fullfile(ex.root,'data',sprintf('s%s_v13/',subject));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/s%s_smooth_pursuit_v13_date%s_fix',subject,num2str(ex.date)), '.mat'));
save(savename,'ex');

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
    savename = fullfile(savedir, strcat(sprintf('/s%s_smooth_pursuit_v13_date%s_fix_eyeDat',subject,num2str(ex.date)), '.mat'));
    save(savename, 'EyeData')
end 