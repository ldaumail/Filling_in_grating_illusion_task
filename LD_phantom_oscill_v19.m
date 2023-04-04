function LD_phantom_oscill_v19(subject, session, debug)

%In this version, we add multiple velocities
% subject = 'Dave';                                                                                                                                                                                                                                                     
% session = 1;                                                                                                                           
% debug = 1;


ex.version = 'v19';
global EyeData rect w xc yc eye_used 
%%%% resolution 
if debug == 1
    % eyetracking on (1) or off (0)
    ET = 0;
    ex.screenWidth = 53.1;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),2880, 1800 ,0); % laptop 1920,1080/ 2880, 1800 ,0
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
    ET = 1;                                                                                                                             
    ex.screenWidth = 53.1;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; %23 in eye tracking                                                                                                                          room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),2880, 1800 ,0); % ET room 1600,900,60
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);

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
ex.stim.spatialFreqDeg = 0.5/2; %0.286;                                          % cycles per degree of visual angle
ex.stim.contrast = 0.3 ;                                                 % in %, maybe??
ex.stim.orientation = [90]; %[90 180];                                                % in degrees
ex.stim.degFromFix = .6;                                              % in degrees of visual angle
ex.stim.gaborHDeg = 8;                                                  % in degrees of visual angle
ex.stim.gaborWDeg = 8;
%ex.stim.rectGaborWDeg = 8;
ex.stim.contrast = 0.15;%linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
ex.stim.contrastMultiplicator = ex.stim.contrast/2;  % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.stim.contrastOffset = .5; %.5 .5 0];                                  % for procedural gabor
ex.stim.cycPerSec = [1.13*1/2,1.13*3/2]; % try multiple speeds
%ex.stim.sumCycPerSec = [0.3,0.7];
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
ex.flipsPerSec = 60;  % 60;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%e.numBlocks = 12;  % 6 for single and 6 for pair...

%%%% conditions & layout (across blocks scale)

ex.conds = {'MeanbgRedDotVel1',...%'RectMinbgSingleVel1'
    'RectMinbgSp8Vel1','RectSumSp8Vel1' ...   %,
    'MeanbgRedDotVel2','RectMinbgSp8Vel2','RectSumSp8Vel2'...%'RectMinbgSingleVel2',
    }; 
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 20;              % repetitions of each condition per run
ex.numBlocks = ex.numConds*ex.repsPerRun;

ex.condShuffle = [];
for i =1:ex.repsPerRun
    ex.condShuffle = [ex.condShuffle, Shuffle([1:ex.numConds/2])];
end
for i =1:ex.repsPerRun
    ex.condShuffle = [ex.condShuffle, Shuffle([ex.numConds/2+1:ex.numConds])];
end

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
%%%% fixation 
ex.fixSizeDeg =  .25;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
ex.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 26;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

ex.onSecs = [ones(1,ex.blockLength(t)) zeros(1,ex.betweenBlocks)];
ex.longFormBlocks = Expand(ex.onSecs,ex.flipsPerSec,1); %1 when block, 0 when between block
length(ex.longFormBlocks)

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
  %load gamma correction file
%   ex.whichCLUT = '/Users/tonglab/Desktop/monitor_calibration/425dell_22-12-09/phase2_photometry_22-12-09.mat';
%   load(ex.whichCLUT);
  load('/Users/tonglab/Desktop/monitor_calibration/425dell_22-12-09/phase2_photometry_22-12-09.mat');
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

%oscillation = sprintf('oscillation%d',t);
%ex.stim.(oscillation) = cos(2*pi*(1/ex.stimDur(t))*flipTimes+ex.stim.tempPhase);
flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.stimDur(1)]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
flipTimes = flipTimes(1:length(flipTimes)-1);
ex.sum.phases1 = nan(length(flipTimes), ex.repsPerRun);
ex.sum.phases2 = nan(length(flipTimes), ex.repsPerRun);
oscillation = 'oscillation1';
for r = 1:ex.repsPerRun
    ex.sum.spatialPhase1 = randi(360); %this makes the grating phase so that the center of the dark stripe is on the center of the screen
    ex.sum.spatialPhase2 = randi(360);
    ex.sum.phases1(:,r) = ex.stim.(oscillation).*180*ex.stim.cycles(1)+ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory
    ex.sum.phases2(:,r) = ex.stim.(oscillation).*180*ex.stim.cycles(2)+ex.stim.spatialPhase;
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
rectLWaveID = nan(length(ex.longFormBlocks),length(ex.stim.cycPerSec));
rectRWaveID = nan(length(ex.longFormBlocks),length(ex.stim.cycPerSec));
clear t

for t =1:length(ex.stim.cycPerSec)
    flipt =sprintf('flipTimes%d',t);
    flipTimes = ex.(flipt);
    
    %rect
    rectLWave = sprintf('rectLWave%d',t);
    rectRWave = sprintf('rectRWave%d',t);

    ex.(rectLWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2);
    ex.(rectRWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2);
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
        ex.(rectLWave)(f,:,:) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,LPhase,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
            ex.ppd);
        ex.(rectRWave)(f,:,:) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,RPhase,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
            ex.ppd);

        %                 figure();
        %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
        %
        tmprectLWaveID(f,1) = Screen('MakeTexture', w, squeeze(ex.(rectLWave)(f,:,:)));
        tmprectRWaveID(f,1) = Screen('MakeTexture', w, squeeze(ex.(rectRWave)(f,:,:)));
 
    end
    %% extend stimulus matrix to include the same total number of flips as the whole experiment
        rectLWaveID(:,t) = [zeros(length(ex.stillDotPhase),1); ...
            repmat(tmprectLWaveID(:),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectLWaveID(:))),1);...
            tmprectLWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectLWaveID(:))));...
            zeros(ex.betweenBlocks*ex.flipsPerSec,1)];
        
        rectRWaveID(:,t) = [zeros(length(ex.stillDotPhase),1);...
            repmat(tmprectRWaveID(:),floor((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec/length(tmprectRWaveID(:))),1);...
            tmprectRWaveID(1:mod((ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec,length(tmprectRWaveID(:))));...
            zeros(ex.betweenBlocks*ex.flipsPerSec,1)];

end
ex.rectLWaveID = rectLWaveID;
ex.rectRWaveID = rectRWaveID;

%% Do the same for the summed gratings condition
rectSWaveID = nan(length(ex.longFormBlocks),ex.repsPerRun);
%rectBWaveID = nan(length(ex.longFormBlocks),1);
clear r
flipTimes = ex.flipTimes1;
ex.rectSWave = nan(ex.repsPerRun,length(flipTimes),ex.gaborHeight,ex.gaborWidth*2);

 for r =1:ex.repsPerRun
    slowPhases = ex.sum.phases1(:,r);
    fastPhases = ex.sum.phases2(:,r);
    %       clear tmprectLWaveID tmprectRWaveID
    for f = 1:length(flipTimes)
        
        phase1 = slowPhases(f);
        phase2 = fastPhases(f);
        %ih in pixels %iw in pixels %spatial freq in cycles per dva
        %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
        %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
        %background color (unused if the grating is not an annulus)
        
        % rect
        ex.rectSWave(r,f,:,:) = (makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,phase1,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
            ex.ppd)+makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,phase2,ex.stim.contrastOffset(1),ex.stim.contrastMultiplicator,...
            ex.ppd))/2;

        %                 figure();
        %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
        %
        tmprectSWaveID(f,1) = Screen('MakeTexture', w, squeeze(ex.rectSWave(r,f,:,:)));
    end
    %% extend stimulus matrix to include the same total number of flips as the whole experiment
        rectSWaveID(:,r) = [zeros(length(ex.stillDotPhase),1); ...
            repmat(tmprectSWaveID(:),floor((ex.blockLength(1)-ex.trialFixation)*ex.flipsPerSec/length(tmprectSWaveID(:))),1);...
            tmprectSWaveID(1:mod((ex.blockLength(1)-ex.trialFixation)*ex.flipsPerSec,length(tmprectSWaveID(:))));...
            zeros(ex.betweenBlocks*ex.flipsPerSec,1)];
 

end
ex.rectSWaveID = rectSWaveID;


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
    longDriftPos = sprintf('longDriftPos%d',t);
    
    ex.(longDriftPos) = [repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:length(ex.(driftPos))/4) ... %drift for 1.25 (5/4) periods of oscilation
        zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor((5/4)*length(ex.(driftPos)))) ...
        zeros(1,ex.betweenBlocks*ex.flipsPerSec)];
    
    % make lonDriftPos for the red dot only condition
    longDriftPosRed = sprintf('longDriftPosRed%d',t);
    ex.(longDriftPosRed) =[repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) repmat(ex.(driftPos),1,4) ...
        ex.(driftPos)(1:ex.flipsPerSec*(ex.blockLength(t)-ex.trialFixation)-length(repmat(ex.(driftPos),1,4))) ... %zeros(1,(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-floor(4.5*length(ex.(driftPos))))...
        zeros(1,ex.betweenBlocks*ex.flipsPerSec)];
end


%% Eyetracking parameters

if ET 
    EyelinkSetup(0);
    eye_used = Eyelink('EyeAvailable');
%     ex.ShowRealTimeGaze = [  ]; % [] or [ something ]
%     ex.nGazetoShow = [ 60 ]; % current~past N fixations
end
%% %%%% initial window - wait for backtick
DrawFormattedText(w,'Follow the oscillating grating or visual phantom within the gap at the center of the screen \n\n as best as you can using the red dot as a guide, even after the red dot is gone. \n\n Do your best not to blink during a trial. \n\n Press Space to start'... % :  '...
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
gray = repmat(mean(squeeze(ex.rectLWave1(1,1,:))), [1,3]);
n=1;
blockCnt = 1;
cnt = 0; %stim onset/ stime offset count
cntS = 0; %count number of summed gratings condition trials
onOffs = [diff(ex.longFormBlocks) 0];
bLength = ex.blockLength(1);
ex.flipTime = nan(length(ex.trialFlips),length(ex.condShuffle));
%%% initial fixation
if n == 1 && blockCnt == 1 %for first block
    ex.tasktstart = clock;
    ex.startRun = GetSecs();
    Screen('FillRect', w, gray);
    Screen(w, 'Flip', 0);
    WaitSecs(ex.initialFixation);
end
%%% Launch the task
for c =1:length(ex.condShuffle)
    cnt = cnt+1;
    thisCond = ex.condShuffle(c);
    condName = conditions(thisCond).name{:};
    %%for each condition, we specify the parameters values before we flip
    %%over the gratings phases
    %screen background color
    if strfind(condName, 'Minbg') | strfind(condName, 'Sum') %contains(condName, 'Minbg')
        gray = repmat(min(min(squeeze(ex.rectLWave1(1,:,:)),[],1)), [1,3]);
        ex.fixCol1Grad = linspace(255,gray(1),90);% make red dot disapear in 90 flips = 1.5 sec ; logspace(log10(255),log10(gray(1)));
        ex.fixCol2Grad = linspace(0,gray(1),90);
        
    elseif strfind(condName, 'Meanbg')
        gray = repmat(mean(squeeze(ex.rectLWave1(1,1,:))), [1,3]);
        ex.fixCol1Grad = linspace(255,gray(1),90);% make red dot disapear in 90 flips = 1.5 sec ; logspace(log10(255),log10(gray(1)));
        ex.fixCol2Grad = linspace(0,gray(1),90);
    end
    %draw guide red dot
    if strfind(condName, 'Vel1')
        stillDotPhase = 'stillDotPhase1';
        driftPos = 'driftPos1';
        longDriftPos = 'longDriftPos1';
        t =1;
    elseif strfind(condName, 'Vel2')
        stillDotPhase = 'stillDotPhase2';
        driftPos = 'driftPos2';
        longDriftPos = 'longDriftPos2';
        t = 2;
    end
    if strfind(condName, 'Sp8') | strfind(condName, 'Sum')% is there gonna be a grating inducers pair? If yes, indicate the location
        ex.rectLRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xL,yL);
        ex.rectRRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xR,yR);
        timeOn = length([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation) ex.(driftPos) ex.(driftPos)(1:length(ex.(driftPos))/4)]);
        guideFlipCnt = 1;
    elseif strfind(condName, 'RedDot')
        timeOn = length([repmat(ex.(stillDotPhase),1,ex.flipsPerSec*ex.trialFixation)  repmat(ex.(driftPos),1,4) ex.(driftPos)(1:(ex.blockLength(t)-ex.trialFixation)*ex.flipsPerSec-length(repmat(ex.(driftPos),1,4)))]);
        longDriftPosRed = sprintf('longDriftPosRed%d',t);
        redFlipCnt = 1;
    end
    
    %flip through the block and following between block time
    while n <= length(ex.trialFlips)
        ex.longFormBlocks(n)
        %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Screen('FillRect', w, gray);
        
        if strfind(condName, 'Sp8')
            if nnz(find(ex.rectLWaveID(n,t)))
                % top stim
                Screen('DrawTexture', w, ex.rectLWaveID(n,t),[],ex.rectLRect);
                % bottom stim
                Screen('DrawTexture', w, ex.rectRWaveID(n,t), [], ex.rectRRect);
            end
            if guideFlipCnt <= timeOn - length(ex.fixCol1Grad)
                xOffset = ex.(longDriftPos)(n);
                Screen('FillOval', w,[255 0 0], [xc+xOffset-round(ex.fixSize/2) yc-round(ex.fixSize/2) xc+xOffset+round(ex.fixSize/2) yc+round(ex.fixSize/2)]);%black fixation solid circle
            elseif guideFlipCnt > timeOn - length(ex.fixCol1Grad) && guideFlipCnt < timeOn
                xOffset = ex.(longDriftPos)(n);
                Screen('FillOval', w,[ex.fixCol1Grad(guideFlipCnt-timeOn+length(ex.fixCol1Grad)) ex.fixCol2Grad(guideFlipCnt-timeOn+length(ex.fixCol1Grad)) ex.fixCol2Grad(guideFlipCnt-timeOn+length(ex.fixCol1Grad))], [xc+xOffset-round(ex.fixSize/2) yc-round(ex.fixSize/2) xc+xOffset+round(ex.fixSize/2) yc+round(ex.fixSize/2)]);% fixation solid circle
            elseif guideFlipCnt == length(ex.trialFlips)%bLength*ex.flipsPerSec+1 %make sure to reset the flipCnt once the trial ended, so that the red dot doe not reappear before the end of the trial
                guideFlipCnt = 0;
            end
            guideFlipCnt = guideFlipCnt+1;
            
        elseif strfind(condName, 'Sum')
            cntS = cntS+1;
            if nnz(find(ex.rectSWaveID(n,cntS)))
                % top stim
                Screen('DrawTexture', w, ex.rectSWaveID(n,cntS),[],ex.rectLRect);
                % bottom stim
                Screen('DrawTexture', w, ex.rectSWaveID(n,cntS), [], ex.rectRRect);
            end
            if guideFlipCnt <= timeOn - length(ex.fixCol1Grad)
                xOffset = ex.(longDriftPos)(n);
                Screen('FillOval', w,[255 0 0], [xc+xOffset-round(ex.fixSize/2) yc-round(ex.fixSize/2) xc+xOffset+round(ex.fixSize/2) yc+round(ex.fixSize/2)]);%black fixation solid circle
            elseif guideFlipCnt > timeOn - length(ex.fixCol1Grad) && guideFlipCnt < timeOn
                xOffset = ex.(longDriftPos)(n);
                Screen('FillOval', w,[ex.fixCol1Grad(guideFlipCnt-timeOn+length(ex.fixCol1Grad)) ex.fixCol2Grad(guideFlipCnt-timeOn+length(ex.fixCol1Grad)) ex.fixCol2Grad(guideFlipCnt-timeOn+length(ex.fixCol1Grad))], [xc+xOffset-round(ex.fixSize/2) yc-round(ex.fixSize/2) xc+xOffset+round(ex.fixSize/2) yc+round(ex.fixSize/2)]);% fixation solid circle
            elseif guideFlipCnt == length(ex.trialFlips)%bLength*ex.flipsPerSec+1 %make sure to reset the flipCnt once the trial ended, so that the red dot doe not reappear before the end of the trial
                guideFlipCnt = 0;
            end
            guideFlipCnt = guideFlipCnt+1;
            
        elseif strfind(condName, 'RedDot')
            if redFlipCnt < timeOn
                xOffset = ex.(longDriftPosRed)(n);
                Screen('FillOval', w,[255 0 0], [xc+xOffset-round(ex.fixSize/2) yc-round(ex.fixSize/2) xc+xOffset+round(ex.fixSize/2) yc+round(ex.fixSize/2)]);%black fixation solid circle
            elseif redFlipCnt >= timeOn && redFlipCnt < (bLength)*ex.flipsPerSec
                Screen('FillOval', w,[gray], [xc+xOffset-round(ex.fixSize/2) yc-round(ex.fixSize/2) xc+xOffset+round(ex.fixSize/2) yc+round(ex.fixSize/2)]);%black fixation solid circle
            elseif redFlipCnt == length(ex.trialFlips)%(bLength)*ex.flipsPerSec+1 %make sure to reset the flipCnt once the trial ended, so that the red dot doe not reappear before the end of the trial
                redFlipCnt =0;
                %clear timeOn
            end
            redFlipCnt = redFlipCnt+1;
        end
        
        %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if n == 1
            [VBLT, ex.startTrial, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
            flipTimes = ex.startTrial;
            %%%%% send messages to Eyelink for event times
            if ET == 1
                Eyelink('Message', char(sprintf('Cond %s', condName)));
                Eyelink('Message', 'TRIALID %d', c);
                Eyelink('Message', 'STIM_ONSET');
            end
        else
            [VBLT,flipTime, FlipT, missed] = Screen(w, 'Flip',ex.startTrial + ex.trialFlips(n) - slack); %,   %%% ex.flipTime(n,c)
            flipTimes = [flipTimes, flipTime];
            if ET == 1 && n == bLength*ex.flipsPerSec+1
                Eyelink('Message', 'STIM_OFFSET');
            end
            
        end
        
        if nnz(onOffs(n)) == 1
            time = GetSecs;
            cnt = cnt+1;
        end
        if (cnt/2 == 1 && GetSecs-time >= 1) && c ~= length(ex.condShuffle)
            DrawFormattedText(w,'Press Space whenever you feel ready'... % : press 1 as soon as letter J appears on the screen,\n\n and press 2 as soon as letter K appears on the screen. \n\n Press Space to start'...
                ,'center', 'center',[0 0 0]);
            Screen(w, 'Flip', 0);
            KbTriggerWait(KbName('Space'), deviceNumber);
            cnt = 0;
        end
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


savedir = fullfile(ex.root,'data',sprintf('s%s_%s/',subject,ex.version));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/s%s_smooth_pursuit_%s_date%s_fix',subject,ex.version,num2str(ex.date)), '.mat'));
%save(savename,'ex');
save(savename,'ex','-v7.3')

ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;                                                                                                                           


if ET
    Eyelink('StopRecording');   
    Eyelink('CloseFile');
        %download edf file
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
    
    Eyelink('Shutdown');
%     savename = fullfile(savedir, strcat(sprintf('/s%s_smooth_pursuit_%s_date%s_fix_eyeDat',subject,ex.version,num2str(ex.date)), '.mat'));
%     save(savename, 'EyeData')
end 