function LD_phantom_grating_compare_v1(subject, session, debug)

%In this version, we add multiple velocities
% subject = 'Dave';                                                                                                                                                                                                                                                     
% session = 1;                                                                                                                           
% debug = 1;


ex.version = 'v1';
% global EyeData rect w xc yc eye_used 
%%%% resolution 
if debug == 1
    % eyetracking on (1) or off (0)
    ET = 0;
    ex.screenWidth = 53.1;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1600,900,60); % laptop 1920,1080/ 2880, 1800 ,0
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
    ET = 1;                                                                                                                             
    ex.screenWidth = 53.1;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; %23 in eye tracking                                                                                                                          room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1600,900,60); % ET room 1600,900,60
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('Return'))=1; % button box 3
responseKeys(KbName('ENTER'))=1; % button box 3

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

gray1 = [108.3751  108.3751  108.3751];%repmat(min(min(squeeze(ex.rectSWave(:,:,1,1)),[],1)), [1,3]);
gray2 = [146.6244  146.6244  146.6244];%repmat(max(max(squeeze(ex.rectSWave(:,:,1,1)),[],1)), [1,3]);
gray3 = [127.4992  127.4992  127.4992];%repmat(mean(squeeze(ex.rectSWave(1,:,1,1))), [1,3]);
grays = [gray1; gray2; gray3];
ex.stim.backgroundLum = grays;
%%%% 2D sine wave grating inducers properties
ex.stim.spatialFreqDeg = 0.5/2; %0.286;                                          % cycles per degree of visual angle
ex.stim.orientation = [90]; %[90 180];                                                % in degrees
ex.stim.degFromFix = .6;                                              % in degrees of visual angle
ex.stim.gaborHDeg = 8;                                                  % in degrees of visual angle
ex.stim.gaborWDeg = 16;
ex.stim.gapSizeDeg = 5;
% ex.stim.contrast = 0.15;%linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
% ex.stim.contrastMultiplicator = ex.stim.contrast/2;  % for sine wave 0.5 = 100% contrast, 0.2 = 40%
% ex.stim.contrastOffset = .5; %.5 .5 0];                                  % for procedural gabor
ex.stim.contrast = 0.15;
ex.stim.contrastOffset = [(ex.stim.backgroundLum(1,1)./255)./(1-ex.stim.contrast), (ex.stim.backgroundLum(2,1)./255)./(1+ex.stim.contrast), ex.stim.backgroundLum(3,1)./255];%+ex.stim.contrast/2;
ex.stim.luminanceRange = 2*ex.stim.contrast*ex.stim.contrastOffset;
ex.stim.contrastMultiplicator = ex.stim.luminanceRange./2;  % for sine wave

ex.stim.maxLum = 255*(ex.stim.contrastOffset+ex.stim.contrastMultiplicator);
ex.stim.minLum = 255*(ex.stim.contrastOffset-ex.stim.contrastMultiplicator);
ex.stim.contrast = (ex.stim.maxLum-ex.stim.minLum)./(ex.stim.maxLum+ex.stim.minLum);
ex.stim.cycPerSec = [1,1.4]; %drifting speed in cycles of grating per sec
ex.stim.cycles = [1 1]; %number of cycles per lap


%%%% sine wave grating timing (within block scale)
ex.stimDur = (ex.stim.cycles./ex.stim.cycPerSec)*2;        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
ex.initialFixation = 6;        % in seconds
ex.finalFixation = 2;          % in seconds
ex.blockLength = 10; %ex.trialFixation+ ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.betweenBlocks = 2;          % in seconds
ex.flipsPerSec = 60;  % 60;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%%%% conditions & layout (across blocks scale)

ex.conds = {'MinbgPairVel3','MaxbgPairVel3','MeanbgPairVel3',...%,'MinbgTopVel3','MinbgBotVel3', 'MeanbgTopVel3','MeanbgBotVel3'
    ... %,'MeanbgRedDotVel2',
    }; 
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 10;              % repetitions of each condition per run
ex.nTrials = ex.numConds*ex.repsPerRun;
ex.numBlocks = ex.numConds*ex.repsPerRun;

ex.condShuffle = [];
ex.locShuffle = [];
for i =1:ex.repsPerRun
    ex.condShuffle = [ex.condShuffle, Shuffle([1:ex.numConds])];
    ex.locShuffle = [ex.locShuffle, Shuffle([1,2])]; %2 locations (left and right)
end
ex.locShuffle = repmat(ex.locShuffle,1,2);
%%%% fixation 
ex.fixSizeDeg =  .25;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
ex.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 26;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

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
  load('/Users/tonglab/Desktop/monitor_calibration/425dell_22-12-09/phase2_photometry_22-12-09.mat');
  Screen('LoadNormalizedGammaTable', w, inverseCLUT);
end
%%%% timing optimization
frameInt = Screen('GetFlipInterval',w);
slack = frameInt/2;
frameRate =  1/frameInt;%Screen('NominalFrameRate',w);

flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.blockLength(1)];
% flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.stimDur(1)];%[0:frameInt*frameRate/ex.flipsPerSec:ex.blockLength(1)-ex.trialFixation]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
flipTimes = flipTimes(1:length(flipTimes)-1);
ex.stim.flipTimes = flipTimes;
ex.stim.tempPhase1 = nan(ex.numConds,1);
ex.stim.tempPhase2 = nan(ex.numConds,1);
ex.stim.oscillation1 = nan(ex.numConds,length(ex.stim.flipTimes));
ex.stim.oscillation2 = nan(ex.numConds,length(ex.stim.flipTimes));
ex.stim.phases = nan(ex.numConds,length(ex.stim.flipTimes));

clear c 
for c =1:ex.numConds

        ex.stim.tempPhase1(c) = rand(1,1)*2*pi;
        ex.stim.tempPhase2(c) = rand(1,1)*2*pi;
        ex.stim.oscillation1(c,:) = cos(2*pi*(1/ex.stimDur(1))*ex.stim.flipTimes+ex.stim.tempPhase1(c));
        ex.stim.oscillation2(c,:) = cos(2*pi*(1/ex.stimDur(2))*ex.stim.flipTimes+ex.stim.tempPhase2(c));
        ex.stim.spatialPhase = 90;
        ex.stim.phases(c,:) = (ex.stim.oscillation1(c,:).*180*ex.stim.cycles(1)+ ex.stim.oscillation2(c,:).*180*ex.stim.cycles(2))/2 + ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory

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
ex.gapSize = round(ex.stim.gapSizeDeg*ex.ppd);
ex.gaborHeight = round(ex.stim.gaborHDeg*ex.ppd);                 % in pixels, the size of our objects
ex.gaborWidth = round(ex.stim.gaborWDeg*ex.ppd);                 % in pixels, the size of our objects 
ex.rawGaborHeight = ex.gaborHeight;
ex.rawGaborWidth = ex.gaborWidth;

%% Create only one big sinewave grating image saved for each repetition and each condition

ex.rectSWave = nan(ex.rawGaborHeight,ex.rawGaborWidth,length(ex.stim.flipTimes), ex.numConds);
ex.rectSWaveID = nan(length(ex.stim.flipTimes),ex.numConds);
clear c r f
for c =1:ex.numConds %-1 %-1 because we only need images for the first 2 conditions
    %     for r = 1:ex.repsPerRun %only save the first image of each trial, that we will move during the trial
    for f = 1:length(ex.stim.flipTimes)
        phase = ex.stim.phases(c,f);
        ex.rectSWave(:,:,f,c) = makeSineGrating(ex.rawGaborHeight,ex.rawGaborWidth,ex.stim.spatialFreqDeg,...
            ex.stim.orientation,phase,ex.stim.contrastOffset(c),ex.stim.contrastMultiplicator(c),...
            ex.ppd);
        ex.rectSWaveID(f,c) = Screen('MakeTexture', w, squeeze(ex.rectSWave(:,:,f,c)));
    end
    %     end
end
%% Sine wave gratings locations (in the task loop since it changes)
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

%% stim location
xL = rect(3)/2-rect(3)/4; % % = stimulus located on the left side of the screen
xR = rect(3)/2+rect(3)/4; % = stimulus located on the right side of the screen

yT = rect(4)/2 - (ex.stim.gapSizeDeg+ex.stim.gaborHDeg)*ex.ppd/2; % stimulus located 4 degrees above screen center
yB = rect(4)/2+ (ex.stim.gapSizeDeg+ex.stim.gaborHDeg)*ex.ppd/2; % stimulus located 4 degrees below screen center
% yC = rect(4)/2; % stimulus located on screen center
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
KbQueueCreate(deviceNumber,responseKeys);
n = 1; %drift phase
c = 1; %condition
% t =1; %trial number
blockCnt = 1;
cntCond = zeros(length(ex.conds),1);

%%% initial fixation
if n == 1 && blockCnt == 1 %for first block
    ex.tasktstart = clock;
    ex.startRun = GetSecs();
    Screen('FillRect', w, gray1);
    Screen(w, 'Flip', 0);
    WaitSecs(ex.initialFixation);
end
%%% Launch the task
while(1) %n <= length(ex.trialFlips)
    KbQueueStart();
    thisCond = ex.condShuffle(c);
    condName = conditions(thisCond).name{:};
    Screen('FillRect', w, grays(thisCond,:))
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ex.locShuffle(c) == 1
        ex.rectTRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xR,yT);
        ex.rectBRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xR,yB);
    elseif ex.locShuffle(c) == 2
        ex.rectTRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xL,yT);
        ex.rectBRect =  CenterRectOnPoint([0 0 ex.rawGaborWidth ex.rawGaborHeight],xL,yB);
    end
    % stim
    Screen('DrawTexture', w, ex.rectSWaveID(n,thisCond),[],ex.rectTRect);
    Screen('DrawTexture', w, ex.rectSWaveID(n,thisCond),[],ex.rectBRect);
    
    
    %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 1
        [VBLT, ex.startTrial, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        flipTimes = ex.startTrial;
        
    else
        [VBLT,flipTime, FlipT, missed] = Screen(w, 'Flip',ex.startTrial + ex.stim.flipTimes(n) - slack); %,   %%% ex.flipTime(n,c)
        flipTimes = [flipTimes, flipTime];
        
    end
    
    
    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    if (pressed && ismember(find(firstPress,1), [KbName('Return') KbName('ENTER')]))
        n = 1;
        c = c+1;
    end
     KbQueueFlush();
    n = n+1;
    if (n == length(ex.rectSWaveID(:,1,1))+1) % for any other block , reset frame index when previous trial ends
        n = 1;
    end
    if c > ex.nTrials
        break;
    end
end



%%%%%%%%%%%%%%%%%%
% done! wrap up  %
%%%%%%%%%%%%%%%%%%

ex.runTime = GetSecs - ex.startRun;


savedir = fullfile(ex.root,'data',sprintf('on_between/s%s_%s/',subject,ex.version));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/s%s_percept_test_%s_date%s_fix',subject,ex.version,num2str(ex.date)), '.mat'));
%save(savename,'ex');
save(savename,'ex','-v7.3')

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;                                                                                                                          


end 