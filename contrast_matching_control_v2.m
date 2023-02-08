function contrast_matching_control_v2(subject, session, debug) 

%%contrast matching control
%Loic 01252023
%In this version, we add multiple velocities
subject = 'Dave';                                                                                                                                                                                                                                                     
session = 1;                                                                                                                           
debug = 1;

ex.version = 'v13';
global EyeData rect w xc yc %eye_used
%%%% resolution 
if debug == 1
    % eyetracking on (1) or off (0)
    ET = 0;
    ex.screenWidth = 53.1;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop 1920,1080
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
    ET = 1;                                                                                                                             
    ex.screenWidth = 53.1;             % in cm; % 16 in eye tracking room 425%laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; %23 in eye tracking                                                                                                                          room 425 3Tb/office=43, miniHelm=57;
    ex.resolution = SetResolution(max(Screen('Screens')),1600,900,60); % scanner
    ex.gammaCorrection = 1;       % make sure this = 1 when you're at the scanner!
end

%%%% keyboard
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);
responseKeys = zeros(1,256);
responseKeys(KbName('1'))=1; % button box 1
responseKeys(KbName('2'))=1; % button box 2
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2
responseKeys(KbName('Return'))=1; % button box 3
responseKeys(KbName('ENTER'))=1; % button box 3

Screen('Preference', 'SkipSyncTests', 0);

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


%%%% 2D sine wave grating properties
ex.stim.spatialFreqDeg = 0.5/2; %0.286;                                          % cycles per degree of visual angle
ex.stim.orientation = [90]; %[90 180];                                                % in degrees
ex.stim.degFromFix = .6;                                              % in degrees of visual angle
ex.stim.gaborHDeg = 8;                                                  % in degrees of visual angle
ex.stim.gaborWDeg = 8;
%ex.stim.rectGaborWDeg = 8;
ex.stim.contrast = 0.15 ;                                                 % in %, maybe??
ex.stim.contrastMultiplicator = ex.stim.contrast/2;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.match.contrastMults = linspace(0,ex.stim.contrastMultiplicator/4,20);%linspace(0,0.075,20);                            %contrast multiplicators for the adjustable sine wave from 0 to 0.6 contrast
ex.stim.contrastOffset = .5;                                  % for procedural gabor
ex.match.contrastOffset = ex.stim.contrastOffset-ex.stim.contrastMultiplicator-ex.match.contrastMults; 
%-ex.stim.contrastMultiplicator because we want the maximum luminance value to match the background luminance (which corresponds to ex.stim.contrastOffset-ex.stim.contrastMultiplicator, the dark stim grating luminance)
%thus our baseline luminance is ex.stim.contrastOffset-ex.stim.contrastMultiplicator
%In addition, we subtract -ex.match.contrastMults to make sure that
%for every matching contrast level, the maximum luminance remains equal to
%the dark stripe luminance, thus we lower the baseline luminance according to each contrast level.
ex.stim.cycPerSec = [1.13*1/2,1.13*3/2]; % try multiple speeds
ex.stim.motionRate = ex.stim.cycPerSec.*360;                                          % 1.13 cycles per second = 360 deg of phase *1.13 per sec
ex.stim.cycles =[1, 3]; %number of cycles shifted per lap  (modified to half the number of cycles per lap)

%%%% sine wave grating timing (within block scale)
ex.initialFixation = 6;        % in seconds
ex.finalFixation = 2;          % in seconds
ex.trialFixation = 0;          % in seconds
ex.stimDur = (ex.stim.cycles./ex.stim.cycPerSec)*2;        % in seconds. 1.77 sec refers to sine wave grating 1.77 = 2cycles/1.13cyc.sec-1 mutiplied by 2 for back and forth
ex.stimsPerBlock = 1;      % number of back-and-forth laps of the stimulus drift
ex.blockLength = ex.trialFixation+ ceil(ex.stimDur*ex.stimsPerBlock);           % in seconds
ex.betweenBlocks = 1;          % in seconds
ex.flipsPerSec = 60;  % 12;         % number of phase changes we want from the visual stimulus, and thus the number of times we want to change visual stimulation on the screen
ex.flipWin = 1/ex.flipsPerSec;         % in seconds then actually in 1 sec the stimuli will change 12 times 
%e.nTrials = 12;  % 6 for single and 6 for pair...

%%%% conditions & layout (across blocks scale)

ex.conds = {'MinbgMatchVel1'; 'MinbgMatchVel2'};
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 20;              % repetitions of each condition per run
ex.nTrials = ex.numConds*ex.repsPerRun;
% condShuffle1 = Shuffle(repmat([1:ex.numConds/2],1,ex.repsPerRun)); % %e.stimsPerBlock make same number of blocks with each condition, randomize order
% condShuffle2 = Shuffle(repmat([ex.numConds/2+1:ex.numConds],1,ex.repsPerRun));
% 
% ex.condShuffle = [condShuffle1 condShuffle2];

%% This section is only required for preallocating data for eye movements
% ex.totalTime = [];
% for t =1:length(ex.blockLength) %there is a different block length for every drifting speed
%     if t == 1
%         ex.totalTime = sum([ex.totalTime, ex.initialFixation + (ex.nTrials/length(ex.blockLength) * (ex.blockLength(t) + ex.betweenBlocks))]);
%     elseif t <length(ex.blockLength) && t > 1
%              ex.totalTime = sum([ex.totalTime, (ex.nTrials/length(ex.blockLength) * (ex.blockLength(t) + ex.betweenBlocks))]); 
%     elseif t == length(ex.blockLength)
%         ex.totalTime = sum([ex.totalTime, ((ex.nTrials/length(ex.blockLength)-1) * (ex.blockLength(t) + ex.betweenBlocks)) + ex.blockLength(t) + ex.finalFixation]);
%     end
% end
% ex.allFlips = (0:ex.flipWin:ex.totalTime);

%% fixation 
ex.fixSizeDeg =  .5;            % in degrees, the size of the biggest white dot in the fixation

%%%% screen
ex.backgroundColor = [127 127 127];%[108.3760 108.3760 108.3760];%;  % color based on minimum gating luminance 
ex.fontSize = 14; %26;

%% %%%%%%%%%%%%%%%%%
   % timing model  %
   %%%%%%%%%%%%%%%%%

ex.onSecs = ones(1,ex.blockLength(1));
ex.longFormBlocks = Expand(ex.onSecs,ex.flipsPerSec,1); %1 when block, 0 when between block
ex.longFormFlicker = repmat(ones(1,1),1,length(ex.longFormBlocks)); %1 all the way to ensure flip at every time selected
length(ex.longFormBlocks)

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
     %  [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
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
for s =1:length(ex.stim.cycPerSec)
    flipTimes = [0:frameInt*frameRate/ex.flipsPerSec:ex.stimDur(s)]; %multiply frameInt by 60/12 = 5 to flip the image every 5 frames
    flipTimes = flipTimes(1:length(flipTimes)-1);
    flipt = sprintf('flipTimes%d',s);
    ex.(flipt) = flipTimes;
    
    %%% Still red dot for 1 sec before trial starts
    ex.stillDotPhase = zeros(1,ex.flipsPerSec*ex.trialFixation);
    
    %%% Stimulus phases
    % (in a linear/triangle scenario) ex.stim.dphasePerFlip = ex.stim.motionRate*frameInt * frameRate/ex.flipsPerSec; %degrees per flip, here we multiply by 60/12 =5 to move the phase 5 times more after each flip since we have 5 times less flips than frame refresh
    ex.stim.tempPhase =pi; % pi/2; %this will make the oscillation start on the fast segment of the oscillation, with red dot on the center of the screen
    
    oscillation = sprintf('oscillation%d',s);
    ex.stim.(oscillation) = cos(2*pi*(1/ex.stimDur(s))*flipTimes+ex.stim.tempPhase);
    ex.stim.spatialPhase = 90; %this makes the grating phase so that the center of the dark stripe is on the center of the screen
    phases = sprintf('phases%d',s);
    ex.stim.(phases) = ex.stim.(oscillation).*180*ex.stim.cycles(s)+ex.stim.spatialPhase; %./ex.stimDur-2*pi*flipTimes./ex.stimDur make it oscillatory
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
for o =1:length(ex.stim.orientation)
    rectLWaveIDO = nan(length(flipTimes),length(ex.stimDur));
    rectRWaveIDO = nan(length(flipTimes),length(ex.stimDur));
    clear t
    
    for s =1:length(ex.stim.cycPerSec)
        flipt =sprintf('flipTimes%d',s);
        flipTimes = ex.(flipt);
        
        %rect
        rectLWave = sprintf('rectLWave%d',s);
        rectRWave = sprintf('rectRWave%d',s);
        
        ex.(rectLWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2,length(ex.stim.orientation));
        ex.(rectRWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2,length(ex.stim.orientation));
        
        phaseNum = sprintf('phases%d',s);
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
                ex.stim.orientation(o),LPhase,ex.stim.contrastOffset,ex.stim.contrastMultiplicator,...
                ex.ppd);
            ex.(rectRWave)(f,:,:,o) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                ex.stim.orientation(o),RPhase,ex.stim.contrastOffset,ex.stim.contrastMultiplicator,...
                ex.ppd);
            %                 figure();
            %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
            %
            tmprectLWaveID(f,o) = Screen('MakeTexture', w, squeeze(ex.(rectLWave)(f,:,:,o)));
            tmprectRWaveID(f,o) = Screen('MakeTexture', w, squeeze(ex.(rectRWave)(f,:,:,o)));
            
        end
        
        %% extend stimulus matrix to include the same total number of flips as the whole experiment
        
        if s<length(ex.stimDur)
            %rect
            rectLWaveIDO(:,s) = tmprectLWaveID(:,o);
            
            rectRWaveIDO(:,s) = tmprectRWaveID(:,o);
            
        elseif s == length(ex.stimDur)
            
            %rect
            rectLWaveIDO(:,s) = tmprectLWaveID(:,o);
            ex.rectLWaveID(1:length(rectLWaveIDO),1:2,o) = rectLWaveIDO;
            
            rectRWaveIDO(:,s) = tmprectRWaveID(:,o);
            ex.rectRWaveID(1:length(rectRWaveIDO),1:2,o) = rectRWaveIDO;
        end
    end
end

%% Create contrast adjustable sine wave grating

clear t
for c =1:length(ex.match.contrastMults)
    rectCWaveIDO = nan(length(flipTimes),length(ex.stimDur)) ;
    for s =1:length(ex.stim.cycPerSec)
        flipt =sprintf('flipTimes%d',s);
        flipTimes = ex.(flipt);
        
        %rect
        rectCWave = sprintf('rectCWave%d',s);
        
        ex.(rectCWave) = nan(length(flipTimes),ex.gaborHeight,ex.gaborWidth*2,length(ex.match.contrastMults));
        
        phaseNum = sprintf('phases%d',s);
        phases = ex.stim.(phaseNum) ;
        
        %       clear tmprectLWaveID tmprectRWaveID
        for f = 1:length(flipTimes)
            
            CPhase = phases(f);
            %ih in pixels %iw in pixels %spatial freq in cycles per dva
            %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
            %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
            %background color (unused if the grating is not an annulus)
             
            % rect
            ex.(rectCWave)(f,:,:,c) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                ex.stim.orientation(1),CPhase,ex.match.contrastOffset(c) ,ex.match.contrastMults(c),...
                ex.ppd);%-min(min(squeeze(ex.rectLWave1(1,:,:)),[],1));  %center contrast level to dark luminance background so that the max luminance matches the dark luminance background
            %                 figure();
            %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
            %
            tmprectCWaveID(f,c) = Screen('MakeTexture', w, squeeze(ex.(rectCWave)(f,:,:,c)));
            
        end
        
        %% extend stimulus matrix to include the same total number of flips as the whole experiment
        
        if s<length(ex.stimDur)
            %rect
            rectCWaveIDO(:,s) = tmprectCWaveID(:,c);            
        elseif s == length(ex.stimDur)
            %rect
            rectCWaveIDO(:,s) = tmprectCWaveID(:,c);
            ex.rectCWaveID(1:length(rectCWaveIDO),1:2,c) = rectCWaveIDO;
        end
    end
end
%% Sine wave gratings locations (in the task loop since it changes)
xc = rect(3)/2; % rect and center, with the flexibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+e.vertOffset;

xL = rect(3)/2-rect(3)/4; % % = stimulus center located on the horizontal center of the screen
xR = rect(3)/2-rect(3)/4; % = stimulus center located on the horizontal center of the screen
xC = rect(3)/2+rect(3)/4; % = stimulus located on the right side of the screen

yL = rect(4)/2 - ex.stim.gaborHDeg*ex.ppd+0*ex.ppd; % stimulus located 4 degrees above screen center
yR = rect(4)/2+ ex.stim.gaborHDeg*ex.ppd-0*ex.ppd; % stimulus located 4 degrees below screen center
yC = rect(4)/2; % stimulus located on screen center


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

DrawFormattedText(w,'Match the contrast level of the oscillating visual phantom within the gap \n\n at the center-left side of the screen \n\n with that of the sinewave grating on the right side of the screen. \n\n To increase contrast press 1, to decrease contrast press 2  \n\n Press Space to start'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
%  WaitSecs(2);  
KbTriggerWait(KbName('Space'), deviceNumber);

% %%%% response listening 

KbQueueCreate(deviceNumber,responseKeys);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% START task TASK/FLIPPING

gray = repmat(mean(squeeze(ex.rectLWave1(1,1,:))), [1,3]);

ex.fixCol1Grad = linspace(255,gray(1),90);% make red dot disapear in 90 flips = 1.5 sec ; logspace(log10(255),log10(gray(1)));
ex.fixCol2Grad = linspace(0,gray(1),90);

Screen('FillRect', w, gray);
    
if ET
    gcnt = 0; 
    EyeData.mx=nan(1,1);
    EyeData.my=nan(1,1);
    EyeData.ma=nan(1,1);
    EyeData.FixDoneT = nan(1,1);
    EyeData.gazeD = nan(1,1);
    EyeData.Fixated = nan(1,1);
end

n=1; %initialize at first flip
cont = 10; %initial contrast level index for adjustable stimulus
ex.cont = [];
cnt = 0;
tstartcnt = 0;
tend = [];
s =1; %stimulus velocity
t =1; %trial number
while(1) %n+1 < length(ex.allFlips)
    KbQueueStart();
    if ET
        run GetEyeDataLoic; %check eyetracker
    end

    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %if ex.longFormBlocks(n+1) == 1 && ex.longFormFlicker(n+1) > 0 % zeros correspond to IBI, in which case we skip this next section
    
    % thisCond = ex.longFormConds(n+1);
    %screen background color
    
    gray = repmat(min(min(squeeze(ex.rectLWave1(1,:,:)),[],1)), [1,3]);
    Screen('FillRect', w, gray);
    
    ex.rectLRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xL,yL);
    ex.rectRRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xR,yR);
    ex.rectCRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xC,yC);
    
    if nnz(find(ex.rectLWaveID(n,1)))
        % top stim
        Screen('DrawTexture', w, ex.rectLWaveID(n,1),[],ex.rectLRect);
        % bottom stim
        Screen('DrawTexture', w, ex.rectRWaveID(n,1 ),[], ex.rectRRect);
        % side stim
        Screen('DrawTexture', w, ex.rectCWaveID(n,s,cont), [], ex.rectCRect);
    end
        %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n == 1
        [VBLT, ex.startRun, FlipT, missed] = Screen(w, 'Flip', 0);%[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed] = Screen('Flip', windowPtr [, when] [, dontclear]...
        ex.flipTime(n) = ex.startRun;
    else
        [VBLT, ex.flipTime(n), FlipT, missed] = Screen(w, 'Flip', ex.startRun + ex.flipTimes1(n) - slack);
    end
    
    KbQueueStop();
    [pressed, firstPress]= KbQueueCheck();
    if (pressed && ismember(find(firstPress,1), [KbName('1') KbName('2') KbName('1!') KbName('2@')]))
        if (find(firstPress,1) == KbName('1')|find(firstPress,1) == KbName('1!'))
            cont = cont+1; % increase contrast
        elseif (find(firstPress,1) == KbName('2')|find(firstPress,1) == KbName('2@'))
            cont = cont-1; % decrease contrast
        end
        if cont == 0 || cont == length(ex.match.contrastMults)+1
            cont = 10;
        end
    elseif (pressed && ismember(find(firstPress,1), [KbName('Return') KbName('ENTER')]))
        ex.cont = [ex.cont ex.match.contrastMults(cont)];
        n = 0;
        t = t+1;
        cont = 10; 
        if  ET == 1%isempty(find(mod(tstartcnt,2)))
            Eyelink('Message', 'STIM_OFFSET');
            Eyelink('Message', 'TRIALID %d', t);
            Eyelink('Message', 'STIM_ONSET');
        end
    end
    if ET == 1 && t == 1
        Eyelink('Message', 'TRIALID %d', t);
        Eyelink('Message', 'STIM_ONSET');
    end
    % end
    KbQueueFlush();
    n = n+1;
    if (n == length(ex.flipTimes1)+1) % for any other block , reset frame index when previous trial ends 
        n = 1;
    end
    if t > ex.nTrials
        break;
    end
  
%             thisCond = ex.longFormConds(n+2);
%             Eyelink('Message', char(sprintf('Cond %s',
%             conditions(thisCond).name{:}))) 

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
