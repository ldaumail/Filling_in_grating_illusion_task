function contrast_matching_control_v3(subject, session, debug) 

%%contrast matching control
%Loic 01252023
%In this version, we add multiple velocities
% subject = 'Dave';                                                                                                                                                                                                                                                     
% session = 1;                                                                                                                           
% debug = 1;

ex.version = 'v13';
global EyeData rect w xc yc %eye_used
%%%% resolution 
if debug == 1
    % eyetracking on (1) or off (0)

    ex.screenWidth = 53.1;             % in cm; %laptop=27.5,office=43, %19=%T3b, miniHelm=39;
    ex.viewingDist = 53.5;             % in cm; 3Tb/office=43, miniHelm=57;
	ex.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop 1920,1080
    ex.gammaCorrection = 0;       % make sure this = 1 when you're at the scanner!
else
    
                                                                                                                         
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
ex.stim.contrast = 0.15;%linspace(0.01,0.20,10);%[0.05, 0.10, 0.15];                                                 % in %, maybe?? %here the number of stimulus contrast levels is the number of different conditions
ex.stim.contrastMultiplicator = ex.stim.contrast/2;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
ex.match.contrastMults = nan(10,length(ex.stim.contrast));
ex.match.contrastOffset = nan(10,length(ex.stim.contrast));
ex.stim.contrastOffset = .5;  % for procedural gabor, 0 = 0 mean luminance, 1 = 100% mean luminance 
for i =1:length(ex.stim.contrast)
    ex.match.contrastMults(:,i) = linspace(0,ex.stim.contrastMultiplicator(i)/4,10);%linspace(0,0.075,20);   %contrast multiplicators for the adjustable sine wave from 0 to 0.6 contrast
    ex.match.contrastOffset(:,i) = ex.stim.contrastOffset-ex.stim.contrastMultiplicator(i)-ex.match.contrastMults(:,i);
end                            
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

ex.conds = {'Cont1Vel1';...%'Cont2Vel1';'Cont3Vel1';'Cont4Vel1';'Cont5Vel1';'Cont6Vel1';'Cont7Vel1';'Cont8Vel1';'Cont9Vel1';'Cont10Vel1';...
    'Cont1Vel2'};%;'Cont2Vel2';'Cont3Vel2';'Cont4Vel2';'Cont5Vel2';'Cont6Vel2';'Cont7Vel2';'Cont8Vel2';'Cont9Vel2';'Cont10Vel2'};%'MinbgMatchVel2'};
ex.numConds = length(ex.conds);
% with line of code below we will have 1 condition per block, randomized. we might need to change that
% later, to have the conditions randomized within each block
ex.repsPerRun = 20;              % repetitions of each condition per run
ex.nTrials = ex.numConds*ex.repsPerRun;

ex.condShuffle = [];
for i =1:ex.repsPerRun
    ex.condShuffle = [ex.condShuffle, Shuffle([1:ex.numConds])];
end
% vel2Idx = find(ex.condShuffle >10);
% ex.condShuffle(vel2Idx) = ex.condShuffle(vel2Idx)-10;

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
      % [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
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
DrawFormattedText(w,'Wait while the task is loading'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
% rectLWaveIDO = nan(length(flipTimes),length(ex.stimDur),length(ex.stim.contrast));
% rectRWaveIDO = nan(length(flipTimes),length(ex.stimDur),length(ex.stim.contrast));
for s =1:length(ex.stim.cycPerSec)
    flipt =sprintf('flipTimes%d',s);
    flipTimes = ex.(flipt);
    
    %rect
    rectLWave = sprintf('rectLWave%d',s);
    rectRWave = sprintf('rectRWave%d',s);
    
    ex.(rectLWave) = nan(ex.gaborHeight,ex.gaborWidth*2,length(flipTimes),length(ex.stim.contrast));
    ex.(rectRWave) = nan(ex.gaborHeight,ex.gaborWidth*2,length(flipTimes),length(ex.stim.contrast));
    for c =1:length(ex.stim.contrast)
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
            ex.(rectLWave)(:,:,f,c) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                ex.stim.orientation,LPhase,ex.stim.contrastOffset,ex.stim.contrastMultiplicator(c),...
                ex.ppd);
            ex.(rectRWave)(:,:,f,c) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                ex.stim.orientation,RPhase,ex.stim.contrastOffset,ex.stim.contrastMultiplicator(c),...
                ex.ppd);
            %                 figure();
            %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
            %
            tmprectLWaveID(f) = Screen('MakeTexture', w, squeeze(ex.(rectLWave)(:,:,f,c)));
            tmprectRWaveID(f) = Screen('MakeTexture', w, squeeze(ex.(rectRWave)(:,:,f,c)));
            
        end        
        ex.rectLWaveID(1:length(tmprectLWaveID(:)),s,c) = tmprectLWaveID(:);
        ex.rectRWaveID(1:length(tmprectLWaveID(:)),s,c) =  tmprectLWaveID(:);

    end
end

%% Create contrast adjustable sine wave grating
DrawFormattedText(w,'The task is still loading'... % :  '...
    ,'center', 'center',[0 0 0]);
Screen(w, 'Flip', 0);
clear c
%rectCWaveIDO = nan(length(flipTimes),length(ex.stimDur),length(ex.stim.contrast)) ;
for s =1:length(ex.stim.cycPerSec)
    flipt =sprintf('flipTimes%d',s);
    flipTimes = ex.(flipt);
    
    %rect
    rectCWave = sprintf('rectCWave%d',s);
    
    ex.(rectCWave) = nan(ex.gaborHeight,ex.gaborWidth*2,length(flipTimes),size(ex.match.contrastMults,1),length(ex.stim.contrast));
    
    phaseNum = sprintf('phases%d',s);
    phases = ex.stim.(phaseNum) ;
    for c = 1:length(ex.stim.contrast)
        for m =1:size(ex.match.contrastMults,1)
            
            
            
            %       clear tmprectLWaveID tmprectRWaveID
            for f = 1:length(flipTimes)
                
                CPhase = phases(f);
                %ih in pixels %iw in pixels %spatial freq in cycles per dva
                %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
                %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
                %background color (unused if the grating is not an annulus)
                
                % rect
                ex.(rectCWave)(:,:,f,m,c) = makeSineGrating(ex.gaborHeight,ex.gaborWidth*2,ex.stim.spatialFreqDeg,...
                    ex.stim.orientation,CPhase,ex.match.contrastOffset(m,c) ,ex.match.contrastMults(m,c),...
                    ex.ppd);%-min(min(squeeze(ex.rectLWave1(1,:,:)),[],1));  %center contrast level to dark luminance background so that the max luminance matches the dark luminance background
                %                 figure();
                %                imshow(squeeze(ex.(rectRWave)(f,:,:,o))./max(squeeze(ex.(rectRWave)(f,:,:,o)),[],'all'));
                %
                tmprectCWaveID(f,m) = Screen('MakeTexture', w, squeeze(ex.(rectCWave)(:,:,f,m,c)));
                
            end
            
            %% save
            ex.rectCWaveID(1:length(tmprectCWaveID(:,m)),s,m,c) = tmprectCWaveID(:,m);
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


%% %%%% initial window - instructions and wait for space

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

gray = repmat(mean(squeeze(ex.rectLWave1(1,:,1,1))), [1,3]);

% ex.fixCol1Grad = linspace(255,gray(1),90);% make red dot disapear in 90 flips = 1.5 sec ; logspace(log10(255),log10(gray(1)));
% ex.fixCol2Grad = linspace(0,gray(1),90);

Screen('FillRect', w, gray);
    
% if ET
%     gcnt = 0; 
%     EyeData.mx=nan(1,1);
%     EyeData.my=nan(1,1);
%     EyeData.ma=nan(1,1);
%     EyeData.FixDoneT = nan(1,1);
%     EyeData.gazeD = nan(1,1);
%     EyeData.Fixated = nan(1,1);
% end

n=1; %initialize at first flip
contM = size(ex.match.contrastMults,1)/2; %initial contrast level index for adjustable stimulus
ex.contM = nan(ex.nTrials,length(ex.stim.contrast));
% tstartcnt = 0;
% tend = [];
t =1; %trial number
while(1) %n+1 < length(ex.allFlips)
    KbQueueStart();
%     if ET
%         run GetEyeDataLoic; %check eyetracker
%     end
    %%%% draw sine wave grating stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thisCond = ex.condShuffle(t);
    %screen background color
    speed = ex.conds{thisCond};
    speed = speed(end);
    s = str2double(speed);
    rectLWave = sprintf('rectLWave%s', speed);
    if thisCond >length(ex.conds)/2 % we need to do this because with n conditions, n/2 are at velocity 1 and n/2 at velocity 2
        condInd = thisCond -length(ex.conds)/2;
    else 
        condInd = thisCond;
    end
    gray = repmat(min(min(squeeze(ex.(rectLWave)(:,:,1,condInd)),[],1)), [1,3]);
    Screen('FillRect', w, gray);
    
    ex.rectLRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xL,yL);
    ex.rectRRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xR,yR);
    ex.rectCRect =  CenterRectOnPoint([0 0 ex.rectGaborWidth ex.gaborHeight],xC,yC);
    
    if nnz(find(ex.rectLWaveID(n,s,condInd)))
        % top stim
        Screen('DrawTexture', w, ex.rectLWaveID(n,s,condInd),[],ex.rectLRect);
        % bottom stim
        Screen('DrawTexture', w, ex.rectRWaveID(n,s,condInd),[], ex.rectRRect);
        % side stim
        Screen('DrawTexture', w, ex.rectCWaveID(n,s,contM,condInd), [], ex.rectCRect);
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
            contM = contM+1; % increase contrast
        elseif (find(firstPress,1) == KbName('2')|find(firstPress,1) == KbName('2@'))
            contM = contM-1; % decrease contrast
        end
        if contM == 0 || contM == length(ex.match.contrastMults)+1
            contM = size(ex.match.contrastMults,1)/2;
        end
    elseif (pressed && ismember(find(firstPress,1), [KbName('Return') KbName('ENTER')]))
        ex.contM(t,thisCond) = ex.match.contrastMults(contM,condInd);
        n = 1;
        t = t+1;
        contM = size(ex.match.contrastMults,1)/2; 
%         if  ET == 1%isempty(find(mod(tstartcnt,2)))
%             Eyelink('Message', 'STIM_OFFSET');
%             Eyelink('Message', 'TRIALID %d', t);
%             Eyelink('Message', 'STIM_ONSET');
%         end
    end
%     if ET == 1 && t == 1
%         Eyelink('Message', 'TRIALID %d', t);
%         Eyelink('Message', 'STIM_ONSET');
%     end
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
ex.rectLWave1 = [];
ex.rectRWave1 = [];
ex.rectCWave1 = [];
ex.rectLWave2 = [];
ex.rectRWave2 = [];
ex.rectCWave2 = [];

savedir = fullfile(ex.root,'data',sprintf('s%s_v3/',subject));
if ~exist(savedir); mkdir(savedir); end
savename = fullfile(savedir, strcat(sprintf('/s%s_contrast_matching_v3_date%s',subject,num2str(ex.date)), '.mat'));
save(savename,'ex','-v7.3')

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
fclose all;                                                                                                                           


% if ET
%     Eyelink('StopRecording');   
%     Eyelink('CloseFile');
%         %download edf file
% %     if downET
%         Eyelink('Command', 'set_idle_mode');
%     try
%         fprintf('Receiving data file ''%s''\n',  EyeData.edfFile);
%         status = Eyelink('ReceiveFile');
%         
%         if status > 0
%             fprintf('ReceiveFile status %d\n', status);
%         end
%         if 2==exist(EyeData.edfFile, 'file')
%             fprintf('Data file ''%s'' can be found in ''%s''\n',  EyeData.edfFile, pwd );
%         end
%     catch
%         fprintf('Problem receiving data file ''%s''\n',  EyeData.edfFile);
%     end
% %     end
%     
%     Eyelink('Shutdown');
%     savename = fullfile(savedir, strcat(sprintf('/s%s_smooth_pursuit_v13_date%s_fix_eyeDat',subject,num2str(ex.date)), '.mat'));
%     save(savename, 'EyeData')
% end
