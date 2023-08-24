function kayPRF_dotTask(subject, session, horiOffset, vertOffset, debug, stim)

% function kayPRF_dotTask(subject, subjectNum, runNum, deviceNumber, vertOffset)
% in this version (mar 30 2016): mulitbar or wedgering is flexibly set
% participant detects red dot at center fixation

% multibar indices are the indices of masks to be shown in this display, at
% 15hz; if the index isn't 0, it points to the correct mask in masks.mat
% from kendrickkay.net: The bar sweeps occur in 32-second chunks, each of
% which consists of a 28-second sweep followed by a 4-second rest. The
% intended frame rate is 15 Hz; thus, the total run duration is 4500/15 =
% 300 seconds = 5 minutes. The overall temporal structure is [16-second
% rest, L-R, D-U, R-L, U-D, 12-second rest, LL-UR, LR-UL, UR-LL, UL-LR,
% 16-second rest] = [16+32*4+12+32*4+16=300].

% 'wedgeringindices' is a 1 x 4500 vector with mask indices, similar to
% 'multibarindices'. This "wedgering" stimulus consists of rotating wedges
% and expanding and contracting rings. The wedges revolve in 32 seconds;
% the rings sweep for 28 seconds and are followed by 4 seconds of rest. The
% total run duration is 4500/15 = 300 seconds = 5 minutes. The overall
% temporal structure is [22-second rest, CCW, CCW, expand, expand, CW, CW,
% contract, contract, 22-second rest] = [22+32*8+22=300].
%
% this version: plays kay movie stimuli rather than flickering
% checkerboards. timing is the same

% clear all;
% input(['Running the ' run.masks ' version. Hit enter to proceed.']);

%%%%%%%%%%%%%%%%%%
% basic setting
%%%%%%%%%%%%%%%%%%

Screen('Preference', 'SkipSyncTests', 0);

%%%% resolution
if debug == 1
    experiment.resolution = SetResolution(max(Screen('Screens')),2880, 1800 ,0); % laptop 2880, 1800 ,0  1024,768,0
else
    experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % scanner1024,768,60
end

%%%% keyboards
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);

%%%% basic naming set-up
run.subject = subject;
run.session = session;       % to keep both of these structs labeled
run.vertOffset = vertOffset;
run.horiOffset = horiOffset; 
run.masks = stim;     % 'Multibar' or 'WedgeRing';
run.scanNum = input('Scan number :');
run.runNum = input('Run number :');

%%%% kay stimuli parameters, hardcoded in kayMasks.mat
eval(['load kay' run.masks '.mat']); 
load kayStims.mat;
if strcmp(run.masks,'WedgeRing')==1
    kay.wedgeLength = 32;   % in seconds
    kay.ringLength = 28;
    kay.postRingOff = 4;    % after each ring sweep only
    kay.initialFixation = 22;
    kay.finalFixation = 22;
    kay.totalSeconds = 300;
    kay.frameRate = 1/15;
    maskInds = wedgeringindices;
elseif strcmp(run.masks,'Multibar')==1
    kay.sweepLength = 28;   % in seconds
    kay.offLength = 4;      % after each sweep
    kay.initialFixation = 16;
    kay.interimOff = 12;
    kay.finalFixation = 16;
    kay.totalSeconds = 300;
    kay.frameRate = 1/15;
    maskInds = multibarindices;
end

%%%% our stimuli parameters
params.degSize = 14;                        % size of our circular stimulus
params.gammaCorrect = 1;                    % make sure this is on at the scanner!
params.fixSizeDeg =  .5;                    % in degrees, the size of the biggest white dot in the fixation
params.littleFixDeg = params.fixSizeDeg* .5;% proportion of the fixSizeDeg occupied by the smaller black dot
params.outerFixPixels = 2;                  % in pixels, the black ring around fixation
params.backgroundColor = [126 126 126];     % color
params.screenWidth = 17;                    % in cm; %laptop=27.5,office=43, %19=%T3b, 7T=17, miniHelm=39;
params.viewingDist = 48;                    % in cm; 3Tb/office=43, 7T=48, miniHelm=57;

%%%% task parameters - detect color change at fixation
task.duration = kay.frameRate *3;           % in s
task.endBuff = 4;                           % in s
task.listen = kay.frameRate *4;             % in s
task.color = [255 0 0];                     % dot change color
task.prob = .1;                             % probabiltiy that each task + buffer interval will have a target

%%% setup task
numTargets = round(kay.totalSeconds/task.duration * task.prob);
task.targetInd = zeros(1,kay.totalSeconds/task.duration);
while 1
    maybeTarget= task.endBuff/task.duration+Randi(length(task.targetInd)-2*task.endBuff/task.duration);
    if sum(task.targetInd(maybeTarget-2:maybeTarget-1)) == 0 && sum(task.targetInd(maybeTarget+1:maybeTarget+2)) == 0
        task.targetInd(maybeTarget) = 1;
    end
    if sum(task.targetInd) == numTargets
        break
    end
end
longFormInd = Expand(task.targetInd,task.duration/kay.frameRate,1);

%%%% set-up rand
% rand('twister', sum(100*clock));
% run.rand = rand;
ex.startTime = clock;
rng(sum(100*ex.startTime));
ex.rand = rng;
%%%% files and things
run.root = pwd; 
run.date = datestr(now,30);

%%%% timing
allFlips = 0:kay.frameRate:(kay.totalSeconds-kay.frameRate);
numTex = size(patterns,4); allTexNum = [randi(numTex) zeros(1,length(allFlips-1))];
for n = 2:length(allFlips)
    while 1
        maybeTex = randi(numTex);
        if maybeTex ~= allTexNum(n-1)
            allTexNum(n) = maybeTex;
            break
        end
    end
end

%%%%%%%%%%%%%%%%
% open screen
%%%%%%%%%%%%%%%%

HideCursor;
Priority(9);

%%%% open screen
screen=max(Screen('Screens'));
[win, rect]=Screen('OpenWindow',screen,180,[],[],[],[],[],kPsychNeed32BPCFloat);
Screen(win, 'TextSize', 20);

%%%% center location
xc = rect(3)/2+horiOffset; 
yc = rect(4)/2+vertOffset;

%%% gamma correction
    if params.gammaCorrect > 0
        load 7T_Sam.mat
        Screen('LoadNormalizedGammaTable', screen, linearizedCLUT);
    end

%%%% timing optimization
flipInt = Screen('GetFlipInterval',win);
slack = flipInt/2;
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%%%% initial window - wait for backtick
Screen('FillRect',win,params.backgroundColor)
Screen(win, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(win, 'Flip', 0);

%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2

%%% set up stims
params.ppd = pi* rect(3) / (atan(params.screenWidth/params.viewingDist/2)) / 360;
params.fixSize = round(params.fixSizeDeg * params.ppd);
params.littleFix = round(params.littleFixDeg * params.ppd);

for n = 1:size(patterns,4);
    patternTex{n} = Screen('MakeTexture',win,patterns(:,:,:,n));
end

%%%% set up bar mask generation
stimRadius = round(params.ppd*params.degSize/2);
maskRadius = round(size(masks,1)/2);
[x,y]=meshgrid(-maskRadius:maskRadius, -maskRadius:maskRadius);
barMask=uint8(ones(2*maskRadius, 2*maskRadius,2) * params.backgroundColor(1));

%%%% get trigger
KbTriggerWait(53, deviceNumber);
KbQueueCreate(deviceNumber,responseKeys);
presses = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PRF MAPPING O'CLOCK                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% start recording the response
KbQueueStart();
task.numCorrect = 0;
task.correctTarget = [];

%%%% START LOC LIPPING
for n = 1:(length(allFlips))
    
    % if we're not in a blank
    if maskInds(n)>0
        
        %%%% checkerboards
        Screen('DrawTexture', win, patternTex{allTexNum(n)}, [], [xc-stimRadius yc-stimRadius xc+stimRadius yc+stimRadius]);
        
        %%%% mask
        barMask(:,:,2)=imcomplement(uint8(masks(:,:,maskInds(n))));
        
        %%%% Build a single transparency mask texture
        maskTex=Screen('MakeTexture', win, barMask);
        Screen('DrawTexture', win, maskTex, [], [xc-stimRadius yc-stimRadius xc+stimRadius yc+stimRadius]);
    end
    
    %%%% draw fixation
    Screen('FillOval', win,[0 0 0], [xc-round(params.fixSize/2+params.outerFixPixels ) yc-round(params.fixSize/2+params.outerFixPixels ) xc+round(params.fixSize/2+params.outerFixPixels ) yc+round(params.fixSize/2+params.outerFixPixels )]); % black fixation ring
    Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
    
    % task - is this dot black or colored?
    if longFormInd(n)==0
    Screen('FillOval', win,[0 0 0], [xc-round(params.littleFix/2) yc-round(params.littleFix/2) xc+round(params.littleFix/2) yc+round(params.littleFix/2)]); % black fixation dot
    else
    Screen('FillOval', win,task.color, [xc-round(params.littleFix/2) yc-round(params.littleFix/2) xc+round(params.littleFix/2) yc+round(params.littleFix/2)]); % black fixation dot
    end
    
    %%%% FLIP 
    if n == 1 [VBLT run.startRun FlipT missed] = Screen(win, 'Flip', 0);
    else [VBLT run.flipTime(n) FlipT missed] = Screen(win, 'Flip', run.startRun + allFlips(n) - slack); end
    
    % listen for response  - correct if you respond to previous 3 letters
    [pressed, firstPress]= KbQueueCheck();
    if n > task.endBuff/kay.frameRate % if the task has started
        if pressed == 1 presses = [presses n]; end
        targetRange = longFormInd(n - task.listen/kay.frameRate-1:n-1);
        if (pressed ==1) && (sum(targetRange)>0)
            task.numCorrect = task.numCorrect + 1;
            task.correctTarget = [task.correctTarget n];
        end % then it's correct
        KbQueueFlush();
    end
    if maskInds(n)>0 Screen('Close', maskTex); end
end

%%%% to show the very last flip screen for its 6ms
[VBLT run.flipTime(n+1) FlipT missed] = Screen(win, 'Flip', run.startRun + allFlips(end) - slack);

%%%%%%%%%%%%%%%%%%
% done! wrap up  
%%%%%%%%%%%%%%%%%%

run.runTime = GetSecs - run.startRun;
task.perf = task.numCorrect/numTargets;
fprintf('Hit rate this run: %.2f%%\n',100*task.perf);

savedir = fullfile(pwd,'data',subject,session,'kayPRF');
if ~exist(savedir); mkdir(savedir); end
save(fullfile(savedir, strcat(subject, '_kayPRF_sn', num2str(run.scanNum), '_rn', num2str(run.runNum), '_', run.date, '.mat')), 'run', 'params', 'kay', 'task');

ShowCursor;
Screen('Close');
Screen('CloseAll');

clear all;
end

