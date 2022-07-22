% LD 7.20.2022
% present two rectangular orientation patches on left side of fixation, with a
% uniform field surrounding them
% RSVP task at fixation throughout
% block design


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
params.initialFixation = 1;%16;     % in seconds
params.finalFixation = 16;          % in seconds
params.phaseFlicker = 1; %.2;       % in seconds (on for Xs, off for Xs, phase changes

%%%% gabor properties
params.stim.spatialFreqDeg = 1.5;                                           % cycles per degree of visual angle
params.stim.contrast =  .3;                                                 % in %, maybe??
params.stim.orientation = 90;                                                % in degrees
%params.stim.guassianSpaceConstant = .4;                                     % approx equal to the number of radians covered by one standard deviation of the radius of the gaussian mask.
params.stim.fromFixation = .6;                                              % in degrees of visual angle
params.stim.gaborHDeg = 4;                                                  % in degrees of visual angle
params.stim.gaborWDeg = 6; 
params.stim.ringPix = 3;                                                    % in pixels, thickness of greyscale ring separating
params.stim.contrastMultiplicator = .2;                                     % for sine wave 0.5 = 100% contrast, 0.2 = 40%
params.stim.contrastOffset = [.5 .5 .5 0];                                  % for procedural gabor
params.stim.motionRate = 1.3*360 ;                                          % 1.3 cycles per second = 360 deg of phase *1.3 per sec

%%%% conditions & layout
params.numMotions = 2;
params.motions = {'left','right'}; %back and forth starts from the left side, back and forth starts from the right side
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

% HideCursor;
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
flipTimes = [0:flipInt*5:params.phaseFlicker]; %multiply flipInt by 5 to flip the image every 5 frames 
flipTimes = flipTimes(1:length(flipTimes)-1);
params.stim.motionPerFlip = params.stim.motionRate * flipInt*5; %degrees per flip, here we multiply by five to move the phase 5 times more after each flip


%%%% scale the stims for the screen
params.ppd = pi* rect(3) / (atan(params.screenWidth/params.viewingDist/2)) / 360; %2pi*(rect(3)/2)= pi*rect(3)
%params.freq =  (params.stim.spatialFreqDeg)*2*pi/params.ppd;                  %converts cycles per degree to cycles *(rad)/ pixel
params.gaborHeight = round(params.stim.gaborHDeg*params.ppd);                 % in pixels, the size of our objects
params.gaborWidth = round(params.stim.gaborWDeg*params.ppd);                 % in pixels, the size of our objects
params.fromFix = round(params.stim.fromFixation*params.ppd);
params.fixSize = round(params.fixSizeDeg*params.ppd);
params.littleFix = round(params.littleFixDeg*params.ppd);


%%% create sine wave gratings and store all phase transitions in structure
%%% along with pointers
topPhase = randi(360);
bottomPhase = topPhase;
params.topWave = nan(length(flipTimes),params.gaborHeight,params.gaborWidth);
params.bottomWave = nan(length(flipTimes),params.gaborHeight,params.gaborWidth);
params.topWaveID = nan(length(flipTimes),1);
params.bottomWaveID = nan(length(flipTimes),1);

for f = 1:length(flipTimes)
    
    if f <= length(flipTimes)/2
        
    topPhase = mod(topPhase + params.stim.motionPerFlip,360);
    bottomPhase = mod(bottomPhase + params.stim.motionPerFlip,360);
    %ih in pixels %iw in pixels %spatial freq in cycles per dva
    %tilt/orientation of the grating in degrees %phase in degrees (not degrees of visual angle)
    %contrast offset in percent    %contrast multiplicator  %ppd = 0 if freq already in cycles per stimulus
    %background color (unused if the grating is not an annulus)
    
    params.topWave(f,:,:) = makeSineGrating(params.gaborHeight,params.gaborWidth,params.stim.spatialFreqDeg,...
        params.stim.orientation,topPhase,params.stim.contrastOffset(1),params.stim.contrastMultiplicator,...
        params.ppd);
    params.bottomWave(f,:,:) = makeSineGrating(params.gaborHeight,params.gaborWidth,params.stim.spatialFreqDeg,...
        params.stim.orientation,bottomPhase,params.stim.contrastOffset(1),params.stim.contrastMultiplicator,...
        params.ppd);
%    figure();
 %   imshow(squeeze(topWave(f,:,:)));
%     
    elseif f > length(flipTimes)/2
        topPhase = mod(topPhase - params.stim.motionPerFlip,360);
        bottomPhase = mod(bottomPhase - params.stim.motionPerFlip,360);   
        params.topWave(f,:,:) = makeSineGrating(params.gaborHeight,params.gaborWidth,params.stim.spatialFreqDeg,...
            params.stim.orientation,topPhase,params.stim.contrastOffset(1),params.stim.contrastMultiplicator,...
            params.ppd);
        params.bottomWave(f,:,:) = makeSineGrating(params.gaborHeight,params.gaborWidth,params.stim.spatialFreqDeg,...
            params.stim.orientation,bottomPhase,params.stim.contrastOffset(1),params.stim.contrastMultiplicator,...
            params.ppd);
    end
    params.topWaveID(f) = Screen('MakeTexture', win, squeeze(params.topWave(f,:,:)));
    params.bottomWaveID(f) = Screen('MakeTexture', win, squeeze(params.bottomWave(f,:,:)));
            
end

%% Sine wave gratings locations
xc = rect(3)/2; % rect and center, with the flixibility to resize & shift center - change vars to zero if not used.
yc = rect(4)/2; %+params.vertOffset;

xtop = rect(3)/4;
ytop = rect(4)/4;
params.topRect =  CenterRectOnPoint([0 0 params.gaborWidth params.gaborHeight],xtop,ytop);

xbottom = rect(3)/4;
ybottom = rect(4)/4*3;
params.bottomRect =  CenterRectOnPoint([0 0 params.gaborWidth params.gaborHeight],xbottom,ybottom);


%%%% initial window - wait for backtick
Screen(win, 'DrawText', 'Waiting for Backtick.', 10,10,[0 0 0]);
Screen(win, 'Flip', 0);

%%% for the loc
save LGNsurroundParams.mat params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         experiment                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% start recording the response

[touch,tpress, keyCode]= KbCheck;
% experiment.numCorrect = 0;
% experiment.correctTarget = [];
flipCount = 1;

%%%%%%% START task TASK/FLIPPING
for n = 1:(length(experiment.allFlips)-1)
    thisCond = experiment.longFormConds(n);
    
    %%%% draw gabors
    if thisCond > 0 && experiment.longFormFlicker(n) > 0% zeros correspond to blanks, in which case we skip this next section
        
%         if experiment.longFormConds(n-1)==0 % if this is the beginning of a block, we need to randomize the phases
%             topPhase = randi(360);
%             if conditions(thisCond).phase == 2 % 1 = oppPhase, 2 = inPhase
%                 bottomPhase = topPhase; 
%             else
%                 bottomPhase = mod(topPhase - 180,360); % opp phase for oppPhase
%             end
%         end
        
        
        % draw & increment stims
        stimFlipCnt = 0;
        for motionFlip = flipTimes
            stimFlipCnt = stimFlipCnt+ 1; 
            % top stim
             Screen('DrawTexture', win, params.topWaveID(stimFlipCnt),[],params.topRect);
             
             % bottom stim
             if strcmp(conditions(thisCond).name, 'double-inPhase') || strcmp(conditions(thisCond).name, 'double-oppPhase')  % draw second stim if it is 'double-inPhase' or 'double-oppPhase' %if it's randomized or incognruent
                 Screen('DrawTexture', win, params.bottomWaveID(stimFlipCnt), [], params.bottomRect);
             end
             

            % draw fixation and RSVP letter
            Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
%             [width,height] = RectSize(Screen('TextBounds',win,experiment.letterSequence{n}));
%             Screen(win, 'DrawText', experiment.letterSequence{n}, xc-width/2, yc-height/2,params.cueColor);
%             
            %%%%%%%%%%% FLIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [VBLT experiment.flipTime(flipCount) FlipT missed] = Screen(win, 'Flip', experiment.startRun + experiment.allFlips(n)+motionFlip - slack);
            flipCount = flipCount+1;
            
            %%%% increment phase to show motion
%             if strfind(conditions(thisCond).name{:},'double') ==1
%                 if strcmp(conditions(thisCond).stimTop,'left') == 1
%                     topPhase = mod(topPhase + params.stim.motionPerFlip,360);
%                     bottomPhase = mod(bottomPhase - params.stim.motionPerFlip,360);
%                 elseif strcmp(conditions(thisCond).stimTop,'right') == 1
%                     topPhase = mod(topPhase - params.stim.motionPerFlip,360);
%                     bottomPhase = mod(bottomPhase + params.stim.motionPerFlip,360);
%                 end
%             elseif strfind(conditions(thisCond).name{:},'single') ==1 % single stimulus condition
%                 if strcmp(conditions(thisCond).stimTop,'left') == 1
%                     topPhase = mod(topPhase + params.stim.motionPerFlip,360);
%                    % bottomPhase = mod(bottomPhase + params.stim.motionPerFlip,360);
%                 elseif strcmp(conditions(thisCond).stimTop,'right') == 1
%                     topPhase = mod(topPhase - params.stim.motionPerFlip,360);
%                    % bottomPhase = mod(bottomPhase - params.stim.motionPerFlip,360);
%                 end
%             end
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


%KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');
%fprintf('Hit rate this run: %.2f%%\n',100*experiment.performance)
fclose all;
clear all;

%end


