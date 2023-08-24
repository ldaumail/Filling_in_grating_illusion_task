function restingState_horiOffset(restingLength, horiOffset, vertOffset, debug)
    
%%%% resolution
if debug == 1
	experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,0); % laptop 2880, 1800 ,0
else
    experiment.resolution = SetResolution(max(Screen('Screens')),1024,768,60); % scanner 1024,768,60
end

% params
params.screenWidth = 17;                    % in cm; laptop=27.5, office=43, 3Tb=19, 7T=17, miniHelm=39;
params.viewingDist = 48;                    % in cm; 3Tb/office=43, 7T=48, miniHelm=57;
params.backgroundColor = [0 0 0];
params.textColor = [255 255 255];
params.fontSize = 20;

%%%% open screen
AssertOpenGL;
Screen('Preference', 'SkipSyncTests', 0);
screens = Screen('Screens');        % determines screen numbers. starts at 0 with main screen
screenNum = max(screens);           % chooses last screen to show stimuli.
[w, rect] = Screen('OpenWindow',screenNum, params.backgroundColor);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
load 7T_Sam.mat;
Screen('LoadNormalizedGammaTable', w, linearizedCLUT);

%%%% layout
params.vertOffset = vertOffset;             % vertical offset
params.horiOffset = horiOffset;  
params.gammaCorrect = 1;                    % make sure this is on at the scanner!
params.ppd = pi * rect(3) / (atan(params.screenWidth/params.viewingDist/2)) / 360; 
params.fixSizeDeg = .5;                     % fixation point size
params.littleFixDeg = params.fixSizeDeg *.5;% proportion of the fixSizeDeg occupied by the smaller black dot
params.outerFixPixels = 2;                  % in pixels, the black ring around fixation
params.fixSize = round(params.fixSizeDeg*params.ppd);
params.littleFix = round(params.littleFixDeg*params.ppd);

HideCursor;

%%%% center location 
xc = rect(3)/2+ params.horiOffset; 
yc = rect(4)/2 + params.vertOffset;
[width,height] = RectSize(Screen('TextBounds',w,'Experiment is starting...'));
Screen(w, 'DrawText', ['Experiment is starting...'], xc-width/2, yc-height/2, params.textColor);
Screen(w, 'Flip', 0);

% wait for trigger
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
deviceNumber = keyboardIndices(1);
KbTriggerWait(53, deviceNumber);

% make screen and fixation
Screen('FillRect', w, params.backgroundColor);
Screen('FillOval', w,[255 255 255], [xc-round(params.fixSize/2+params.outerFixPixels) yc-round(params.fixSize/2+params.outerFixPixels) xc+round(params.fixSize/2+params.outerFixPixels) yc+round(params.fixSize/2+params.outerFixPixels)]); % black fixation ring
Screen('FillOval', w,[0 0 0], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
Screen('FillOval', w,[255 255 255], [xc-round(params.littleFix/2) yc-round(params.littleFix/2) xc+round(params.littleFix/2) yc+round(params.littleFix/2)]); % black fixation dot   
Screen(w, 'Flip', 0);

% wait for duration of scan
WaitSecs(restingLength);

% finish
Screen('Close');
Screen('CloseAll');
ShowCursor;
