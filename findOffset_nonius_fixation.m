function [vertOffsets,horiOffsets] = findOffset_nonius_fixation(subject, session, debug)

%%%%%%%%%%%%%%%%%%%
% BASIC SETTINGS
%%%%%%%%%%%%%%%%%%%

%%%% resolution
if debug == 1
	SetResolution(max(Screen('Screens')), 1280,1024,85);
else
    SetResolution(max(Screen('Screens')),1280,1024,85); 
end
experiment.runNum = input('Run number :');

%%%% keyboards
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);

%%%% input this at the beginning of the scan session for 7T
ex.vertOffsetL = 0;                      % vertical offset 
ex.horiOffsetL = 0;
ex.vertOffsetR = 0;                      % vertical offset 
ex.horiOffsetR = 0;

% scSize.screenDeg = 14;                      % pRF screen size
ex.screenDegY = 8; %figSizeDeg;
ex.screenDegX = 16; %figSizeDeg*3;%9;                      % pRF screen size
% scSize.centerDeg = 4;                      % center circle size in diameter
% scSize.surroundDeg = 12;                      % surround circle size in diameter

%%%% scales all of the stimuli in DVA to the screensize
ex.screenWidth = 40;                    % in cm; laptop=27.5, office=43, 3Tb=19, 7T=17, miniHelm=39;
ex.viewingDist = 46;                    % in cm; 3Tb/office=43, 7T=48, miniHelm=57;
ex.backgroundColor = [128 128 128];     % color
ex.boxColor = ex.backgroundColor * .7;

%%%% Fixation
ex.littleFixDeg =  .2;            % in degrees, the size of the biggest white dot in the fixation
ex.bigFixSizeDeg = 0.5;
ex.outerFixPixels = 2;          % in pixels, the black ring around fixation

%%%% Nonius lines 
ex.lineHdeg = 0.4;
ex.lineWdeg = 0.06;

%%% horizontal line 
ex.horiLineWdeg = 0.7;

%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('1'))=1; % button box 1
responseKeys(KbName('2'))=1; % button box 2
responseKeys(KbName('3'))=1; % button box 3
responseKeys(KbName('4'))=1; % button box 3
responseKeys(KbName('Enter'))=1; % button box 3
%%%%%%%%%%%%%%%%
% OPEN SCREEN 
%%%%%%%%%%%%%%%%

HideCursor;
Priority(9); % highest priority 

%%%% open screen
Screen('Preference', 'SkipSyncTests', 1);
screen=max(Screen('Screens'));
[w, rect]=Screen('OpenWindow',screen,180,[],[],[],[],[],kPsychNeed32BPCFloat);
Screen(w, 'TextSize', 32);


% xc = rect(3)/2; 
% yc = rect(4)/2;

%Load the gamma table of the specified screen
load("phase2_photometry.mat");
Screen('LoadNormalizedGammaTable', w, inverseCLUT);
%%%% scale the stims for the screen
ex.ppd = pi* rect(3) / (atan(ex.screenWidth/ex.viewingDist/2)) / 360;
ex.fixSize = round(ex.littleFixDeg*ex.ppd);
ex.screenX = ex.screenDegX*ex.ppd;
ex.screenY = ex.screenDegY*ex.ppd;

%%%% Fixation dot
% ex.yFixOffset = 3*ex.ppd;

KbQueueCreate(deviceNumber,responseKeys);

%%%% find center Y
Screen('FillRect',w,ex.backgroundColor)
text = 'Find CENTER Y. \n\n 1 = Left Up, 2 = Left Down, \n\n 3 = Right Up, 4: Right Down, \n\n Enter: Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',w,text));
DrawFormattedText(w, text, 'center', 'center');
Screen(w, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

%%%%%%%%
% RUN 
%%%%%%%%

while 1
    
%     xS = rect(3)/2+ex.horiOffset;
%     yS = rect(4)/2+ex.vertOffset;
    
    xcL = rect(3)/2+ex.horiOffsetL;
    xcR = rect(3)/2+ex.horiOffsetR;
    ycL = rect(4)/2+ex.vertOffsetL; 
    ycR = rect(4)/2+ex.vertOffsetR; 
    
    KbQueueStart();    
    
    Screen('DrawDots', w, [xcR*3/2 ycR], ex.fixSize, [255 255 255], [], 2); %Right dot
    
%     Screen('DrawLines', w, [xlineLl, xlineLr; ylineLBot, ylineLTop], ex.lineW, [255 255 255]);
%     Screen('DrawLines', w, [xlineRl, xlineRr; ylineRBot, ylineRTop], ex.lineW, [255 255 255]);
%     Screen('DrawLines', w, [xhorilineLl, xhorilineLr; yhoriline, yhoriline],ex.lineW, [255 255 255]);
%     Screen('DrawLines', w, [xhorilineRl, xhorilineRr; yhoriline, yhoriline], ex.lineW, [255 255 255]);

    %%%% FLIP 
    Screen(w, 'Flip', 0);
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if find(firstPress,1) == KbName('1')
        ex.vertOffsetL = ex.vertOffsetL - 5;
        KbQueueFlush();
        Screen('DrawDots', w, [xcL/2 ycL], ex.fixSize, [255 255 255], [], 2); %Left dot
    elseif find(firstPress,1) == KbName('2')
        ex.vertOffsetL = ex.vertOffsetL + 5;
        KbQueueFlush();
        Screen('DrawDots', w, [xcL/2 ycL], ex.fixSize, [255 255 255], [], 2); %Left dot
    elseif find(firstPress,1) == KbName('3')
        ex.vertOffsetR = ex.vertOffsetR - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('4')
        ex.vertOffsetR = ex.vertOffsetR + 5;
        KbQueueFlush();
    
    elseif find(firstPress,1) == KbName('Enter')
        break;
        KbQueueFlush();
    end
end

vertOffsets = [ex.vertOffsetL; ex.vertOffsetR];

Screen('FillRect',w,ex.backgroundColor)
text = 'Find CENTER X. \n\n 1 = Left left, 2 = Left right, \n\n 3 = Right left, 4: Right right,\n\n Enter: Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',w,text));
DrawFormattedText(w, text, 'center', 'center');
Screen(w, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

while 1
    
    xcL = rect(3)/2+ex.horiOffsetL;
    ycL = rect(4)/2+ex.vertOffsetL;
    xcR = rect(3)/2+ex.horiOffsetR;
    ycR = rect(4)/2+ex.vertOffsetR;
    KbQueueStart();
    %%%% draw screen square

    Screen('DrawDots', w, [xcR*3/2 ycR], ex.fixSize, [255 255 255], [], 2);
    %%%% FLIP 
    Screen(w, 'Flip', 0);
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if find(firstPress,1) == KbName('1')
        ex.horiOffsetL = ex.horiOffsetL - 5;
        KbQueueFlush();
        Screen('DrawDots', w, [xcL/2 ycL], ex.fixSize, [255 255 255], [], 2);
    elseif find(firstPress,1) == KbName('2')
        ex.horiOffsetL = ex.horiOffsetL + 5;
        KbQueueFlush();
        Screen('DrawDots', w, [xcL/2 ycL], ex.fixSize, [255 255 255], [], 2);
    elseif find(firstPress,1) == KbName('3')
        ex.horiOffsetR = ex.horiOffsetR - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('4')
        ex.horiOffsetR = ex.horiOffsetR + 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('Enter')
        break;
        KbQueueFlush();
    end
end

horiOffsets = [ex.horiOffsetL; ex.horiOffsetR];

%%%%%%%%%%%%%%%%%%
% DONE! WRAP UP
%%%%%%%%%%%%%%%%%%

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');

% save
% savedir = fullfile(pwd,'data',subject,session,'findOffset_figureGround');
% if ~exist(savedir); mkdir(savedir); end
% savename=fullfile(savedir, strcat(subject,'_screenSize.mat'));   
% save(savename,'ex')
% 
% end

