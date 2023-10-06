function [vertOffsets,horiOffsets] = findOffset_nonius_fixationV3(subject, session, debug)

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
%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('2'))=1; % button box 1
responseKeys(KbName('4'))=1; % button box 2
responseKeys(KbName('6'))=1; % button box 3
responseKeys(KbName('8'))=1; % button box 4
responseKeys(KbName('LeftArrow'))=1; % button box 5
responseKeys(KbName('RightArrow'))=1; % button box 5
responseKeys(KbName('Enter'))=1; % button box 3

Screen('Preference', 'SkipSyncTests', 0);


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

%%%% Lines
ex.lineHdeg = 1/2;%1/sqrt(2);
ex.lineThdeg = 0.1;
ex.lineHdegSq = 5; %square size
%%% horizontal line
ex.horiLineWdeg = 0.7;


%%%%%%%%%%%%%%%%
% OPEN SCREEN
%%%%%%%%%%%%%%%%

HideCursor;
Priority(9); % highest priority

%%%% open screen
% Screen('Preference', 'SkipSyncTests', 1);
screen=max(Screen('Screens'));
[w, rect]=Screen('OpenWindow',screen,180,[],[],[],[],[],kPsychNeed32BPCFloat);
%     [w, rect]=Screen('OpenWindow',screen,ex.backgroundColor,[],[],[],[],[],kPsychNeed32BPCFloat); %might need to switch 900 and 600 by 1600 and 1200 for room 425
Screen(w, 'TextSize', 32);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % Set up alpha-blending for smooth (anti-aliased) lines


% xc = rect(3)/2;
% yc = rect(4)/2;

%Load the gamma table of the specified screen
load("phase2_photometry.mat");
Screen('LoadNormalizedGammaTable', w, inverseCLUT);
%%%% scale the stims for the screen
ex.ppd = pi* rect(3) / (atan(ex.screenWidth/ex.viewingDist/2)) / 360;
ex.fixSize = round(ex.littleFixDeg*ex.ppd);

%%%% Fixation dot
% ex.yFixOffset = 3*ex.ppd;
ex.lineH = ex.lineHdeg*ex.ppd;
ex.lineW = ex.lineHdeg*ex.ppd;
ex.lineTh = ex.lineThdeg*ex.ppd;
ex.lineHsq = ex.lineHdegSq*ex.ppd;
ex.lineColor = [0 0 0];


%%%% find center Y
Screen('FillRect',w,ex.backgroundColor)
text = 'Find CENTER X,Y for left eye. \n\n 8 =  Up, 2 = Down, \n\n 6 = Right , 4: Left, \n\n Enter: Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',w,text));
xc =rect(3)/2;
yc = rect(4)/2;
DrawFormattedText(w, text, xc/5, yc/2);
DrawFormattedText(w, text, xc+xc/5, yc/2);
Screen(w, 'Flip', 0);
% WaitSecs(1); KbPressWait; KbQueueFlush();
% KbTriggerWait(KbName('Space'), deviceNumber);
KbWait(deviceNumber,2)
text = 'Press left or right arrow \n\n to select a square to move.';   
DrawFormattedText(w, text, xc/5, yc/2);
DrawFormattedText(w, text, xc+xc/5, yc/2);
Screen(w, 'Flip', 0);
%%%%%%%%
% RUN
%%%%%%%%
KbQueueCreate(deviceNumber,responseKeys);
KbQueueStart();
KbWait(deviceNumber,2)

% % [secs, keyCode, ~] = KbWait(deviceNumber,2); %KbCheck; 

%  [pressed, firstPress]= KbQueueCheck();
% if find(firstPress,1) == KbName('LeftArrow')%
%     chosenSquare = 'Left';
%     %     break;
%     KbQueueFlush();
% elseif find(firstPress,1) == KbName('RightArrow')%keyCode(KbName('RightArrow'))
%     chosenSquare = 'Right';
%     %     break;
%     KbQueueFlush();
% end
% KbQueueStop();
chosenSquare = 'Left';
while 1
    [keyIsDown, secs, keyCode] = KbCheck;
    if keyIsDown
        if keyCode(KbName('LeftArrow'))%
            chosenSquare = 'Left';
            
        elseif keyCode(KbName('RightArrow'))%keyCode(KbName('RightArrow'))
            chosenSquare = 'Right';
            
        end
    end
    
    KbQueueStart();
    %Dot coordinates
    %     xcL = rect(3)*3/2+ex.horiOffsetL;
    %     ycL = rect(4)/2+ex.vertOffsetL;
    %% RIGHT SQUARE    
    % Grid Line coordinates
    % 45 deg line coordinates
    xRforLl =  3 * xc / 2-ex.lineW+ex.horiOffsetR;
    xRforLr =  3 * xc / 2+ex.lineW+ex.horiOffsetR;
    yRforLBot = yc + ex.lineH+ex.vertOffsetR;
    yRforLTop = yc - ex.lineH+ex.vertOffsetR;
    
    % 135 deg line coordinates
    xRbackRl = 3 * xc / 2-ex.lineW+ex.horiOffsetR;
    xRbackRr = 3 * xc / 2+ex.lineW+ex.horiOffsetR;
    yRbackRTop = yc- ex.lineH+ex.vertOffsetR;
    yRbackRBot = yc +ex.lineH+ex.vertOffsetR;
    
%    

    % Grid lines
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower right cross/grid
            Screen('DrawLines', w, [xRforLl+xshift, xRforLr+xshift; yRforLBot+yshift, yRforLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl+xshift, xRbackRr+xshift; yRbackRTop+yshift, yRbackRBot+yshift], ex.lineTh, ex.lineColor);
            %Upper right grid
            Screen('DrawLines', w, [xRforLl+xshift, xRforLr+xshift; yRforLBot-yshift, yRforLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl+xshift, xRbackRr+xshift; yRbackRTop-yshift, yRbackRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower left cross/grid
            Screen('DrawLines', w, [xRforLl-xshift, xRforLr-xshift; yRforLBot+yshift, yRforLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl-xshift, xRbackRr-xshift; yRbackRTop+yshift, yRbackRBot+yshift], ex.lineTh, ex.lineColor);
            
            %Upper left grid
            Screen('DrawLines', w, [xRforLl-xshift, xRforLr-xshift; yRforLBot-yshift, yRforLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl-xshift, xRbackRr-xshift; yRbackRTop-yshift, yRbackRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
%     %% square
    % Define the square parameters (position and size) 
    squareRectRight = [3 * xc / 2-ex.lineHsq/2+ex.horiOffsetR, yc - ex.lineHsq/2+ex.vertOffsetR, 3 * xc / 2+ex.lineHsq/2+ex.horiOffsetR, yc + ex.lineHsq/2+ex.vertOffsetR]; 
    % Draw the square 
    Screen('FrameRect', w, [0 0 0], squareRectRight, ex.lineTh); 
    
    %%%% FLIP
    Screen(w, 'Flip', 0);
    WaitSecs(0.5);
    %% LEFT SQUARE
    
    % Line coordinates
    % 45 deg line coordinates
    xLforl = xc / 2-ex.lineW+ex.horiOffsetL;
    xLforr = xc / 2+ex.lineW+ex.horiOffsetL;
    yLforBot = yc + ex.lineH+ex.vertOffsetL;
    yLforTop = yc - ex.lineH+ex.vertOffsetL;
    
    % 135 deg line coordinates
    xLbackl = xc / 2-ex.lineW+ex.horiOffsetL;
    xLbackr = xc / 2+ex.lineW+ex.horiOffsetL;
    yLbackTop = yc- ex.lineH+ex.vertOffsetL;
    yLbackBot = yc +ex.lineH+ex.vertOffsetL;
    
    % Define the square parameters (position and size)
    squareRectLeft = [xc / 2-ex.lineHsq/2+ex.horiOffsetL, yc - ex.lineHsq/2+ex.vertOffsetL, xc / 2+ex.lineHsq/2+ex.horiOffsetL, yc + ex.lineHsq/2+ex.vertOffsetL];
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower right cross/grid
            Screen('DrawLines', w, [xLforl+xshift, xLforr+xshift; yLforBot+yshift, yLforTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xLbackl+xshift, xLbackr+xshift; yLbackTop+yshift, yLbackBot+yshift], ex.lineTh, ex.lineColor);
            %Upper right grid
            Screen('DrawLines', w, [xLforl+xshift, xLforr+xshift; yLforBot-yshift, yLforTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xLbackl+xshift, xLbackr+xshift; yLbackTop-yshift, yLbackBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower left cross/grid
            Screen('DrawLines', w, [xLforl-xshift, xLforr-xshift; yLforBot+yshift, yLforTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xLbackl-xshift, xLbackr-xshift; yLbackTop+yshift, yLbackBot+yshift], ex.lineTh, ex.lineColor);
            
            %Upper left grid
            Screen('DrawLines', w, [xLforl-xshift, xLforr-xshift; yLforBot-yshift, yLforTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xLbackl-xshift, xLbackr-xshift; yLbackTop-yshift, yLbackBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    %right grid lines 
        % Grid lines
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower right cross/grid
            Screen('DrawLines', w, [xRforLl+xshift, xRforLr+xshift; yRforLBot+yshift, yRforLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl+xshift, xRbackRr+xshift; yRbackRTop+yshift, yRbackRBot+yshift], ex.lineTh, ex.lineColor);
            %Upper right grid
            Screen('DrawLines', w, [xRforLl+xshift, xRforLr+xshift; yRforLBot-yshift, yRforLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl+xshift, xRbackRr+xshift; yRbackRTop-yshift, yRbackRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower left cross/grid
            Screen('DrawLines', w, [xRforLl-xshift, xRforLr-xshift; yRforLBot+yshift, yRforLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl-xshift, xRbackRr-xshift; yRbackRTop+yshift, yRbackRBot+yshift], ex.lineTh, ex.lineColor);
            
            %Upper left grid
            Screen('DrawLines', w, [xRforLl-xshift, xRforLr-xshift; yRforLBot-yshift, yRforLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xRbackRl-xshift, xRbackRr-xshift; yRbackRTop-yshift, yRbackRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    % Draw the square
    Screen('FrameRect', w, [0 0 0], squareRectLeft, ex.lineTh);
    Screen('FrameRect', w, [0 0 0], squareRectRight, ex.lineTh); 
    %%%% FLIP
    Screen(w, 'Flip', 0);
    WaitSecs(0.5);
    Screen(w, 'Flip', 0);
    KbWait(deviceNumber,2)
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if contains(chosenSquare, 'Right')
        if find(firstPress,1) == KbName('8')
            ex.vertOffsetR = ex.vertOffsetR - 5;
            KbQueueFlush();
        elseif find(firstPress,1) == KbName('2')
            ex.vertOffsetR = ex.vertOffsetR + 5;
            KbQueueFlush();
        elseif find(firstPress,1) == KbName('4')
            ex.horiOffsetR = ex.horiOffsetR - 5;
            KbQueueFlush();
        elseif find(firstPress,1) == KbName('6')
            ex.horiOffsetR = ex.horiOffsetR + 5;
        elseif  find(firstPress,1) == KbName('Enter')
            break;
            KbQueueFlush();
        end
    elseif contains(chosenSquare, 'Left')
        if find(firstPress,1) == KbName('8')
            ex.vertOffsetL = ex.vertOffsetL - 5;
            KbQueueFlush();
        elseif find(firstPress,1) == KbName('2')
            ex.vertOffsetL = ex.vertOffsetL + 5;
            KbQueueFlush();
        elseif find(firstPress,1) == KbName('4')
            ex.horiOffsetL = ex.horiOffsetL - 5;
            KbQueueFlush();
        elseif find(firstPress,1) == KbName('6')
            ex.horiOffsetL = ex.horiOffsetL + 5;
            KbQueueFlush();
        elseif  find(firstPress,1) == KbName('Enter')
            break;
            KbQueueFlush();
        end
    end
    
end

vertOffsets = [ex.vertOffsetL; ex.vertOffsetR];
horiOffsets = [ex.horiOffsetL;ex.horiOffsetR];
KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');