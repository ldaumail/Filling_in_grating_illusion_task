function [vertOffsets,horiOffsets] = findOffset_nonius_fixationV2(subject, session, debug)

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

%%%% Lines
ex.lineHdeg = 1/2;%1/sqrt(2);
ex.lineThdeg = 0.1;
ex.lineHdegSq = 5; %square size
%%% horizontal line
ex.horiLineWdeg = 0.7;

%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('2'))=1; % button box 1
responseKeys(KbName('4'))=1; % button box 2
responseKeys(KbName('6'))=1; % button box 3
responseKeys(KbName('8'))=1; % button box 3
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

%%%% Fixation dot
% ex.yFixOffset = 3*ex.ppd;
ex.lineH = ex.lineHdeg*ex.ppd;
ex.lineW = ex.lineHdeg*ex.ppd;
ex.lineTh = ex.lineThdeg*ex.ppd;
ex.lineHsq = ex.lineHdegSq*ex.ppd;
ex.lineColor = [0 0 0];
KbQueueCreate(deviceNumber,responseKeys);

%%%% find center Y
Screen('FillRect',w,ex.backgroundColor)
text = 'Find CENTER X,Y for left eye. \n\n 8 =  Up, 2 = Down, \n\n 6 = Right , 4: Left, \n\n Enter: Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',w,text));
xc =rect(3)/2;
yc = rect(4)/2;
DrawFormattedText(w, text, xc/5, yc/2);
DrawFormattedText(w, text, xc+xc/5, yc/2);
Screen(w, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

%%%%%%%%
% RUN
%%%%%%%%

while 1
    
    %Dot coordinates
    %     xcL = rect(3)*3/2+ex.horiOffsetL;
    %     ycL = rect(4)/2+ex.vertOffsetL;
    
    % Grid Line coordinates
    % 45 deg line coordinates
    xlineLl =  3 * xc / 2-ex.lineW+ex.horiOffsetR;
    xlineLr =  3 * xc / 2+ex.lineW+ex.horiOffsetR;
    ylineLBot = yc + ex.lineH+ex.vertOffsetR;
    ylineLTop = yc - ex.lineH+ex.vertOffsetR;
    
    % 135 deg line coordinates
    xlineRl = 3 * xc / 2-ex.lineW+ex.horiOffsetR;
    xlineRr = 3 * xc / 2+ex.lineW+ex.horiOffsetR;
    ylineRTop = yc- ex.lineH+ex.vertOffsetR;
    ylineRBot = yc +ex.lineH+ex.vertOffsetR;
    
    KbQueueStart();
    % Grid lines
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower right cross/grid
            Screen('DrawLines', w, [xlineLl+xshift, xlineLr+xshift; ylineLBot+yshift, ylineLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl+xshift, xlineRr+xshift; ylineRTop+yshift, ylineRBot+yshift], ex.lineTh, ex.lineColor);
            %Upper right grid
            Screen('DrawLines', w, [xlineLl+xshift, xlineLr+xshift; ylineLBot-yshift, ylineLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl+xshift, xlineRr+xshift; ylineRTop-yshift, ylineRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower left cross/grid
            Screen('DrawLines', w, [xlineLl-xshift, xlineLr-xshift; ylineLBot+yshift, ylineLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl-xshift, xlineRr-xshift; ylineRTop+yshift, ylineRBot+yshift], ex.lineTh, ex.lineColor);
            
            %Upper left grid
            Screen('DrawLines', w, [xlineLl-xshift, xlineLr-xshift; ylineLBot-yshift, ylineLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl-xshift, xlineRr-xshift; ylineRTop-yshift, ylineRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    %% square
    %left edge
    xlineLsq =  3 * xc / 2-ex.lineHsq/2+ex.horiOffsetR;
    ylineLBotsq = yc + ex.lineHsq/2+ex.vertOffsetR;
    ylineLTopsq = yc - ex.lineHsq/2+ex.vertOffsetR;
    Screen('DrawLines', w, [xlineLsq, xlineLsq; ylineLBotsq, ylineLTopsq], ex.lineTh, ex.lineColor);
    %right edge
    xlineRsq =  3 * xc / 2+ex.lineHsq/2+ex.horiOffsetR;
    ylineRBotsq = yc + ex.lineHsq/2+ex.vertOffsetR;
    ylineRTopsq = yc - ex.lineHsq/2+ex.vertOffsetR;
    Screen('DrawLines', w, [xlineRsq, xlineRsq; ylineRBotsq, ylineRTopsq], ex.lineTh, ex.lineColor);
    %Lower edge
    xlineBlsq =  3 * xc / 2-ex.lineHsq/2+ex.horiOffsetR;
    xlineBrsq =  3 * xc / 2+ex.lineHsq/2+ex.horiOffsetR;
    ylineBsq = yc + ex.lineHsq/2+ex.vertOffsetR;
    Screen('DrawLines', w, [xlineBlsq, xlineBrsq; ylineBsq, ylineBsq], ex.lineTh, ex.lineColor);
    %Upper edge
    xlineTlsq =  3 * xc / 2-ex.lineHsq/2+ex.horiOffsetR;
    xlineTrsq =  3 * xc / 2+ex.lineHsq/2+ex.horiOffsetR;
    ylineTsq = yc - ex.lineHsq/2+ex.vertOffsetR;
    Screen('DrawLines', w, [xlineTlsq, xlineTrsq; ylineTsq, ylineTsq], ex.lineTh, ex.lineColor);
    
    
    %%%% FLIP
    Screen(w, 'Flip', 0);
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
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
        KbQueueFlush();
        
    elseif find(firstPress,1) == KbName('Enter')
        break;
        KbQueueFlush();
    end
end

vertOffsets = ex.vertOffsetR;
horiOffsets = ex.horiOffsetR;

Screen('FillRect',w,ex.backgroundColor)
text = 'Find CENTER X,Y for right eye. \n\n 8 =  Up, 2 = Down, \n\n 6 = Right , 4: Left, \n\n Enter: Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',w,text));
DrawFormattedText(w, text, xc/5, yc/2);
DrawFormattedText(w, text, xc+xc/5, yc/2);
Screen(w, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

while 1
    
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower right cross/grid
            Screen('DrawLines', w, [xlineLl+xshift, xlineLr+xshift; ylineLBot+yshift, ylineLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl+xshift, xlineRr+xshift; ylineRTop+yshift, ylineRBot+yshift], ex.lineTh, ex.lineColor);
            %Upper right grid
            Screen('DrawLines', w, [xlineLl+xshift, xlineLr+xshift; ylineLBot-yshift, ylineLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl+xshift, xlineRr+xshift; ylineRTop-yshift, ylineRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    for c = 1:3
        xshift = 2*ex.lineH*(c-1);
        for r=1:3
            yshift = 2*ex.lineW*(r-1);
            %lower left cross/grid
            Screen('DrawLines', w, [xlineLl-xshift, xlineLr-xshift; ylineLBot+yshift, ylineLTop+yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl-xshift, xlineRr-xshift; ylineRTop+yshift, ylineRBot+yshift], ex.lineTh, ex.lineColor);
            
            %Upper left grid
            Screen('DrawLines', w, [xlineLl-xshift, xlineLr-xshift; ylineLBot-yshift, ylineLTop-yshift], ex.lineTh, ex.lineColor);
            Screen('DrawLines', w, [xlineRl-xshift, xlineRr-xshift; ylineRTop-yshift, ylineRBot-yshift], ex.lineTh, ex.lineColor);
            
        end
    end
    %% square
    %left edge
    
    Screen('DrawLines', w, [xlineLsq, xlineLsq; ylineLBotsq, ylineLTopsq], ex.lineTh, ex.lineColor);
    %right edge
    Screen('DrawLines', w, [xlineRsq, xlineRsq; ylineRBotsq, ylineRTopsq], ex.lineTh, ex.lineColor);
    %Lower edge
    Screen('DrawLines', w, [xlineBlsq, xlineBrsq; ylineBsq, ylineBsq], ex.lineTh, ex.lineColor);
    %Upper edge
    Screen('DrawLines', w, [xlineTlsq, xlineTrsq; ylineTsq, ylineTsq], ex.lineTh, ex.lineColor);
    
    
    % Line coordinates
    % 45 deg line coordinates
    xlineMl = xc / 2-ex.lineW+ex.horiOffsetL;
    xlineMr = xc / 2+ex.lineW+ex.horiOffsetL;
    ylineMBot = yc + ex.lineH+ex.vertOffsetL;
    ylineMTop = yc - ex.lineH+ex.vertOffsetL;
    
    % 135 deg line coordinates
    xlineNl = xc / 2-ex.lineW+ex.horiOffsetL;
    xlineNr = xc / 2+ex.lineW+ex.horiOffsetL;
    ylineNTop = yc- ex.lineH+ex.vertOffsetL;
    ylineNBot = yc +ex.lineH+ex.vertOffsetL;
    
    
    %left edge
    xlineLLsq =  xc / 2-ex.lineHsq/2+ex.horiOffsetL;
    ylineLLBotsq = yc + ex.lineHsq/2+ex.vertOffsetL;
    ylineLLTopsq = yc - ex.lineHsq/2+ex.vertOffsetL;
    %right edge
    xlineLRsq =   xc / 2+ex.lineHsq/2+ex.horiOffsetL;
    ylineLRBotsq = yc + ex.lineHsq/2+ex.vertOffsetL;
    ylineLRTopsq = yc - ex.lineHsq/2+ex.vertOffsetL;
    %Lower edge
    xlineLBlsq =  xc / 2-ex.lineHsq/2+ex.horiOffsetL;
    xlineLBrsq =  xc / 2+ex.lineHsq/2+ex.horiOffsetL;
    ylineLBsq = yc + ex.lineHsq/2+ex.vertOffsetL;
    %Upper edge
    xlineLTlsq =  xc / 2-ex.lineHsq/2+ex.horiOffsetL;
    xlineLTrsq =  xc / 2+ex.lineHsq/2+ex.horiOffsetL;
    ylineLTsq = yc - ex.lineHsq/2+ex.vertOffsetL;
    
    KbQueueStart();
    %     Screen('DrawDots', w, [xcL ycL], ex.fixSize, [255 255 255], [], 2); %Left dot
    
    
    %     Screen('DrawDots', w, [xcL*3/2 ycL], ex.fixSize, [255 255 255], [], 2); %Draw Left Dot
    
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if find(firstPress,1) == KbName('8')
        ex.vertOffsetL = ex.vertOffsetL - 5;
        KbQueueFlush();
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower right cross/grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                %Upper right grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower left cross/grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                
                %Upper left grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        Screen('DrawLines', w, [xlineLLsq, xlineLLsq; ylineLLBotsq, ylineLLTopsq], ex.lineTh, ex.lineColor);
        %right edge
        Screen('DrawLines', w, [xlineLRsq, xlineLRsq; ylineLRBotsq, ylineLRTopsq], ex.lineTh, ex.lineColor);
        %Lower edge
        Screen('DrawLines', w, [xlineLBlsq, xlineLBrsq; ylineLBsq, ylineLBsq], ex.lineTh, ex.lineColor);
        %Upper edge
        Screen('DrawLines', w, [xlineLTlsq, xlineLTrsq; ylineLTsq, ylineLTsq], ex.lineTh, ex.lineColor);
        
        %         Screen('DrawDots', w, [xcR/2 ycR], ex.fixSize, [255 255 255], [], 2);
    elseif find(firstPress,1) == KbName('2')
        ex.vertOffsetL = ex.vertOffsetL + 5;
        KbQueueFlush();
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower right cross/grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                %Upper right grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower left cross/grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                
                %Upper left grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        Screen('DrawLines', w, [xlineLLsq, xlineLLsq; ylineLLBotsq, ylineLLTopsq], ex.lineTh, ex.lineColor);
        %right edge
        Screen('DrawLines', w, [xlineLRsq, xlineLRsq; ylineLRBotsq, ylineLRTopsq], ex.lineTh, ex.lineColor);
        %Lower edge
        Screen('DrawLines', w, [xlineLBlsq, xlineLBrsq; ylineLBsq, ylineLBsq], ex.lineTh, ex.lineColor);
        %Upper edge
        Screen('DrawLines', w, [xlineLTlsq, xlineLTrsq; ylineLTsq, ylineLTsq], ex.lineTh, ex.lineColor);
        
        %         Screen('DrawDots', w, [xcR/2 ycR], ex.fixSize, [255 255 255], [], 2);
    elseif find(firstPress,1) == KbName('4')
        ex.horiOffsetL = ex.horiOffsetL - 5;
        KbQueueFlush();
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower right cross/grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                %Upper right grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower left cross/grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                
                %Upper left grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        Screen('DrawLines', w, [xlineLLsq, xlineLLsq; ylineLLBotsq, ylineLLTopsq], ex.lineTh, ex.lineColor);
        %right edge
        Screen('DrawLines', w, [xlineLRsq, xlineLRsq; ylineLRBotsq, ylineLRTopsq], ex.lineTh, ex.lineColor);
        %Lower edge
        Screen('DrawLines', w, [xlineLBlsq, xlineLBrsq; ylineLBsq, ylineLBsq], ex.lineTh, ex.lineColor);
        %Upper edge
        Screen('DrawLines', w, [xlineLTlsq, xlineLTrsq; ylineLTsq, ylineLTsq], ex.lineTh, ex.lineColor);
        
        %         Screen('DrawDots', w, [xcR/2 ycR], ex.fixSize, [255 255 255], [], 2);
    elseif find(firstPress,1) == KbName('6')
        ex.horiOffsetL = ex.horiOffsetL + 5;
        KbQueueFlush();
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower right cross/grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                %Upper right grid
                Screen('DrawLines', w, [xlineMl+xshift, xlineMr+xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl+xshift, xlineNr+xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        for c = 1:3
            xshift = 2*ex.lineH*(c-1);
            for r=1:3
                yshift = 2*ex.lineW*(r-1);
                %lower left cross/grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot+yshift, ylineMTop+yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop+yshift, ylineNBot+yshift], ex.lineTh, ex.lineColor);
                
                %Upper left grid
                Screen('DrawLines', w, [xlineMl-xshift, xlineMr-xshift; ylineMBot-yshift, ylineMTop-yshift], ex.lineTh, ex.lineColor);
                Screen('DrawLines', w, [xlineNl-xshift, xlineNr-xshift; ylineNTop-yshift, ylineNBot-yshift], ex.lineTh, ex.lineColor);
                
            end
        end
        Screen('DrawLines', w, [xlineLLsq, xlineLLsq; ylineLLBotsq, ylineLLTopsq], ex.lineTh, ex.lineColor);
        %right edge
        Screen('DrawLines', w, [xlineLRsq, xlineLRsq; ylineLRBotsq, ylineLRTopsq], ex.lineTh, ex.lineColor);
        %Lower edge
        Screen('DrawLines', w, [xlineLBlsq, xlineLBrsq; ylineLBsq, ylineLBsq], ex.lineTh, ex.lineColor);
        %Upper edge
        Screen('DrawLines', w, [xlineLTlsq, xlineLTrsq; ylineLTsq, ylineLTsq], ex.lineTh, ex.lineColor);
        
        %         Screen('DrawDots', w, [xcR/2 ycR], ex.fixSize, [255 255 255], [], 2);
    elseif find(firstPress,1) == KbName('Enter')
        break;
        KbQueueFlush();
    end
    %%%% FLIP
    Screen(w, 'Flip', 0);
end
vertOffsets = [vertOffsets; ex.vertOffsetL];
horiOffsets = [horiOffsets; ex.horiOffsetL];

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

