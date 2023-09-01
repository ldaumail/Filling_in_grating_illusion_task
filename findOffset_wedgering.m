function [vertOffset,horiOffset] = findOffset_wedgering(subject, session, debug)

%%%%%%%%%%%%%%%%%%%
% BASIC SETTINGS
%%%%%%%%%%%%%%%%%%%

%%%% resolution
if debug == 1
	SetResolution(max(Screen('Screens')), 2880, 1800, 0); % laptop
else
    SetResolution(max(Screen('Screens')),1024,768,60); % scanner
end
experiment.runNum = input('Run number :');

%%%% keyboards
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);

%%%% input this at the beginning of the scan session for 7T
scSize.vertOffset = 0;                      % vertical offset from FindScreenSize.m
scSize.horiOffset = 0;


scSize.screenDeg = 14;                      % pRF screen size
% scSize.screenDegY = 8; %figSizeDeg;
% scSize.screenDegX = 16; %figSizeDeg*3;%9;                      % pRF screen size
% scSize.centerDeg = 4;                      % center circle size in diameter
% scSize.surroundDeg = 12;                      % surround circle size in diameter

%%%% scales all of the stimuli in DVA to the screensize
params.screenWidth = 17;                    % in cm; laptop=27.5, office=43, 3Tb=19, 7T=17, miniHelm=39;
params.viewingDist = 48;                    % in cm; 3Tb/office=43, 7T=48, miniHelm=57;
params.fixSizeDeg =  .5;                    % in degrees, the size of the biggest white dot in the fixation
params.littleFixDeg = params.fixSizeDeg* .5;% proportion of the fixSizeDeg occupied by the smaller black dot
params.outerFixPixels = 2;                  % in pixels, the black ring around fixation
params.backgroundColor = [128 128 128];     % color
params.boxColor = params.backgroundColor * .7;

%%%% response listening - so we don't read in scanner triggers!
responseKeys = zeros(1,256);
responseKeys(KbName('1!'))=1; % button box 1
responseKeys(KbName('2@'))=1; % button box 2
responseKeys(KbName('3#'))=1; % button box 3


%%%%%%%%%%%%%%%%
% OPEN SCREEN 
%%%%%%%%%%%%%%%%

HideCursor;
Priority(9); % highest priority 

%%%% open screen
Screen('Preference', 'SkipSyncTests', 1);
screen=max(Screen('Screens'));
[win, rect]=Screen('OpenWindow',screen,180,[],[],[],[],[],kPsychNeed32BPCFloat);
Screen(win, 'TextSize', 32);


% xc = rect(3)/2; 
% yc = rect(4)/2;

%Load the gamma table of the specified screen
load 7T_Sam.mat
Screen('LoadNormalizedGammaTable', win, linearizedCLUT);

%%%% scale the stims for the screen
params.ppd = pi* rect(3) / (atan(params.screenWidth/params.viewingDist/2)) / 360;
params.fixSize = round(params.fixSizeDeg*params.ppd);
params.littleFix = round(params.littleFixDeg*params.ppd);
scSize.screenRadius = round(scSize.screenDeg*params.ppd/2);
% scSize.screenY = scSize.screenDeg*params.ppd;

%%%% Fixation dot
% params.yFixOffset = 3*params.ppd;

KbQueueCreate(deviceNumber,responseKeys);

%%%% find center Y
Screen('FillRect',win,params.backgroundColor)
text = 'Find CENTER Y. 1 = Up, 2 = Down, 3 = Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',win,text));
DrawFormattedText(win, text, 'center', 'center');
Screen(win, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

%%%%%%%%
% RUN 
%%%%%%%%

while 1
    
%     xS = rect(3)/2+scSize.horiOffset;
%     yS = rect(4)/2+scSize.vertOffset;
    
    xc = rect(3)/2+scSize.horiOffset;
    yc = rect(4)/2+scSize.vertOffset; %+params.yFixOffset;
    
    KbQueueStart();    
    %%%% draw screen circle
%     Screen('FillRect', win,params.boxColor, [xS-round(scSize.screenX/2) yS-round(scSize.screenY/2) xS+round(scSize.screenX/2) yS+round(scSize.screenY/2)]);    
    Screen('FrameOval', win,params.boxColor, [xc-scSize.screenRadius yc-scSize.screenRadius xc+scSize.screenRadius yc+scSize.screenRadius], 10 , 10);    
%     Screen('DrawTexture', win, patternTex{allTexNum(n)}, [], [xc-stimRadius yc-stimRadius xc+stimRadius yc+stimRadius]);
        
    %%%% draw fixation
    Screen('FillOval', win,[0 0 0], [xc-round(params.fixSize/2+params.outerFixPixels ) yc-round(params.fixSize/2+params.outerFixPixels ) xc+round(params.fixSize/2+params.outerFixPixels ) yc+round(params.fixSize/2+params.outerFixPixels )]); % black fixation ring
    Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
    Screen('FillOval', win,[0 0 0], [xc-round(params.littleFix/2) yc-round(params.littleFix/2) xc+round(params.littleFix/2) yc+round(params.littleFix/2)]); % black fixation dot

    
    %%%% FLIP 
    Screen(win, 'Flip', 0);
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if find(firstPress,1) == KbName('1!')
        scSize.vertOffset = scSize.vertOffset - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('2@')
        scSize.vertOffset = scSize.vertOffset + 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('3#')
        break;
        KbQueueFlush();
    end
end

vertOffset = scSize.vertOffset;

Screen('FillRect',win,params.backgroundColor)
text = 'Find CENTER X. 1 = Left, 2 = Right, 3 = Accept\n\nPress any button to start';
width = RectWidth(Screen('TextBounds',win,text));
DrawFormattedText(win, text, 'center', 'center');
Screen(win, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

while 1
    
    xc = rect(3)/2+scSize.horiOffset;
    yc = rect(4)/2+scSize.vertOffset;
    
    KbQueueStart();
    %%%% draw screen square
    Screen('FillRect', win,params.boxColor, [xc-round(scSize.screenX/2) yc-round(scSize.screenY/2) xc+round(scSize.screenX/2) yc+round(scSize.screenY/2)]);    
    
    %%%% draw fixation
    Screen('FillOval', win,[0 0 0], [xc-round(params.fixSize/2+params.outerFixPixels ) yc-round(params.fixSize/2+params.outerFixPixels ) xc+round(params.fixSize/2+params.outerFixPixels ) yc+round(params.fixSize/2+params.outerFixPixels )]); % black fixation ring
    Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
    Screen('FillOval', win,[0 0 0], [xc-round(params.littleFix/2) yc-round(params.littleFix/2) xc+round(params.littleFix/2) yc+round(params.littleFix/2)]); % black fixation dot

    %%%% FLIP 
    Screen(win, 'Flip', 0);
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if find(firstPress,1) == KbName('1!')
        scSize.horiOffset = scSize.horiOffset - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('2@')
        scSize.horiOffset = scSize.horiOffset + 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('3#')
        break;
        KbQueueFlush();
    end
end

horiOffset = scSize.horiOffset;

%%%%%%%%%%%%%%%%%%
% DONE! WRAP UP
%%%%%%%%%%%%%%%%%%

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');

% save
savedir = fullfile(pwd,'data',subject,session,'findOffset_figureGround');
if ~exist(savedir); mkdir(savedir); end
savename=fullfile(savedir, strcat(subject,'_screenSize.mat'));   
save(savename,'scSize')

end
