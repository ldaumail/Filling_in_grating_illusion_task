function scrSz = findScreenSize_LD(subject, session, debug)

%%%%%%%%%%%%%%%%%%%
% BASIC SETTINGS
%%%%%%%%%%%%%%%%%%%

%%%% resolution
if debug == 1
	SetResolution(max(Screen('Screens')),1920,1080,0); % laptop
else
    SetResolution(max(Screen('Screens')),1024,768,60); % scanner
end

%%%% keyboards
[keyboardIndices, productNames, ~] = GetKeyboardIndices ;
deviceNumber = keyboardIndices(1);

%%%% input this at the beginning of the scan session for 7T
scSize.topDegY = 4; 
scSize.bottomDegY = 4; 
scSize.leftDegX = 8; 
scSize.rightDegX = 8; 

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
responseKeys(KbName('4$'))=1; % button box 4
responseKeys(KbName('5%'))=1; % button box 5
responseKeys(KbName('6^'))=1; % button box 6
responseKeys(KbName('7&'))=1; % button box 7
responseKeys(KbName('8*'))=1; % button box 8
responseKeys(KbName('9('))=1; % button box 9

%%%%%%%%%%%%%%%%
% OPEN SCREEN 
%%%%%%%%%%%%%%%%

HideCursor;
Priority(9); % highest priority 

%%%% open screen
Screen('Preference', 'SkipSyncTests', 1);
screen=max(Screen('Screens'));
%[win, rect]=Screen('OpenWindow',screen,180,[100 100 600 400],[],[],[],[],kPsychNeed32BPCFloat); 
[win, rect]=Screen('OpenWindow',screen,180,[],[],[],[],[],kPsychNeed32BPCFloat);
Screen(win, 'TextSize', 32);
xc = rect(3)/2; 
yc = rect(4)/2;

%%%% scale the stims for the screen
params.ppd = pi* rect(3) / (atan(params.screenWidth/params.viewingDist/2)) / 360;
params.fixSize = round(params.fixSizeDeg*params.ppd);
params.littleFix = round(params.littleFixDeg*params.ppd);
scSize.topY = scSize.topDegY*params.ppd;
scSize.bottomY = scSize.bottomDegY*params.ppd; 
scSize.leftX = scSize.leftDegX*params.ppd; 
scSize.rightX = scSize.rightDegX*params.ppd; 

KbQueueCreate(deviceNumber,responseKeys);

%%%% find center Y
Screen('FillRect',win,params.backgroundColor)
text = 'Adjust screen size to visible field of view.\n\n 1/2 = Left, 3/4 = Right, 5/6 = Top, 7/8 = Bottom, 9 = Accept\n\nPress any button to start';
%width = RectWidth(Screen('TextBounds',win,text));
DrawFormattedText(win, text, 'center', 'center');
Screen(win, 'Flip', 0);
WaitSecs(1); KbPressWait; KbQueueFlush();

%%%%%%%%
% RUN 
%%%%%%%%

while 1

    KbQueueStart();
    %%%% draw screen square
    Screen('FillRect', win,params.boxColor, [xc-round(scSize.leftX) yc-round(scSize.topY) xc+round(scSize.rightX) yc+round(scSize.bottomY)]);    
    
    %%%% draw fixation
    Screen('FillOval', win,[0 0 0], [xc-round(params.fixSize/2+params.outerFixPixels ) yc-round(params.fixSize/2+params.outerFixPixels ) xc+round(params.fixSize/2+params.outerFixPixels ) yc+round(params.fixSize/2+params.outerFixPixels )]); % black fixation ring
    Screen('FillOval', win,[255 255 255], [xc-round(params.fixSize/2) yc-round(params.fixSize/2) xc+round(params.fixSize/2) yc+round(params.fixSize/2)]); % white fixation ring
    Screen('FillOval', win,[0 0 0], [xc-round(params.littleFix/2) yc-round(params.littleFix/2) xc+round(params.littleFix/2) yc+round(params.littleFix/2)]); % black fixation dot
    
    %%%% FLIP 
    Screen(win, 'Flip', 0);
    
    %%%% check for responses
    [pressed, firstPress]= KbQueueCheck();
    if find(firstPress,1) == KbName('1!')
        scSize.leftX = scSize.leftX - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('2@')
        scSize.leftX = scSize.leftX + 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('3#')
        scSize.rightX = scSize.rightX - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('4$')
        scSize.rightX = scSize.rightX + 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('5%')
        scSize.topY = scSize.topY - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('6^')
        scSize.topY = scSize.topY + 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('7&')
        scSize.bottomY = scSize.bottomY - 5;
        KbQueueFlush();
    elseif find(firstPress,1) == KbName('8*')
        scSize.bottomY = scSize.bottomY + 5;
        KbQueueFlush();
     elseif find(firstPress,1) == KbName('9(')
        break;
        KbQueueFlush();
    end
end

scrSz = [xc-round(scSize.leftX) yc-round(scSize.topY) xc+round(scSize.rightX) yc+round(scSize.bottomY)];%./params.ppd;

%%%%%%%%%%%%%%%%%%
% DONE! WRAP UP
%%%%%%%%%%%%%%%%%%

KbQueueRelease();
ShowCursor;
Screen('Close');
Screen('CloseAll');

% save
savedir = fullfile(pwd,'data',subject,session,'findVertOffset');
if ~exist(savedir); mkdir(savedir); end
savename=fullfile(savedir, strcat(subject,'_screenSize.mat'));   
save(savename,'scSize')

end

