global el EyeData

if ET
    errormes=Eyelink('CheckRecording');
    %         if(errormes~=0)
    %         end
    if Eyelink('NewFloatSampleAvailable') > 0
        % get the sample in the form of an event structure
        evt = Eyelink('NewestFloatSample');
        x = evt.gx(eye_used+1);
        y = evt.gy(eye_used+1);
        a = evt.pa(eye_used+1);
        gcnt = gcnt + 1;
        
        % If there is no missing data assign to mx, my
        if x~=el.MISSING_DATA && y~=el.MISSING_DATA && a>0
            EyeData.mx(gcnt)= x;
            EyeData.my(gcnt)= y;
            EyeData.ma(gcnt)= a;
            EyeData.FixDoneT(gcnt) = GetSecs(); % record each time we get valid eyetracking data point
            EyeData.gazeD(gcnt) = sqrt( sum( ([x y]-[xc yc]).^2 ) );
            
        else
            EyeData.mx(gcnt)= NaN;%These are already preallocated with NaNs; but will be extended when gcnt exceeds preallocated space
            EyeData.my(gcnt)= NaN;
            EyeData.ma(gcnt)= a;
            EyeData.FixDoneT(gcnt) = GetSecs(); % record each time we get valid eyetracking data point
            EyeData.gazeD(gcnt) = NaN;
        end
        % check whether fix the center
        if ~ (isempty(x) || isempty(y) || isempty(a))
            if x~=el.MISSING_DATA && y~=el.MISSING_DATA && a>0 % do we have valid data and is the pupil visible?
                dist_from_fix = round(( (x-xc)^2 + (y-yc)^2 ) ^ (1/2));
                
                if dist_from_fix <= 50
                    EyeData.Fixated(gcnt) = 1;
                    
                    %                     if currPosID==1 % until the motion turns on, keep noting the difference in time from inital start to 'now'
                    %                         InitialFixDur(iTrial) = fixT-EyeStart(iTrial);
                    %                     end
                else
                    EyeData.Fixated(gcnt) = 0;
                end
            end
        else
            EyeData.Fixated(gcnt) = 0;
        end
        
       % Plot real-time gaze positions
if ~isempty(ex.ShowRealTimeGaze) % length(gazex)>0 &&
    
	if gcnt>=params.nGazetoShow %  length(gazex)>=s.nGazetoShow
        Screen('DrawDots', win, [EyeData.mx(gcnt-ex.nGazetoShow+1:gcnt); EyeData.my(gcnt-ex.nGazetoShow+1:gcnt)], 10, [ repmat(cmap', 1, ex.nGazetoShow);  linspace(10,255, ex.nGazetoShow) ], [], 0);                           
%     	Screen('DrawDots', w, [gazex(end-s.nGazetoShow+1:end); gazey(end-s.nGazetoShow+1:end)], 1, [repmat([255 0 0]', 1, s.nGazetoShow);  linspace(10,255, s.nGazetoShow) ], [], 0);                           
	elseif gcnt>0 % length(gazex)>0
        Screen('DrawDots', win, [EyeData.mx(1:gcnt); EyeData.my(1:gcnt)], 10, [ repmat(cmap', 1, gcnt);  linspace(.2,255, gcnt) ], [], 0);            
%         Screen('DrawDots', w, [gazex; gazey], 1, [repmat([255 0 0]', 1, length(gazex));  linspace(.2,255, length(gazex)) ], [], 0);            
    end       
    
end % if ~isempty(ShowRealTimeGaze) 

 
        %             if ~motionOn && t>(TrialStart(iTrial)+2)
        %                 % fixation break when...
        %                 % motion hasn't started and the program has waited
        %                 % for a really long time
        %                 if ~Fixated
        %                     fixBreak = 1;
        %                 end
        %             else
        %                 fixBreak = 0;
        %
        %             end
    end
else
    fixBreak=0;
end
