global el EyeData

if ET
    errormes=Eyelink('CheckRecording');
%         if(errormes~=0)
%         end
    if Eyelink('NewFloatSampleAvailable') > 0
        % get the sample in the form of an event structure
        evt = Eyelink('NewestFloatSample');
        x = evt.gx(evt.gx~=el.MISSING_DATA);
        y = evt.gy(evt.gy~=el.MISSING_DATA);
        a = evt.pa(evt.pa~=-el.MISSING_DATA);
        fixT = GetSecs(); %doesn't account for blinks, so fix this to get accurate timestamps
        EyeData.mx{1}=[x];
        EyeData.my{1}=[y];
        EyeData.ma{1}=[a];
        EyeData.FixDoneT{1} = [fixT]; % record each time we get valid eyetracking data point

        % check whether fix the center
        if ~ (isempty(x) || isempty(y) || isempty(a))
            if x~=el.MISSING_DATA && y~=el.MISSING_DATA && a>0 % do we have valid data and is the pupil visible?
                dist_from_fix = round(( (x-(rect(3)/2))^2 + (y-(rect(4)/2))^2 ) ^ (1/2));

                if dist_from_fix <= 50
                    Fixated = 1;
                    
%                     if currPosID==1 % until the motion turns on, keep noting the difference in time from inital start to 'now'
%                         InitialFixDur(iTrial) = fixT-EyeStart(iTrial);
%                     end
                else
                    Fixated = 0;
                end
            end
        else
            Fixated = 0;
        end


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
