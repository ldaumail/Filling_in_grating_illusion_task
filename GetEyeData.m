
if ~isempty(ET_ON)

% Get eye data
if Eyelink('NewFloatSampleAvailable') > 0;
	% Make the new sample into an event struct
    sample = Eyelink('NewestFloatSample');
    % Find current x and y
    gx = sample.gx(eye_used+1);
    gy = sample.gy(eye_used+1);
    pa = sample.pa(eye_used+1);
    
    g_cnt = g_cnt + 1;
    
    % If there is no missing data assign to mx, my
    if gx~=el.MISSING_DATA && gy~=el.MISSING_DATA && sample.pa(eye_used+1)>0;
    	gazex(g_cnt) = gx;
        gazey(g_cnt) = gy;
        gazepa(g_cnt) = pa;
        gazeTS(g_cnt) = GetSecs;                    
        gazeD(g_cnt) = sqrt( sum( ([gx gy]-[centerX centerY]).^2 ) );
    else
    	gazex(g_cnt) = NaN; % These are already preallocated with NaNs; but will be extended when g_cnt exceeds preallocated space
        gazey(g_cnt) = NaN;
        gazepa(g_cnt) = pa;
        gazeTS(g_cnt) = GetSecs;                    
        gazeD(g_cnt) = NaN;    
    end    
end
            
    cmap=cmap_gaze;
    if g_cnt>=params.nGazeAllowed   %   length(gazex)>=s.nGazeAllowed
        % if recent N trials were ALL outside the circle, abort the trial
        % and show Blue dots outside the circle
        % I.e. if there's no frames where the gaze was Inside the box
        if ( sum( gazeD(g_cnt-params.nGazeAllowed+1:end) <= params.ETbox_size_pix(1)/2 ) == 0);% s.nGazeAllowed )  % mean( gazeD(end-s.nGazeAllowed+1:end) ) > s.ETbox_size_pix        
            cmap = cmap_gaze_exceed; gaze_exceeded=1;
%             if ~isempty(AbortTrialWhenExceeded); abort=1; 
%                 break % exceeded=1;
%             end % ~isempty(s.AbortTrialWhenExceeded)
        else 
        end % if ( sum(gazeD(end-s.nGazeAllowed+1:end) >
    end % length(gazex)>=s.nGazeAllowed   
    
% Plot real-time gaze positions
if ~isempty(params.ShowRealTimeGaze) % length(gazex)>0 &&
    
	if g_cnt>=params.nGazetoShow %  length(gazex)>=s.nGazetoShow
        Screen('DrawDots', win, [gazex(g_cnt-params.nGazetoShow+1:g_cnt); gazey(g_cnt-params.nGazetoShow+1:g_cnt)], 10, [ repmat(cmap', 1, params.nGazetoShow);  linspace(10,255, params.nGazetoShow) ], [], 0);                           
%     	Screen('DrawDots', w, [gazex(end-s.nGazetoShow+1:end); gazey(end-s.nGazetoShow+1:end)], 1, [repmat([255 0 0]', 1, s.nGazetoShow);  linspace(10,255, s.nGazetoShow) ], [], 0);                           
	elseif g_cnt>0 % length(gazex)>0
        Screen('DrawDots', win, [gazex(1:g_cnt); gazey(1:g_cnt)], 10, [ repmat(cmap', 1, g_cnt);  linspace(.2,255, g_cnt) ], [], 0);            
%         Screen('DrawDots', w, [gazex; gazey], 1, [repmat([255 0 0]', 1, length(gazex));  linspace(.2,255, length(gazex)) ], [], 0);            
    end       
    
end % if ~isempty(ShowRealTimeGaze) 

end % if ~isempty(ET_ON)            