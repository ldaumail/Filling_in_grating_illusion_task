function pursuitDV

% Import DV sample report
[~]=menu('Choose a DV sample report location','...');
h = waitbar(0/10, 'Please wait...');
dat = tdfread;
delete(h)

% Open output text file and write header
[~]=menu('Choose where to save output report','...');
DirName = uigetdir;
f = fopen(fullfile(DirName,'results.txt'),'w');
header=[{'SESSION_LABEL'},{'TRIAL_INDEX'},{'EYE_USED'},{'FIX_INDEX'},{'START_EDF_TIME'},{'END_EDF_TIME'},{'DURATION'},{'VEL_GAIN'},{'RMSE_GAZE'}];
fprintf(f,'%-10s\t',header{1:end-1});
fprintf(f,'%-10s\n',header{end});

indFixL = find(~isnan(dat.LEFT_FIX_INDEX)==1);
indFixR = find(~isnan(dat.RIGHT_FIX_INDEX)==1);

if ~isempty(indFixL)
    % indAll column 1 = first index, column 2 = second index, column 3 = fixation number
    indAllL=indFixL([1,find(diff(indFixL)~=1)'+1]);
    indAllL(:,2)=indFixL([find(diff(indFixL)~=1)',length(indFixL)]);
    indAllL(:,3)=dat.LEFT_FIX_INDEX(indAllL(:,1)')';
    
        % loop through number of fixations
    for i = 1: length(indAllL) 
        %velGain: divide gaze angular velocity by target angular velocity
        targAvgVelX = nanmean(dat.TARGET_VELOCITY_X_TARG1(indAllL(i,1):indAllL(i,2)));
        targAvgVelY = nanmean(dat.TARGET_VELOCITY_Y_TARG1(indAllL(i,1):indAllL(i,2)));
        [~, targAvgVel]=cart2pol(targAvgVelX, targAvgVelY);
        eyeAvgVelX = nanmean(dat.LEFT_VELOCITY_X(indAllL(i,1):indAllL(i,2)));
        eyeAvgVelY = nanmean(dat.LEFT_VELOCITY_Y(indAllL(i,1):indAllL(i,2)));
        [~, eyeAvgVel]=cart2pol(eyeAvgVelX, eyeAvgVelY);
        velGain = eyeAvgVel/targAvgVel;
        
        %RMSE_Gaze:sqrt of eucledian eye pos - targ pos in degrees ^2
        targPosX = dat.TARGET_X_TARG1(indAllL(i,1):indAllL(i,2));
        targPosY = dat.TARGET_Y_TARG1(indAllL(i,1):indAllL(i,2));
        eyePosX = dat.LEFT_GAZE_X(indAllL(i,1):indAllL(i,2));
        eyePosY = dat.LEFT_GAZE_Y(indAllL(i,1):indAllL(i,2));
        diffPosX = 1/nanmean(dat.RESOLUTION_X(indAllL(i,1):indAllL(i,2))).*(eyePosX - targPosX);
        diffPosY = 1/nanmean(dat.RESOLUTION_X(indAllL(i,1):indAllL(i,2))).*(eyePosY - targPosY);
        [~, diffPos]=cart2pol(diffPosX, diffPosY);
        RMSE_Gaze = sqrt(nanmean(diffPos).^2);
        
        %EDF time
        startEDF = dat.TIMESTAMP(indAllL(i,1));
        endEDF = dat.TIMESTAMP(indAllL(i,2));
        dur = (endEDF - startEDF) + (dat.TIMESTAMP(indAllL(i,1)+1) - dat.TIMESTAMP(indAllL(i,1)));
                               
        %save to text file
        sessName = dat.DATA_FILE(indAllL(i,1),:);        
        trialIn = dat.TRIAL_INDEX(indAllL(i,1),:);
        eye='L';
        fprintf(f,'%s\t',sessName);
        fprintf(f,'%d\t',trialIn);
        fprintf(f,'%s\t',eye);
        fprintf(f,'%d\t',indAllL(i,3));
        
        fprintf(f,'%d\t',startEDF);
        fprintf(f,'%d\t',endEDF);
        fprintf(f,'%d\t',dur);
        
        fprintf(f,'%-10.4f\t',velGain);
        fprintf(f,'%-10.4f\t\n',RMSE_Gaze);
    end
end

if ~isempty(indFixR)
    % indAll column 1 = first index, column 2 = second index, column 3 = fixation number
    indAllR=indFixR([1,find(diff(indFixR)~=1)'+1]);
    indAllR(:,2)=indFixR([find(diff(indFixR)~=1)',length(indFixR)]);
    indAllR(:,3)=dat.RIGHT_FIX_INDEX(indAllR(:,1)')';
    
     % loop through number of fixations
    for i = 1: length(indAllR) 
        %velGain: divide gaze angular velocity by target angular velocity
        targAvgVelX = nanmean(dat.TARGET_VELOCITY_X_TARG1(indAllR(i,1):indAllR(i,2)));
        targAvgVelY = nanmean(dat.TARGET_VELOCITY_Y_TARG1(indAllR(i,1):indAllR(i,2)));
        [~, targAvgVel]=cart2pol(targAvgVelX, targAvgVelY);
        eyeAvgVelX = nanmean(dat.RIGHT_VELOCITY_X(indAllR(i,1):indAllR(i,2)));
        eyeAvgVelY = nanmean(dat.RIGHT_VELOCITY_Y(indAllR(i,1):indAllR(i,2)));
        [~, eyeAvgVel]=cart2pol(eyeAvgVelX, eyeAvgVelY);
        velGain = eyeAvgVel/targAvgVel;
        
        %RMSE_Gaze:sqrt of eucledian eye pos - targ pos in degrees ^2
        targPosX = dat.TARGET_X_TARG1(indAllR(i,1):indAllR(i,2));
        targPosY = dat.TARGET_Y_TARG1(indAllR(i,1):indAllR(i,2));
        eyePosX = dat.RIGHT_GAZE_X(indAllR(i,1):indAllR(i,2));
        eyePosY = dat.RIGHT_GAZE_Y(indAllR(i,1):indAllR(i,2));
        diffPosX = 1/nanmean(dat.RESOLUTION_X(indAllR(i,1):indAllR(i,2))).*(eyePosX - targPosX);
        diffPosY = 1/nanmean(dat.RESOLUTION_X(indAllR(i,1):indAllR(i,2))).*(eyePosY - targPosY);
        [~, diffPos]=cart2pol(diffPosX, diffPosY);
        RMSE_Gaze = sqrt(nanmean(diffPos).^2);        
        
        %EDF time
        startEDF = dat.TIMESTAMP(indAllR(i,1));
        endEDF = dat.TIMESTAMP(indAllR(i,2));
        dur = (endEDF - startEDF) + (dat.TIMESTAMP(indAllR(i,1)+1) - dat.TIMESTAMP(indAllR(i,1)));
                
        %save to text file
        sessName = dat.DATA_FILE(indAllR(i,1),:);        
        trialIn = dat.TRIAL_INDEX(indAllR(i,1),:);
        eye='R';
        fprintf(f,'%s\t',sessName);
        fprintf(f,'%d\t',trialIn);
        fprintf(f,'%s\t',eye);
        fprintf(f,'%d\t',indAllR(i,3));
        
        fprintf(f,'%d\t',startEDF);
        fprintf(f,'%d\t',endEDF);
        fprintf(f,'%d\t',dur);        
        
        fprintf(f,'%-10.4f\t',velGain);
        fprintf(f,'%-10.4f\t\n',RMSE_Gaze);
    end
end

fclose(f);
mtitle = ['Output file "results.txt" saved in ', DirName];
[~]=menu(mtitle ,'OK');


