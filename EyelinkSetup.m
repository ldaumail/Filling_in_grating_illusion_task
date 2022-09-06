function EyelinkSetup(ini)
%Basic Eyelink Eyetracker setup. Assumes that the eyetracker is
%connected, and instantiates the eyelink object as 'el'.
%Modified from IdenLocGabor 01/22/2015

global el w rect xc yc EyeData

if ini == 0
tmp = EyelinkInit(0);
end
el = EyelinkInitDefaults(w);

el.backgroundcolour = 255;
el.foregroundcolour = 255;  %what's the differenece btwn foregroundcolour and calibtargetcolour?
el.calibrationtargetcolour= 0;
el.msgfontcolour  = 0;
el.targetbeep = 0;

EyelinkUpdateDefaults(el);

Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0,0,rect(3)-1, rect(4)-1);
Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0,0,rect(3)-1, rect(4)-1);

Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA'); 


EyeData.edfFile=sprintf('demo.edf'); % overwrite the edf each time. or it would full the hard drive of eyelink computer
Eyelink('Openfile', EyeData.edfFile);

% Calibrate the eye tracker
EyelinkDoTrackerSetup(el);

% do a final check of calibration using driftcorrection
% EyelinkDoDriftCorrection(el);

WaitSecs(0.1);
%begin_record_time = GetSecs();
Eyelink('StartRecording');

el.eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
if el.eye_used == el.BINOCULAR % if both eyes are tracked
    el.eye_used = el.LEFT_EYE; % use left eye
end

WaitSetMouse(xc,yc,0); % set cursor and wait for it to take effect

HideCursor;

% Wait until all keys on keyboard are released:
while KbCheck; WaitSecs(0.1); end

end
