% -------------------------------------------------------------------------
% Scan plan.
% 1. Reboot Matlab after all plug-in.
% 2. Check subject, session and debug variables
% 3. Run findVertOffset()
% 4. Run functional scans as needed
% -------------------------------------------------------------------------

%% run this on matlab start up
clc; clear; close all;
cd('/Users/tonglab/Desktop/Loic')
addpath(genpath('/Users/tonglab/Desktop/Loic'));

%% run this if scanning
debug = 0;
subject = 'M015';
session = datestr(now,'YYmmDD');
figSizeDeg = 4;

%% run this if debugging
debug = 1;
subject = 'M099';
session = '';
vertOffset = 0;
horiOffset = 0;
figSizeDeg = 4;

%% Screen offset for phantom, first vertical, then horizontal
[vertOffset,horiOffset] = findOffset_phantom(subject, session, debug);

%% Screen dimensions
%scrSz = findScreenSize_LD(subject, session, debug);

%% Localizer task
%alternate check 12s - fixation 12s -- (150 TR, 5:00)
LD_localizer_horiOffset(subject, session, horiOffset, vertOffset, figSizeDeg, debug);

%% phantom task
%%Wait for backtick before every new block/ color discrimination  (147 TR, 4:54)
````````````sca
sca

%%backtick at exp start only/ color discrimination (147 TR, 4:54)
%LD_phantom_mri_V8(subject, session, vertOffset, debug, figSizeDeg);

%%backtick at exp start only/ color detection (147 TR, 4:54)
%LD_phantom_mri_V7(subject, session, vertOffset, debug, figSizeDeg);

%% PRF (150 TR, 5:00) multibar
% stim = 'Multibar';
% kayPRF_dotTask(subject, session, vertOffset, debug, stim)

%% Screen offset for pRF mapping wedgering, first vertical, then horizontal
[vertOffset,horiOffset] = findOffset_wedgering(subject, session, debug);

%% PRF (150 TR, 5:00) wedgering
stim = 'WedgeRing';
kayPRF_dotTask_horiOffset(subject, session, horiOffset, vertOffset, debug, stim)

%% PRF (150 TR, 5:00)
% rot = 'CW';
% kayPRF_dotTask_wedgeOnly(subject, session, vertOffset, debug, rot)
% 
%% Resting state 
restingLength = 120; % length of resting state in seconds
restingState_horiOffset(restingLength, horiOffset, vertOffset, debug);
