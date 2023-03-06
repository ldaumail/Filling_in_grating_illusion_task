%% run this on matlab start up
clc; clear; close all;
cd('/Users/tonglab/Desktop/Loic')
addpath(genpath('/Users/tonglab/Desktop/Loic'));

%% run this if debugging
debug = 1;
subject = 'M099';
session = '';
%% Visual phantom task

%% contrast matching control
contrast_matching_control_v3(subject, session, debug) 

contrast_matching_control_v3_laptop(subject, session, debug) 
