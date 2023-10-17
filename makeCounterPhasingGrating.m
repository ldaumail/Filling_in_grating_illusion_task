function [img] = makeCounterPhasingGrating(ih,iw,spfreq,angle,phase,mean,amp,tmpPhase,pixPerDeg,r1,r2,bcolor);

% usage: [img] = makeSineGrating(ih,iw,freq,angle,phase,mean,amp,pixPerDeg,r1,r2,bcolor)
%

% Function makes sine-wave grating and presents it either in a rectangle, circle or annulus aperature
%
% ih = image height, iw = image width (in pixels)
% freq = grating frequency in either: cycles/degree (if pixPerDeg specified) or cyc/stimulus width (iw)
% angle = angle of grating, 0=horizontal, 90=vertical
% mean = mean intensity level
% amp = amplitude of sinewave
% r1 and r2 = inner and outer radii of annulus (r1=0 makes circle), program makes rectangle if not specified
%
% created by Frank Tong on 2000/01/10

%disp('In makeSineGrating')
if nargin < 12; bcolor = mean; end		% bcolor = background color, default value = mean intensity
	
if nargin < 9; pixPerDeg = 0; end		% set pixPerDeg to 0 to specify freq in cycles/stimulus
if pixPerDeg > 0
	spfreq = spfreq*iw/pixPerDeg;			% convert from cycles/degree to cycles per stimulus
end

% ih = 512;
% iw = 512; % width
% freq = 4; % spatial frequency
% pixPerDeg = 51;
% angle = 90;
% phase = 0; % phase is in absolute units of the space
% mean = 0.5;
% amp = 0.02;
% r1 =5;
% r2 =6;
% bcolor = experiment.backgroundColor(1);
%  
[X,Y] = meshgrid(0:(iw-1),0:(ih-1));	% specify range of meshgrid
%[X,Y] = meshgrid(linspace(-0.5,0.5, iw),linspace(-0.5,0.5,ih));
img = sin(spfreq*2*pi/iw*(X.*sin(angle*2*pi/360)+Y.*cos(angle*2*pi/360))-phase*2*pi/360)*sin(2*pi-tmpPhase*2*pi/360);% make sine wave, values range from -1 and 1
img = (img*amp+mean)*255;		% values range from mean-amp to mean+amp, between 0 and 1 for 100% contrast multiply by 255 for display function Screen()
%img(img>255) = 255;
%img(img<0) = 0;

if mean > 1	|| amp >1					% likely for ranges 0-255 rather than 0-1
	img = round(img);					% round img values
end


if nargin > 9							% make circle or annulus aperature
	r=sqrt((X-iw/2).^2+(Y-ih/2).^2); 	% calculate eccentricity of each point in grid relative to center
	img(r>r2 | r<r1) = bcolor;
end
