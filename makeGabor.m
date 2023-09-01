function [img] = makeGabor(ih,iw,freq,angle,phase,mean,amp,pixPerDeg);

% ih = 512;
% iw = 512; % width
% freq = 4; % spatial frequency
% pixPerDeg = 51;
% angle = 90;
% phase = 0; % phase is in absolute units of the space
% mean = 0.5;
% amp = 0.02;

if pixPerDeg > 0
	freq = freq*iw/pixPerDeg;			% convert from cycles/degree to cycles per stimulus
end
% meshgrid(0:(iw-1),0:(ih-1));
[x y] = meshgrid(linspace(-0.5,0.5,iw),linspace(-0.5,0.5,ih)); % work with range of 1 

phase = phase/freq; % convert phase, relative to frequency of the grating/Gabor
theta = angle*pi/180;
%theta = pi/2; % angle of grating or Gabor filter

sd = 0.2; % standard deviation of Gaussian for specifying Gabor, should be 0.2 or less ideally 

gauss = exp(-(x.^2+y.^2)/(2*sd^2)); % assume circular 2D Gaussian
Gs = sin(2*pi*freq * (x.*sin(theta) + y.*cos(theta) - phase)); % sine grating, if phz = 0
Gc = cos(2*pi*freq * (x.*sin(theta) + y.*cos(theta) - phase)); % cosine grating, if phz = 0

% note, phz refers to phase-shifting the pattern in desired direction, so it is subtracted here 
gaborS = gauss.*Gs; % Gabor filter, odd-symmetric if phz = 0
gaborC = gauss.*Gc; % Gabor filter, even-symmetric if phz = 0

Gc = cos(2*pi*freq * (x.*sin(theta) + y.*cos(theta) - phase)); % cosine grating, if phz = 0
gaborC = gauss.*Gc; % Gabor filter, even-symmetric if phz = 0

img = (gaborC*amp+mean)*255;		% values range from mean-amp to mean+amp, between 0 and 1 for 100% contrast multiply by 255 for display function Screen() 

 

% figure();
% imagesc(img)
% % imagesc(img/max(img,[],'all'))
% % imagesc(gaborC)
% colormap(gray)


% min(gaborC,[],'all')
% max(gaborC,[],'all')
% 
% min(img,[],'all')
% max(img,[],'all')
% 
% min(img/max(img,[],'all'),[],'all')
% max(img/max(img,[],'all'),[],'all')

% title(sprintf('Even-symmetric Gabor for theta =%s',th{1}))

 
%Odd-symmetric

% Gs = sin(2*pi*freq * (x.*sin(theta(i)) + y.*cos(theta(i)) - phase)); % sine grating, if phz = 0
% 
% gaborS = gauss.*Gs; % Gabor filter, odd-symmetric if phz = 0
% 
% imagesc(gaborS)
% 
% colormap(gray)
% 
% title(sprintf('Odd-symmetric Gabor for theta =%s',th{i}))

 

% subplot(4,3,3+(i-1)*3)
% 
% amp = abs(fftshift(fft2(gaborC)));
% 
% imagesc(amp)
% 
% colormap(gray)
% 
% title(sprintf('Amplitude spectrum of the even-symmetric Gabor for theta =%s',th{i}))

