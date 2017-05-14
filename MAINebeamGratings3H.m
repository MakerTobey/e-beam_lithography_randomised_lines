%%% MAINebeamGratings3H
% Run this script to generate ebeam-lithography pattern files and to keep
% a convenient overview over your settings.

%% Inputs:
% individual fieldsize in mm (edgelength of square):
fieldsize = 0.6;
% Nr of total pixels
grid = 20000;
% period length of structure in um
period = 1.3;
% width of grating teeth in um
width = 0.73;
% standard deviation of structure element displacements
sigma_d = 0.29*2; % disorder Hibiscus
% standard deviation of grating teeth width
sigma_a = 0.16*2; % disorder Hibiscus
% percentage of lines to be height modified, keep below 50
PHeightChanged = 40;% 10; 20; 40;

%%% Matrix settings:
% total field width in mm (edgelength of square):
NrFieldsWidth = 33; % 33*0.6 = 19.8mm
% if it should not be a square, change this value (dose trial etc.):
NrFieldsLength = NrFieldsWidth;%6;
% small stripes (optical testing) 4 x 6

%% Generate a grating that can be used in the ebeam software
% Nr of periods that fit in the grating
nx = floor(fieldsize*1000/period); % disp('Nr of periods that fit in the grating:'); disp(nx);
% period of structure
d = period;%round(period*resolution);
% width of grating teeth
a = width;%round(width*resolution);

disp('resolution: size of pixel in nm'); disp(round((fieldsize*1000000)/grid));

%%% Run and write whole grating matrix
randomlinesMATRIX3H(d,sigma_d, a,sigma_a, nx, grid, fieldsize, NrFieldsWidth, NrFieldsLength, PHeightChanged);

