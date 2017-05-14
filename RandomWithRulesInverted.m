function [vx2] = RandomWithRulesInverted(d,sigma_d, a,sigma_a, nx, grid, resolution, pitch, dose)
% Create randomnes with border conditions and invert values to write
% valeys of grating into polymer

% For the ebeam output
y1=0; % 0 ~ no offset
y2=grid;

vx=zeros(6,nx+1); % initialise
vx2=zeros(6,nx);

for i = 0:1:nx; % for each period (plus one extra to generate the line needed for the inversion (see below))
    
    % Generate random offset from a fixed lattice. Trunkate gaussian at 2sigma
    [x, width] = MakeRandom(i,d,sigma_d, a,sigma_a, resolution);
    % COORDINATES of structure in pixel
    x1 = round((x-width/2) *resolution); % the first one (negative) here is not used anyways because of the inversion below
    x2 = round((x+width/2) *resolution);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Make sure the lines don't run into each other %%%
    if(i>0.5) %%% not the first time (i=0) since previous element is not existant
    while(vx(3,i)>x1) % if they run into each other, run until they don't bump anymore
        [x, width] = MakeRandom(i,d,sigma_d, a,sigma_a, resolution); % Generate random offset from a fixed lattice. Trunkate gaussian at 2sigma
        % COORDINATES of structure in pixel
        x1 = round((x-width/2) *resolution); x2 = round((x+width/2) *resolution);
    end
    end
    
    if(i==nx) %%% make sure the left edge of the last element does not overrun the last period base
    while( round(((i*d)+(d/4))*resolution) < x1) % while the left edge of the last element overruns the edge
        [x, width] = MakeRandom(i,d,sigma_d, a,sigma_a, resolution); % Generate random offset from a fixed lattice. Trunkate gaussian at 2sigma
        % COORDINATES of structure in pixel
        x1 = round((x-width/2) *resolution); x2 = round((x+width/2) *resolution);
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% create the data matrix column
    vx(1,i+1)=x1;
    vx(2,i+1)=y1;
    vx(3,i+1)=x2;      
    vx(4,i+1)=y2; 
    vx(5,i+1)=pitch;
    vx(6,i+1)=dose;
end


%% invert x1 and x2 in order to write the right boxes and not the
% invers, that is not centered on the grid anymore because of the
% randomnes.
for i=1:1:nx;
    vx2(1,i) = vx(3,i); % take the right line of the first grating line
    vx2(2,i) = vx(2,i);
    vx2(3,i) = vx(1,i+1); % then take the left one of i+1 to invert structure effectively
    vx2(4,i) = vx(4,i);
    vx2(5,i) = vx(5,i);
    vx2(6,i) = vx(6,i);
end