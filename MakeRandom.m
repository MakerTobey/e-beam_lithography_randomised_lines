function [x, width] = MakeRandom(i,d,sigma_d, a,sigma_a, resolution)
% Actual randomisation of grating lines in period and width excluding
% any border conditions

    % POSITION
    x=0; %initialize position variable   
    % randomize (total) position:
    while x > i*d+2*sigma_d || x <i*d-2*sigma_d; %repeat if x is not in range of 2*sigma around the mean
        nd= random('Normal',0,sigma_d); % Normal Distribution, mean, standard deviation
        x = i*d+nd; %create spatial variable with random shift
    end
   
    % WIDTH
    width=a + 3*sigma_a; %initialize
    if sigma_a==0 %if no width variation is introduced
        width =a;
    else
    % here we randomise the width of the each element. Again trunkated at 2sigma 
    %  and in the boundaries of the period length and the resolution limit
    while width > a + 2*sigma_a || width < a-2*sigma_a || width < (1/resolution) || width > d;
        width= random('Normal',a,sigma_a);
    end
    end