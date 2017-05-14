function randomlinesMATRIX3H(d,sigma_d, a,sigma_a, nxInput, grid, fieldsize, NrFieldsWidth, NrFieldsLength, PHeightChanged)

%randomlinesMATRIX3(d,sigma_d, a,sigma_a)%,max_delta_d, a, sigma_a, h, delta_h, dose,gridx, gridy)
% It generates a large array that corresponds to one row in the
% large pattern and then splits this into individual ccc files. It then
% generates the corresponding CON file with references to all of the CCC
% files.

% To set again later in the lithography control:
pitch = 1;
dose = 0.33;
doseMiddle = 0.04;


%% Create pixel data matrix with randomnes including rules
nxL = (NrFieldsWidth * nxInput) +1; % uneven number to compensate for file loop -1

resolution = grid/(fieldsize*1000);
[VX] = RandomWithRulesInverted(d,sigma_d, a,sigma_a, nxL, grid, resolution, pitch, dose);
% plot(1:nxL,mod(VX(1,:),grid));hold on; plot(1:nxL,mod(VX(3,:),grid));


%% Randomise height (dose)

NrTop = 10; NrBttm = 1; %initialise
NrH = 2*nxL -3; % possible elements to change the eights of (tops and valeys)

while((NrBttm/(NrTop + NrBttm))<0.2); % at least 20% of modifications should happen on top of structures
    IndHeight = unidrnd(100,1,NrH) <= PHeightChanged; % logical derived from random numbers up to 100
    % Delete neighbouring modifications (pick top to delete)
    for p = 1:1:ceil(NrH/2)-1; % for all grating teeth (tops)
        if IndHeight(p*2)==1 && ( IndHeight(p*2 -1)==1 || IndHeight(p*2 +1)==1 );% if neighboruing surfaces are selected for height change,..
            IndHeight(p*2) = 0; %...remove the top one
        end
    end
    NrTop = sum(IndHeight( not(mod(1:NrH,2)) ));
    NrBttm = sum(IndHeight( not(not(mod(1:NrH,2))) ));
end

% incorporate heights
VXwithheights=zeros(6,length(VX(2,:))+NrTop); %initialise

IndTops = not(mod(1:NrH,2)) & IndHeight;
IndBttms = mod(1:NrH,2) & IndHeight;
countVXh = 1; % initialise
countVX = 1; % initialise
y1=0; y2=grid;
for p = 1:1:NrH;
    if IndTops(p);%TOPs
        VXwithheights(1,countVXh) = VX(3,countVX-1)+1; % in pixel
        VXwithheights(2,countVXh) = y1;
        VXwithheights(3,countVXh) = VX(1,countVX)-1; % in pixel
        VXwithheights(4,countVXh) = y2;
        VXwithheights(5,countVXh) = pitch;
        VXwithheights(6,countVXh) = doseMiddle;
        countVXh = countVXh+1;
    elseif IndBttms(p); %BTTMs
        VXwithheights(1:5,countVXh) = VX(1:5,countVX);
        VXwithheights(6,countVXh) = doseMiddle;
        countVX = countVX+1;
        countVXh = countVXh+1;
    elseif mod(p,2); %unchanged
        VXwithheights(:,countVXh) = VX(:,countVX);
        countVX = countVX+1;
        countVXh = countVXh+1;
    end
end


%% Output for ebeam

%%% CCC files
output_counter = 1;
Writestartingbit = 0;    

NrToWrite = length(VXwithheights(2,:));
for i = 1:1:NrToWrite; % for each grating element

%%%-- Possible CCC file openings

    if i == 1; % first run only; Start a fresh CCC file
        fOut= sprintf('%i.ccc',output_counter);
        %tr = strvcat(s1, vv,s3)
        fid = fopen(fOut, 'w');
        %fprintf(fid,'DWSL, num2str(vx)');
        fprintf(fid,'/*--- %i.ccc ---*/\n',output_counter);
        fprintf(fid,['/* CZ' num2str(fieldsize) ',' num2str(grid) ' */\n']);
        fprintf(fid,'PATTERN\n');
    end
    
    if Writestartingbit == 1; % Start a new CCC file; 'w' end of grating
        % element that was not finished in the last block
        fOut= sprintf('%i.ccc',output_counter);
        %tr = strvcat(s1, vv,s3)
        fid = fopen(fOut, 'w');
        %fprintf(fid,'DWSL, num2str(vx)');
        fprintf(fid,'/*--- %i.ccc ---*/\n',output_counter);
        fprintf(fid,['/* CZ' num2str(fieldsize) ',' num2str(grid) ' */\n']);
        fprintf(fid,'PATTERN\n');
        x1 = 0;
        x2 = mod(VXwithheights(3, i-1), grid); % end of previous element
        fprintf(fid,'DWSL(%1.0f,%1.0f,%1.0f,%4.0f,%1.0f,%4.2f);3\n',x1, VXwithheights(2,i), x2, VXwithheights(4,i), VXwithheights(5,i), VXwithheights(6,i));
        Writestartingbit = 0;
    end
    
    if Writestartingbit == 2; % Start a new CCC file after a block has been
        % finished without residues
        fOut= sprintf('%i.ccc',output_counter);
        %tr = strvcat(s1, vv,s3)
        fid = fopen(fOut, 'w');
        %fprintf(fid,'DWSL, num2str(vx)');
        fprintf(fid,'/*--- %i.ccc ---*/\n',output_counter);
        fprintf(fid,['/* CZ' num2str(fieldsize) ',' num2str(grid) ' */\n']);
        fprintf(fid,'PATTERN\n');    
        Writestartingbit = 0;
    end
    
%%%-- Sort each grating element into into individual small grating blocks for each CCC file

    x1 = mod(VXwithheights(1, i), grid); % remainder after clustering into small grating blocks
    x2 = mod(VXwithheights(3, i), grid);
    
    % Add line to CCC file with grating element
    if x2 >= x1; % normally, the end of the grating element has a greater remainder than the beginning ...
        fprintf(fid,'DWSL(%1.0f,%1.0f,%1.0f,%4.0f,%1.0f,%4.2f);3\n',x1, VXwithheights(2,i), x2, VXwithheights(4,i), VXwithheights(5,i), VXwithheights(6,i));
    else % ... but if it already reaches the next small grating block, then...
        x1 = mod(VXwithheights(1, i), grid);
        x2 = grid; %... let it stopp at the end of this block ...
        fprintf(fid,'DWSL(%1.0f,%1.0f,%1.0f,%4.0f,%1.0f,%4.2f);3\n',x1, VXwithheights(2,i), x2, VXwithheights(4,i), VXwithheights(5,i), VXwithheights(6,i));
        fprintf(fid,'!END'); % ... finish the file ...
        fclose(fid);
        output_counter = output_counter + 1;
        Writestartingbit = 1; % ... and starts a new CCC file where the
        % end of the started grating block is written down immediatelly
    end

    if i == (NrToWrite); % at the very end print end loop, finish file
        fprintf(fid,'!END');
        fclose(fid); 
    else % for all others, do nothing, except when ... 
        %... this grating was completely within one block but the next one will start a new block...
        if  mod(VXwithheights(1, i+1), grid) < mod(VXwithheights(3, i), grid); 
            output_counter = output_counter + 1;
            Writestartingbit = 2; % ... then open a new empty file next
        end 
    end

end


%%% CON file
fOut2= sprintf('largearea.con');
fid = fopen(fOut2, 'w');
fprintf(fid,'/*--- largearea.CON ---*/\n');
fprintf(fid,['CZ' num2str(fieldsize) ',' num2str(grid) '\n']);
for j=0:1:NrFieldsLength-1;
    for i=0:1:NrFieldsWidth-1;
        fprintf(fid,'PC%i;\n',i+1);
        fprintf(fid,'%f,%f;\n',i*fieldsize,j*fieldsize);
        fprintf(fid,'PP%i;\n',i+1);
        fprintf(fid ,'%f,%f;\n',i*fieldsize,j*fieldsize);
    end
end
fprintf(fid,'!END');

