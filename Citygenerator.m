function Citygenerator(mu_inp, sigma_inp, PAF_inp)

% mu_inp = 10;
% sigma_inp = 5;
% PAF_inp = 0.1;

%% City specifications
%%%%% Domain and block sizes %%%%%%%
worldDim.x = 1000; %210;
worldDim.y = 1000; %110
block.x = 100;
block.y = 50;
street_width = 7;
%%%%% Number of Blocks %%%%
Nblocks_x = 1 + floor((worldDim.x - block.x)/(block.x + street_width));    %Number of blocks on the x direction
Nblocks_y = 1 + floor((worldDim.y - block.y)/(block.y + street_width));    %Number of blocks on the y direction
Nblocks = Nblocks_x*Nblocks_y;                                             %total number of blocks
%%%% Recalculation of the domain size without edges %%%%%%
worldDim.x = (Nblocks_x - 1)*(block.x + street_width) + block.x;           %To cut put domain edges without buildings
worldDim.y = (Nblocks_y - 1)*(block.y + street_width) + block.y;
%%%% Writing QU_simparams.inp %%%%%%
worldDim.z = 100;
d.x = 5;
d.y = 5;
d.z = 5;
worldDim.x = d.x*(ceil(worldDim.x/d.x));
worldDim.y = d.y*(ceil(worldDim.y/d.y));
worldDim.z = d.z*(ceil(worldDim.z/d.z));
write_QU_simparams(worldDim, d);                                           %Write .inp file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_city = worldDim.x*worldDim.y;                                            %New domain without edges

%% Parameters
% PAF = 0.1;                                                                 %Plan Area Fraction = A_roof/A_city
PAF = PAF_inp/100;
mu_DT = 0;                                                                %Mean value and standard deviation
sigma_DT = 0;                                                             %for downtown and outskirt region
DTfrac = 0;                                                              %Fraction of Downtown [0 1]
% mu_OS = 10;
% sigma_OS = 5;
mu_OS = int32(mu_inp);
sigma_OS = int32(sigma_inp);
A_DT = DTfrac^2 * A_city;                                                  %DT area
A_OS = (1-DTfrac^2) * A_city;                                              %OS area
if mu_DT == 0 || sigma_DT == 0
    mu_ave = mu_OS;                                
    sigma_ave = sigma_OS; 
else
    mu_ave = (mu_DT*A_DT + mu_OS*A_OS) / A_city;                                
    sigma_ave = (sigma_DT*A_DT + sigma_OS*A_OS) / A_city;
end
A_variation = 0.15;                                                        %Variation of rooftop area in the "block._fun"

%%%% Rooftop Area %%%%%
A_roof = PAF*A_city;                                                       %Total rooftop area
A_roof_block_in = A_roof/Nblocks;                                          %A_roof per block
%%%% Domain Radius %%%%%
domainRad = sqrt(worldDim.x^2 + worldDim.y^2)/2;

%%%%%%%%%%%%%%%%%%%%%%
%Matric_city structure
%%%%%%%%%%%%%%%%%%%%%%
%      matrix_city
%       X   Y   Z
%bldg1  x1  y1  z1
%       x2  y2  z2
%bldg2  x1  y1  z1
%       x2  y2  z2
%bldg3  x1  y1  z1
%       x2  y2  z2
%  .    .   .    .
%  .    .   .    . 
%  .    .   .    .
%%%%%%%%%%%%%%%%%%%

%% Build the city
%%%% Parametrs inizializations %%%%%
g = 1;
k = 1;
step_block.y = 0;
matrix_city = zeros(3);
A_roof_tot = 0;
A_roof_block_out = 0;
%%%%%% Block vectors %%%%%%%
n_dt = [1 3 4 6 7 8];         %Block types in downtown
n_os = [1 2 3 4 5 6 7 8];       %Block types in outskirt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  While loop for Y  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while step_block.y <= (worldDim.y - block.y)
    step_block.x = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  While loop for X  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while step_block.x <= (worldDim.x - block.x)
            A_roof = A_roof - double(A_roof_block_out) ;                   %Total rooftop area left to create buildings
            A_roof_block_in = A_roof/Nblocks ;                             %Rooftop area left per block
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %When it comes to the last block, we want to use the
            %remaining rooftop area, without randomizing it
            if Nblocks == 1;
                flag = 1;                                                  %This counter is used to identify the last block in the city
            else
                flag = 0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculation of the radial distance from the domain center, X
            % distance and Y distance
            DTRad(g) = sqrt((worldDim.y/2 - step_block.y - block.y/2)^2 + (worldDim.x/2 - step_block.x - block.x/2)^2); 
            DTcoord.x = abs(worldDim.x/2 - step_block.x);
            DTcoord.y = abs(worldDim.y/2 - step_block.y);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If the block is out of the DT fraction area, outskirt
            % parameters are used, otherwise DT parameters are employed
%             if DTcoord.x/(worldDim.x/2) <= DTfrac && DTcoord.y/(worldDim.y/2) <= DTfrac
%                 block_type = randsample(n_dt,1);
%                 mu = mu_DT
%                 sigma = mu_DT
%             else
                block_type = randsample(n_os,1);
%                 mu = mu_OS;
%                 sigma = mu_OS;
%             end
%                         
%                    block_type = 8;%randi([1 4]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now the buildings for picked block type are displaced in the block
            %with the "block_fun" function. Then the coordinates of each bldng are added to the
            %city matrix.
            %If any building height is lower or higher than the limits, the
            %block is discarded and "block_fun" is called again
            h_check = 1;
%             while h_check < 2 %h_check > 100

                [matrix_block, block, A_roof_block_out, h_check] = block_fun(block, block_type, A_roof_block_in, flag, mu_OS, sigma_OS, A_variation); 
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %The block matrix from "block_fun" is read and the buildings
            %are displaced according to previus buildings in the city
                for j=1:size(matrix_block,1)
                    %X
                    matrix_city(k,1) = step_block.x + matrix_block(j,1);
                    %Y
                    matrix_city(k,2) = step_block.y + matrix_block(j,2);
                    %Z
                    matrix_city(k,3) = matrix_block(j,3);
                    k = k+1; 
                end
                g = g + 1;
                A_roof_tot = A_roof_tot + double(A_roof_block_out);        %Displaced rooftop area increases at every step
                Nblocks = Nblocks - 1;                                     %The number of available blocks to displace decreases at each step
                step_block.x = step_block.x + block.x + street_width;      %This variable is used to move the X coordinate to the position of the coming block
        end  %end while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step_block.y = step_block.y + block.y + street_width;                      %to pass to the following line of blocks on Y
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_QU_buildings(matrix_city);                                           %Write .inp file

% 
% 
% 
h_vector = matrix_city(2:2:end,3);
% % A_city = A_city;
% % A_roof_tot = A_roof_tot;
PAF2 = A_roof_tot/A_city                                                  %Plan area fraction
% % matrix_city;
% % size(matrix_city,1);
% aveh = sum(matrix_city(:,3))/(size(matrix_city,1)/2)
aveh = mean(h_vector);
% size(matrix_city,1)/2
% size(h_vector)
% h_vector';
stdv_h = std(h_vector);

exit
% toc
