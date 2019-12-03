function [matrix_block,block, A_roof, h_check] = block_fun(block, type, A_roof_block, flag, mu, sigma, A_variation)
% "Oblong blocks range considerably in width and length.
% The standard block in Manhattan is about 264 by 900 feet (80 m × 274 m);
% and in some U.S. cities standard blocks are as wide as 660 feet (200 m).
% The blocks in Edmonton, Canada are 330 by 560 feet (100 m × 170 m)."
%Wikipedia


% bld.x = randi([5 30]);
% bld.y = randi([5 30]);
% bld.z = randi([20 80]);
%
% block.x = randi([190 220]);;
% block.y = randi([80 100]);;

% block.x = 200;
% block.y = 100;
% flag
% A_roof_block
% type

h_min = 1;
h_max = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%When it comes to the last building, we want to use the remaining rooftop
%area, without randomizing it. Otherwise the rooftop area per block is
%randomized not to have each block with the same area
if flag == 1;
    A_roof = A_roof_block;
else
    A_roof = (1-A_variation)*A_roof_block + rand(1)*(A_variation*2)*A_roof_block;                 % +- 15% A_roof_block
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_A = A_roof/(block.x*block.y);                                            %New PAF, if A_variation = 0, PAF = A_A
flowerbed = randi([3 6]);                                                  %Size of flowerbed
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
A_block = 0;                                                               %Inizialization of the rooftop area per block
%% Block type 1
%1 single building per block. Typically skyscrapers
if type == 1;
    %X
    matrix_block(1,1) =  double(block.x/2 - block.x/2*sqrt(A_A));
    matrix_block(2,1) =  double(block.x/2 + block.x/2*sqrt(A_A));
    %Y
    matrix_block(1,2) =  double(block.y/2 - block.y/2*sqrt(A_A));
    matrix_block(2,2) =  double(block.y/2 + block.y/2*sqrt(A_A));
    %Z
    matrix_block(1,3) =  0;
    matrix_block(2,3) =  height;
end

%% Block type 2
%8 single family houses --> villas
if type == 2;
    k=1;
    y_step = 0;
    for j=1:4
        x_step = 0;
           for i=1:4
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
%             height = normrnd(mu, sigma);
            %X
            matrix_block(k,1) = x_step + block.x/8*(1 - sqrt(A_A));
            matrix_block(k+1,1) = x_step + block.x/8*(1 + sqrt(A_A));
            %Y
            matrix_block(k,2) = y_step + block.y/8*(1 - sqrt(A_A));
            matrix_block(k+1,2) = y_step + block.y/8*(1 + sqrt(A_A));
            %Z
            matrix_block(k,3) = 0;
            matrix_block(k+1,3) = height;
            %
            k=k+2;
            x_step = x_step + block.x/4;
            end  %end For
            y_step = y_step + block.y/4;
    end  %end While
end

%% Block type 3
%single large building with courtyard in the middle
if type == 3;

    %Bld 1 (west)
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
%     height = normrnd(mu, sigma);
    matrix_block(1,1) = block.x/4 - (A_roof/4/(block.y - 2*flowerbed));
    matrix_block(1,2) = flowerbed;
    matrix_block(1,3) = 0;

    matrix_block(2,1) = block.x/4;
    matrix_block(2,2) = block.y - flowerbed;
    matrix_block(2,3) = height;
    
    %Bld 2 (North)
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
%     height = normrnd(mu, sigma);
    matrix_block(3,1) = block.x/4;
    matrix_block(3,2) = block.y - flowerbed - A_roof/4/(block.x/2);
    matrix_block(3,3) = 0;

    matrix_block(4,1) = block.x*3/4;
    matrix_block(4,2) = block.y - flowerbed;
    matrix_block(4,3) = height;
    
    %Bld 3 (South)
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
%     height = normrnd(mu, sigma);
    matrix_block(5,1) = block.x/4;
    matrix_block(5,2) = flowerbed;
    matrix_block(5,3) = 0;

    matrix_block(6,1) = block.x*3/4;
    matrix_block(6,2) = flowerbed + A_roof/4/(block.x/2);
    matrix_block(6,3) = height;

    %Bld 4 (East)
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
%     height = normrnd(mu, sigma);
    matrix_block(7,1) = block.x*3/4;
    matrix_block(7,2) = flowerbed;
    matrix_block(7,3) = 0;

    matrix_block(8,1) = block.x*3/4  + (A_roof/4/(block.y - 2*flowerbed));
    matrix_block(8,2) = block.y - flowerbed;
    matrix_block(8,3) = height;
end

%% Block type 4
%1 single building per block. Typically deparment stores
if type == 4;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %X
    matrix_block(1,1) =  flowerbed;
    matrix_block(2,1) =  double(flowerbed + A_roof/(block.y - 2*flowerbed));
    %Y
    matrix_block(1,2) =  flowerbed;
    matrix_block(2,2) =  block.y - flowerbed;
    %Z
    matrix_block(1,3) =  0;
    matrix_block(2,3) =  height;
end

%% Block type 5
%2 lines of random stores
if type == 5;
    %Width
    w.min = int32(0.05*block.x);
    if w.min < 5
        w.min = 5;
    end
    w.max = int32(0.15*block.x);
    if w.max > 20
        w.max = 20;
    end
    width = randi([w.min w.max]);
    %Length
    l.min = int32(0.35*block.y - flowerbed);
    l.max = int32(0.45*block.y - flowerbed);
    length = randi([l.min l.max]);
    edge.y = flowerbed;
    k = 2;
    n = 1;
    matrix_block(1,1) = width;
    matrix_block(1,2) = edge.y;
    matrix_block(1,3) = 0;
    check = 0;
    for n = 1:2;
        i = 2;
        X(1) = width;
        while matrix_block(k-1,1) < (block.x-w.max) 
            %Random sizes
            width = randi([w.min w.max]);
            length = randi([l.min l.max]);
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
%             height = normrnd(mu, sigma);
            %X
            X(i) = X(i-1) + width;
            X(i+1) = X(i);

            %Bld1
            matrix_block(k,1) = X(i);
            matrix_block(k,2) = length + edge.y;
            matrix_block(k,3) = height;

            A_bldg = (matrix_block(k,1) - matrix_block(k-1,1))*(matrix_block(k,2) - matrix_block(k-1,2));
            A_roof = A_roof - A_bldg;
            
            if A_roof < 10
               check = 1;
               break
            else       
                %Bld2
                matrix_block(k+1,1) = matrix_block(k,1);
                matrix_block(k+1,2) = edge.y;
                matrix_block(k+1,3) = 0;
                %
                k = k + 2;
                i = i + 2;
            end
        end %end while loop
        
        if check == 1;
            break
        else        
            edge.y = 0.55*block.y;
            matrix_block(k-1,1) = width;
            matrix_block(k-1,2) = edge.y;
            matrix_block(k-1,3) = 0;
            n = 2;
        end
    end %end for loop
    
        if check == 1;
        else
        matrix_block = matrix_block(1:end-1,:);
        end
end    %First "if" for the building type   

%% Block type 6
%2 buildings with alley in between
if type == 6;
    %Bldg 1
    %X
    matrix_block(1,1) =  flowerbed;
    matrix_block(2,1) =  double(flowerbed + A_roof/(2*(block.y - 2*flowerbed)));
    %Y
    matrix_block(1,2) =  flowerbed;
    matrix_block(2,2) =  block.y - flowerbed;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %Z
    matrix_block(1,3) =  0;
    matrix_block(2,3) =  height;
    
    %Bldg 2
    %X
    matrix_block(3,1) =  double(block.x - flowerbed - A_roof/(2*(block.y - 2*flowerbed)));
    matrix_block(4,1) =  block.x - flowerbed;
    %Y
    matrix_block(3,2) =  flowerbed;
    matrix_block(4,2) =  block.y - flowerbed;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %Z
    matrix_block(3,3) =  0;
    matrix_block(4,3) =  height;
end

%% Block type 7
%2 buildings attached with no alley in between
if type == 7;
    %Bldg 1
    %X
    matrix_block(1,1) =  flowerbed;
    matrix_block(2,1) =  double(flowerbed + A_roof/(2*(block.y - 2*flowerbed)));
    %Y
    matrix_block(1,2) =  flowerbed;
    matrix_block(2,2) =  block.y - flowerbed;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %Z
    matrix_block(1,3) =  0;
    matrix_block(2,3) =  height;
    
    %Bldg 2
    %X
    matrix_block(3,1) =  matrix_block(2,1);
    matrix_block(4,1) =  double(matrix_block(2,1) + A_roof/(2*(block.y - 2*flowerbed)));
    %Y
    matrix_block(3,2) =  flowerbed;
    matrix_block(4,2) =  block.y - flowerbed;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %Z
    matrix_block(3,3) =  0;
    matrix_block(4,3) =  height;
end

%% Block type 8
%2 buildings attached with no alley in between
if type == 8;
    %Bldg 1
    %X
    matrix_block(1,1) =  flowerbed;
    matrix_block(2,1) =  double(flowerbed + A_roof/(2*(block.y - 2*flowerbed)));
    %Y
    matrix_block(1,2) =  flowerbed;
    matrix_block(2,2) =  block.y - flowerbed;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %Z
    matrix_block(1,3) =  0;
    matrix_block(2,3) =  height;
    
    %Bldg 2
    %X
    matrix_block(3,1) =  double(block.x - flowerbed - A_roof/(2*(block.y - 4*flowerbed)));
    matrix_block(4,1) =  block.x - flowerbed;
    %Y
    matrix_block(3,2) =  2*flowerbed;
    matrix_block(4,2) =  block.y - 2*flowerbed;
height = 0;
while (height < h_min) || (height > h_max)
    height = 5 + normrnd(mu, sigma);                                               %Height picked according to the corresponding normal distribution
end
    %Z
    matrix_block(3,3) =  0;
    matrix_block(4,3) =  height;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix_block = int32(matrix_block);                                        %Round any coordinate to the closer integer

%Check if any height is out of the limits
h_check = 0;
% for i = 2:2:size(matrix_block,1)
%     if matrix_block(i,3) < 5 || matrix_block(i,3) > 120
%         h_check = h_check + 1;
%     else 
%     end
% end

%Calculate the rooftop area
A_roof = 0;
p = 1;
for i = 1:size(matrix_block)/2
    A_roof = (matrix_block(p+1,1) - matrix_block(p,1))*(matrix_block(p+1,2) - matrix_block(p,2)) + A_roof;
    p = p+2;
end

