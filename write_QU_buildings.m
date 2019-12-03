function write_QU_buildings(matrix_city)
%this function writes the building input file for Quic-Urb

fid = fopen(fullfile('QU_buildings.inp'),'wt+');
% [matrix_city,PAF,sigma] = Citygenerator;

Nbldg = size(matrix_city,1)/2;      %Number of buildings = matrix size divided by 2
Nnodes = 0;
Version = 6.01;
rough = 0.1;
Geometry = 1;
GroupID = 1;
BldgType = 1;
Atten = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !QUIC 6.01
% 0.1		!Wall roughness length (m)
% 2			!Number of Buildings
% 0			!Number of Polygon Building Nodes

% !Start Building 1
% 1			!Group ID
% 1			!Geometry = Rectangular
% 1			!Building Type = Solid
% 12		!Height [m]
% 0			!Base Height (Zfo) [m]
% 702.5			!Centroid X [m]
% 665			!Centroid Y [m]
% 680			!Xfo [m]
% 665			!Yfo [m]
% 45			!Length [m]
% 50			!Width [m]
% 0			!Rotation [deg]
% !End Building 1

% !Start Building 103
% 49			 !Group ID
% 2			 !Geometry = Elliptical
% 2			 !Building Type = Canopy
% 4.03			 !Attenuation Coefficient
% 4.1			 !Height [m]
% 1.4			 !Base Height (Zfo) [m]
% 233.0			 !Centroid X [m]
% 185.0			 !Centroid Y [m]
% 230.0			 !Xfo [m]
% 185.0			 !Yfo [m]
% 6.0			 !Length [m]
% 6.0			 !Width [m]
% 0			 !Rotation [deg]
% !End Building 103

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'!QUIC %g\n',Version);
fprintf(fid,'%g\t\t\t!Wall roughness length (m)\n', rough);
fprintf(fid,'%g\t\t\t!Number of Buildings\n',Nbldg);
fprintf(fid,'%g\t\t\t!Number of Polygon Building Nodes\n',Nnodes);

j = 1;
for i = 1:size(matrix_city,1)/2;
fprintf(fid,'!Start Building %g\n', i);
fprintf(fid,'%g\t\t\t!Group ID\n',GroupID);
fprintf(fid,'%g\t\t\t!Geometry = Rectangular\n',Geometry);
fprintf(fid,'%g\t\t\t!Building Type = Solid\n',BldgType);
% fprintf(fid,'%g\t\t\t!Attenuation Coefficient\n',Atten);
%Calculations
h = matrix_city(j+1,3);
Zfo = 0;
CentX = (matrix_city(j,1) + matrix_city(j+1,1))/2;
CentY = (matrix_city(j,2) + matrix_city(j+1,2))/2;

Xfo = matrix_city(j,1);
Yfo = CentY;
L = matrix_city(j+1,1) - matrix_city(j,1);
W = matrix_city(j+1,2) - matrix_city(j,2);
Rot = 0;
%
fprintf(fid,'%g\t\t\t!Height [m]\n',h);
fprintf(fid,'%g\t\t\t!Base Height (Zfo) [m]\n',Zfo);
fprintf(fid,'%g\t\t\t!Centroid X [m]\n',CentX);
fprintf(fid,'%g\t\t\t!Centroid Y [m]\n',CentY);
fprintf(fid,'%g\t\t\t!Xfo [m]\n',Xfo);
fprintf(fid,'%g\t\t\t!Yfo [m]\n',Yfo);
fprintf(fid,'%g\t\t\t!Length [m]\n',L);
fprintf(fid,'%g\t\t\t!Width [m]\n',W);
fprintf(fid,'%g\t\t\t!Rotation [deg]\n',Rot);
fprintf(fid,'!End Building %g\n', i);
%
j = j + 2;
end


% fprintf(fid,'!%g\t\t\tPlan Area Fraction\n',PAF);
% fprintf(fid,'!%g\t\t\tSigma\n',sigma);

fclose(fid); %closing fid so it will not be written over


