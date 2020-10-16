% PURPOSE
%     Compute element centroids of a specially ordered
%     regular, 2D mesh
% DEPENDENCIES:
% 
% RELATED SCRIPTS:
%     
% Date:
%     Oct-14-2020
%  ----------------------------------------------------------------
clearvars; close all; clc;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath); addpath('functions') 

%% Specify dimmensions:
% Use consistent units in computations:
% MASS	LENGTH	TIME   FORCE   STRESS	ENERGY	 		
%    g	    mm	  ms	   N	  MPa	  N-mm	
%
% steel      steel           
% DENSITY    YOUNG's    GRAVITY   56.33 km/h = 35 mph
% 7.83e-03  2.07e+05  9.806e-03   15.65 mm/ms
%
%      Y
%      ^
%      |                                |
%      |------------------------.       -
%      | 26 | 27 | 28 | 29 | 30 |       |
%      |------------------------|       |
%      | 21 | 22 | 23 | 24 | 25 |       |
%      |------------------------|       |
%      | 16 | 17 | 18 | 19 | 20 |      Ly   (ny)
%      |------------------------|       |
%      | 11 | 12 | 13 | 14 | 15 |       |
%      |------------------------|       |
%      |  6 |  7 |  8 |  9 | 10 |       |
%      |------------------------|       |
%      |  1 |  2 |  3 |  4 |  5 |       |
%      -----------------------------------------> X
%                            (nx * ny)
%                  
%     -/-------- Lx -----------/--
%    (nx - number of elements in X-direction)
%

% Dimmnesions:
Lx = 1000; %mm
Ly = Lx; %mm, same dimmensions
nx = 20;
ny = nx; %same number of subdivsions
n_elems = nx * ny;

% Elements ids are:
element_ids = [1:1:n_elems]';   % we number the element centroids
coordZ = zeros(n_elems,1);  % our tests are 2D so out-of-plane coordinate is zero

% Compute the centroids of the elements first such that we can compute
% the distance between the elements in the next step:
[coordX, coordY] = elemCoordVector(Lx,Ly,nx,ny);

%% Save centroid coordiantes and element_ids to text file:
elem_data = [element_ids, coordX, coordY, coordZ];

% - mopen the file with write permission
propFileName = strcat('elem_centroids','.txt');
fid = fopen(propFileName, 'w');  % open file identifier (handle)
fprintf(fid, ' elem_id      x-coord      y-coord      z-coord\n');
fprintf(fid, '%8d %12g %12g %12g\n', elem_data');  % TRANSPOSE matrix
fclose(fid);   %close file identifier

%% FUNCTIONS:

function [coordX, coordY] = elemCoordVector(Lx,Ly,nx,ny)
    n_elements = nx * ny;
    deltaX = Lx / nx;  % size of element in x-direction
    deltaY = Ly / ny;  % size of element in y-direction
    
    coordX = zeros(n_elements,1); % placeholder for the element coordinates
    coordY = zeros(n_elements,1); % 
    
    %index of the element position in the grid
    ind_gridX =1;  % start at the first element in position (1,1)
    ind_gridY =1;
    
    %   ind_gridY
    %      ^
    %      |                             
    %      |-----------------------------.   
    %      | 6,1 | 6,2 | 6,3 | 6,4 | 6,5 | 
    %      |-----------------------------| 
    %      | 5,1 | 5,2 | 5,3 | 5,4 | 5,5 | 
    %      |-----------------------------| 
    %      | 4,1 | 4,2 | 4,3 | 4,4 | 4,5 | 
    %      |-----------------------------|       
    %      | 3,1 | 3,2 | 3,3 | 3,4 | 3,5 |  
    %      |-----------------------------|  
    %      | 2,1 | 2,2 | 2,3 | 2,4 | 2,5 |       
    %      |-----------------------------|       
    %      | 1,1 | 1,2 | 1,3 | 1,4 | 1,5 |       
    %      ------------------------------------> ind_gridX
    %                            (nx * ny)
    %                  

    for elem_id=1:n_elements  % march element by element:
        % as we move from element to the next element
        coordX(elem_id,1) = deltaX/2 + (ind_gridX -1)* deltaX; % increment x-coordinate
        
        % as we move from the row to the next row of elements
        coordY(elem_id,:) = deltaY/2 + (ind_gridY -1)* deltaY; % increment y-coordinate
    
        ind_gridX = ind_gridX+1;  % increment index x as we progress
        if ind_gridX > nx  % if we exceed the number of subdivisions in x-direction
            % then reset the ind_gridX=1 and increment ind_gridY by +1
            ind_gridX = 1;
            ind_gridY = ind_gridY +1;
        end
    end
end
