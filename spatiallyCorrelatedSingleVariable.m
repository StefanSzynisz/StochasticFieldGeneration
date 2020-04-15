% PURPOSE
%     Generate stochastic, correlated variables 
%     two variables are generated but the principle works 
%     for any number of variables
% DEPENDENCIES:
% 
% RELATED SCRIPTS:
%     steelfoam_rand_props.m
% Date:
%     March-31-2020
%  ----------------------------------------------------------------
clearvars; close all; clc;
[thisPath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(thisPath);
%% Specify key parameters for FE solver function:
% Specify dimmensions:
% Use consistent units in computations:
% MASS	LENGTH	TIME   FORCE   STRESS	ENERGY	 		
%    g	    mm	  ms	   N	  MPa	  N-mm	
%
% steel      steel           
% DENSITY    YOUNG's    GRAVITY   56.33 km/h = 35 mph
% 7.83e-03  2.07e+05  9.806e-03   15.65 mm/ms
%
% Statistical input:
% -------------------
variableMEAN = 0; %mean of the random field
variableSTD = 1.0 ;  %standard deviation of the random field

gamma = 150;  %Lx/nx;  %mm, spatial correl ation length
%==========

calc_correl_flag = 1;  % if 1 then the correlation matrix is computed and saved,
% if 0 then it is loaded from previuosly computed .mat variable
infile = 'correlationMatrix';  %name of file with the variables

% 'cholesky' is the fastest methods BUT it seems to work/succeed for
% small correlation length(s) only, 'eigen' is second fasters, and 
% 'svd' single value decomp (SVD) the slowest method.
decomposition_flag = 'cholesky';  % 'cholesky','eigen','svd'
%-------------------------------

precision_flag = 'single'; %'double' or 'single
plot_flag = 1;   % =1 means plot, =0 do not plot. 
sparse_flag = 0; %=1 means sprase matrix is used for Cholesky decomposition

% In case of memory issues see this website:
% https://uk.mathworks.com/help/matlab/matlab_prog/resolving-out-of-memory-errors.html 
% gives some tips for reducing memory, in particular the bit at the bottom about 
% the swap memory size might be interesting. 
% Also, you need to remove the array size limit in Matlab preferences.

%% Mesh and sample geometry:

Lx = 1000; %mm 
Ly = 1000; %mm
nx = 50;  % 250 x 250 (=62,500) elements are OK.
ny = 50;  % 40^3 = 64,000
n_elements = nx * ny;
%-------------------

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

%% Compute the correlation matrix, K
if calc_correl_flag == 1
    
    % 1) Compute the centroids of the elements first such that we can compute
    % the distance between the elements in the next step:
    coordCentroids = elemCoordVector(Lx,Ly,nx,ny,precision_flag);

    % 2) using centroid coordinates of each element, compute the correlation 
    % matrix based on distance between each element and spatial correlation
    % parameter, gamma:
    if sparse_flag==1
        K = correlMatrix(coordCentroids, gamma,'sparse',precision_flag);  % 'sparse' matrix is generated
        figure(1);
        spy(K); % show the sparsity of the correlation matrix
        title('Sparsity of K');
    else
        K = correlMatrix(coordCentroids, gamma,'full',precision_flag);  % 'full' matrix is generated
        disp('Correlation matrix, K was succesfully assembled.');
    end
    
    % Eigenanalysis is the only working option for sparse matrices.
    if sparse_flag == 1 
          
        % Eigenvalues:
        disp('Eigenvalue diagonalization for the sparse matrices was used because it is the only option with sparse matrices.');
        [Veig,D] = eigs(K,n_elements);
        B = Veig * sqrt( abs(D) );
        disp('Eigenvector decomposition completed.'); 
        
    else %for non-sparse matrices:
        
        % Go over the decomposition options for 'full' matrices
        if strcmp(decomposition_flag, 'cholesky')
            disp('Cholesky decomposition has started.');
            [Btransp,chol_flag] = chol(real(K));
            % flag = 1; % simulate situation when cholesky decomposition failed
            if ~chol_flag % If flag = 0 then the input matrix is symmetric positive definite 
               %and the factorization was successful
                B = transpose(Btransp);
                clear Btransp K; % clear variables to free memory
                disp('Cholesky decomposition was succesful!');
            else % Cholesky was not succesfull so use SVD decomposition
                % which takes more time
                disp('Cholesky decomposition did not succeed. We can try eigen decomposition instead.');
                decomposition_flag = 'eigen'; % change the decomposition flag
            end
        end
        
        if strcmp(decomposition_flag, 'eigen')
            disp('Eigenvalue diagonalization has initiated.');
            if strcmp(precision_flag,'single')
                [Veig,D,~] = eig(K); % works on both single and double precision.
            elseif strcmp(precision_flag,'double')
                [Veig,D] = eigs(K,n_elements); %numerical eigenanalysis
                % works only on double precision
            end
            clear K;  % To release memory
            B = Veig * sqrt( abs(D) );
            clear Veig D; % again clear variables to release memory
            disp('Eigenvector diagonalization was succesful!'); 
        end
                
        if strcmp(decomposition_flag, 'svd')
            disp('Singe Value Decomposition (SVD) was used.');
            [U, S, ~] = svd(K);   % SVD decomposition 
            % Matrix with the correlation structure:
            B = U * sqrt(S);
            disp('SVD was succesful!')
        end
    
    end

    %save correlation matrix if instructed to do so.
    save(strcat(infile,'_',decomposition_flag),'B','coordCentroids','-v7.3');
    
elseif calc_correl_flag == 0
    %load correlation matrices from a file named infile_svd if
    %instructed to do so
    load(strcat(infile,'_svd'));
end

clear Btransp

%% Spatially correlated random variable:
n_elements = nx * ny;
phi_uncorrel = rand(n_elements,1); % Generate vector of uncorelated random numbers
phi_correl = B *phi_uncorrel;
timestamp = datestr(now,'HHMMSS_FFF'); % get_miliseconds to create unique realizations

if sparse_flag ==1
    figure(2);
    spy(B);  % Show the sparsity of B.
    title('Sparsity of B');
end

clear phi_uncorrel
%% Generate random variable:
variableRand = phi_correl*variableSTD + ones(n_elements,1)*variableMEAN;
min_Var = min(variableRand); max_Var = max(variableRand);

%% Plot Random variable to asses its values:
% ****************************************
resolution = 1.0 * min(Lx/nx, Ly/ny);  %mm
saveDataPath = strcat(thisPath,'\data\',timestamp,'\');
mkdir(saveDataPath);
%-------------------------------------

% random Variable:
if plot_flag == 1
    plotTitle = strcat('Random Variable--Correl-lenght--',' ',num2str(gamma),'mm');
    figure(3); left_pos = 3; bott_pos = 12; %cm
    % double() precisin for plotting:
    surfPlot(double(coordCentroids(:,1)),double(coordCentroids(:,2)),double(variableRand), plotTitle,'', ...
        saveDataPath,double(resolution),left_pos,bott_pos, ...
        double(min_Var), double(max_Var), 'colorbar_on');
end

%% Save the random variables into .mat file;
dataFileName = strcat('data\','randomVariable_', timestamp, '.mat');
save(dataFileName,'variableRand', '-v7.3');
% Version -v7.3 or higher will save your variable ( > 2GB ).

%% FUNCTIONS:
%============

function matrixK = correlMatrix(coordCentroids, gamma,sparse_flag,precision_flag)
    n_elements = size(coordCentroids,1);
    % Correlation matrix (Turn-off when not enough memory):
    if strcmp(sparse_flag,'sparse') % if the flag is on, then use sparse matrix
        matrixK = sparse(n_elements,n_elements);
    else
        if strcmp(precision_flag,'single')
            matrixK = zeros(n_elements,n_elements,'single'); % use full matrix
        elseif strcmp(precision_flag,'double')
            matrixK = zeros(n_elements,n_elements); % use full matrix
        end
    end
    
    h1 = waitbar(0,'Correlation matrix | Overal progress');
    set(h1, 'Units', 'centimeters', 'Position', [15 15 10 2]);
    % it is a symmetric matric, so we can add up things only over halfo the
    % matrix, and using transpose to get the full matrix:
    for i=1:n_elements-1
        waitbar(i/n_elements,h1);  % show progress of the correlation matrix generation
        %--------------
        skip_rows = 1;
        %-------------
        if strcmp(sparse_flag,'sparse') && rem(i, skip_rows)==0
            h2 = waitbar(0,sprintf('Progress for row %d out of %d.',i,n_elements));
        end
        for j=i+1:n_elements
            % Plot waitbar every waitNrows:
            if strcmp(sparse_flag,'sparse') && rem(i, skip_rows)==0
                if rem(j,n_elements/20) ==0
                    waitbar(j/n_elements,h2);
                end
            end
        
            % distance between the element centorids:
            dist_ij = norm(coordCentroids(i,:)-coordCentroids(j,:));
            
            %calculate the correlation matrix to be used in generating the correlated
            %random field.  This matrix defines the correlation coefficient of the
            %value of the material property field at every pair of element centroids.
            K_ij = exp( -dist_ij^2 / gamma^2 );
             
             % round small numbers to zero
            tolerance = 1e-4;
            if (K_ij>tolerance)
                matrixK(i,j) = K_ij;  % filter out small correlations.
            else
                % just leave it blank, which means zero.
            end
        end
        
        if strcmp(sparse_flag,'sparse') && rem(i, skip_rows)==0
            close(h2);
        end
    end
    close(h1); % close the waitbar
    
    % Add the missing other half of the matrix, and ones on the
    % diagonal:
    TOL = 1e-06; % add small tolerance to the diagonal to improve convergence of cholesky decomposition
    matrixK = matrixK + matrixK' + eye(n_elements)*(1+TOL);
end


function [coordCentroids] = elemCoordVector(Lx,Ly,nx,ny,precision_flag)
    n_elements = nx * ny;
    deltaX = Lx / nx;  % size of element in x-direction
    deltaY = Ly / ny;  % size of element in y-direction
    
    % single precision:
    if strcmp(precision_flag,'single')
        % X, Y and Z coordinates:
        coordCentroids = zeros(n_elements,3,'single'); % placeholder for the element coordinates
    elseif strcmp(precision_flag,'double')
        coordCentroids = zeros(n_elements,3); % placeholder for the element coordinates
    end
    
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
        % X-coord:
        coordCentroids(elem_id,1) = deltaX/2 + (ind_gridX -1)* deltaX; % increment x-coordinate
        
        % Y-coord, as we move from the row to the next row of elements
        coordCentroids(elem_id,2) = deltaY/2 + (ind_gridY -1)* deltaY; % increment y-coordinate
    
        ind_gridX = ind_gridX+1;  % increment index x as we progress
        if ind_gridX > nx  % if we exceed the number of subdivisions in x-direction
            % then reset the ind_gridX=1 and increment ind_gridY by +1
            ind_gridX = 1;
            ind_gridY = ind_gridY +1;
        end
    end
end

