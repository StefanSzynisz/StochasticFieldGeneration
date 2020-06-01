%This script shows how to correlate two arrays using the following
%methods: svd decomposition, eigendecomposition and cholesky decomposition
clc; close all; clear all

%%
%To be defined by the user:

%select correlation: from -1 to 1
corr = 0.8;

%select decomposition: svd = 0, eigen = 1, cholesky = 2,
decomposition = 1;

%%
% covariance matrix
cov_y1_y2 = [1 corr; corr 1];

%Decomposition
if (decomposition==0)
    [U,S,~] = svd(cov_y1_y2); %svd
    B = U*sqrt(S);
    disp('You are using SVD decomposition')
    
elseif (decomposition==1)
    [V,D] = eig(cov_y1_y2); %eigen
    B = V * sqrt( abs(D) );
    disp('You are using Eigen decomposition')
    
elseif (decomposition==2)
    L = chol(cov_y1_y2); %cholesky
    B = L';
    disp('You are using Cholesky decomposition')
    
end

%create two normal and uncorrelated random arrays
y1 = randn(10000,1);
y2 = randn(10000,1);

%check their correlation
corrplot([y1 y2])

%append them
y = [y1,y2];

%correlate arrays based on selected decomposition
correlate = (B*y')';

%check correlation between arrays
corrplot(correlate)

%check covariance matrix between arrays
cov(correlate)
disp('This covariance matrix should have the correlation you selected!')



