clearvars; close all; clc;

nx = 50;  
ny = 50 ; 

n_elements = nx * ny;

random_vector = rand(n_elements,1); % Generate vector of uncorelated random numbers

csvwrite('phi_uncorrel.csv',random_vector)