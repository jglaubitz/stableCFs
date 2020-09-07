%% Test script 

%% Setting up the script 
clc, clear 

dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
points = 'Halton'; % points (equid, uniform, Halton, Legendre)

M = 8^dim; % number of points 

%% Generate and illustrate the data points 
weightFun = '1'; 
omega = generate_weightFun( weightFun, dim);
Sample = generate_points( points, domain, dim, omega, M );

figure(1)
axis equal
f1 = plot( Sample.coord(:,1), Sample.coord(:,2), 'ro', 'MarkerSize', 8); % plot points
set(f1, 'markerfacecolor', get(f1, 'color')); % use same color to fill in markers
set(gca,'FontSize',18)
