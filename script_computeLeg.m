%% Script to compute (transformed) product Legendre rules 

%% Setting up the script 
clc, clear 

dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 

for dim=1:3

for counter3 = 1:2 
    if counter3 == 1 
        domain = 'cube';
    else 
        domain = 'ball';
    end

%% set up weight function
omega = generate_weightFun( '1', dim);

if dim == 1 
    n = 20;
    n_max = 400; 
    n_lenght = 20;
elseif dim == 2 
    n = 2;
    n_max = 40; 
    n_lenght = 20;
else 
    n = 1;
    n_max = 16; 
    n_lenght = 16;
end

while n <= n_max 
        
    [X, w_Leg, d_Leg, K_Leg ] = compute_LegendreRule( dim, domain, n );
    [dim, counter3, n ]
    
    % save points, weights, and d and K in a matrix
    CF_Leg = zeros(n^dim,dim+2); % initiate 
    CF_Leg(:,1:dim) = X; % store data points
    CF_Leg(:,dim+1) = w_Leg; % store cubature weights 
    CF_Leg(1,dim+2) = d_Leg; % store d
    CF_Leg(2,dim+2) = d_Leg; % store K      
    save( ['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_n=',num2str(n),'.mat'], 'CF_Leg' ); % safe matrix   
        
    % increase n
    if dim == 1 
        n = n + 20;
    elseif dim == 2 
        n = n + 2;
    else 
        n = n + 1;
    end
    
end

end
end