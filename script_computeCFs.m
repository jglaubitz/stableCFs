%% Test script 

%% Setting up the script 
clc, clear 

dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
weightFun = '1'; % weight function - 1, C2k, sqrt(r)
points = 'equid'; % points (equid, r, Halton) 
init_basis = 'monomials'; % monomials
%approach = 'LS'; % LS, l1

for dim=1:1

for counter3 = 1:2 
    if counter3 == 1 
        domain = 'cube';
    else 
        domain = 'ball';
    end

for counter2 = 1:2 
    if counter2 == 1 
        weightFun = '1'; 
    else 
        if strcmp( domain, 'cube')
            weightFun = 'C2k'; 
        else 
            weightFun = 'sqrt(r)';
        end
    end
            
for counter1 = 1:4 
    if counter1 == 1 
        points = 'equid';
    elseif counter1 == 2 
        points = 'uniform';
    elseif counter1 == 3 
        points = 'semi-uniform';    
    else 
        points = 'Halton';
    end

%% set up weight function
omega = generate_weightFun( weightFun, dim);
d_LS = 0; 
d_l1 = 0;

if dim == 1 
    n = 20;
    n_max = 400; 
elseif dim == 2 
    n = 2;
    n_max = 40; 
else 
    n = 1;
    n_max = 16; 
end

while n <= n_max 
    
    M = n^dim; % number of points
    
    %% generate the data points and discrete weights 
    Sample = generate_points( points, domain, dim, omega, M);
    
    if Sample.N > 0 
        %% compute the weights 
        d_start = 0;
        [ w_LS, d_LS, K_LS] = compute_cubatureWeights( Sample, domain, weightFun, 'LS', d_start); 
        %d_start = d_LS_opt;
        [ w_l1, d_l1, K_l1] = compute_cubatureWeights( Sample, domain, weightFun, 'l1', d_start); 

        [dim, counter3, counter2, counter1, n, d_LS, d_l1]
        e = -1e-14; % tollerance to allow small rounding errors
        if min(w_LS) <= e || min(w_l1) <= e 
            error('negative weight!')
        end 
        
        % save points, weights, and d and K in a matrix
        CF_LS = zeros(Sample.N,dim+2); CF_l1 = zeros(Sample.N,dim+2); % initiate
        CF_LS(:,1:dim) = Sample.coord; CF_l1(:,1:dim) = Sample.coord; % store data points 
        CF_LS(:,dim+1) = w_LS; CF_l1(:,dim+1) = w_l1; % store cubature weights 
        CF_LS(1,dim+2) = d_LS; CF_l1(1,dim+2) = d_l1; % store d
        CF_LS(2,dim+2) = K_LS; CF_l1(2,dim+2) = K_l1; % store K
        
        save( ['CFs/CF_LS_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat'], 'CF_LS' ); % safe matrix 
        save( ['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat'], 'CF_l1' ); % safe matrix
    end
       
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
end
end