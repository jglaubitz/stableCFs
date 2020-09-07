%% Test script 

%% Setting up the script 
clc, clear 

dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
weightFun = '1'; % weight function - 1, C2k, sqrt(r)
points = 'semi-uniform'; % points (equid, uniform, semi-uniform, Halton) 
init_basis = 'monomials'; % monomials
%approach = 'LS_direct'; % LS_direct, LS_opt, l1

for dim=1:3

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
d_LS_dir = 0; 
d_LS_opt = 0;
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

uniform_aux = 2*rand(n_max^dim,dim)-1; % generate uniformly distrubuted points
d_start = 0;

while n <= n_max 
    
    M = n^dim; % number of points
    
    %% generate the data points and discrete weights 
    Sample = generate_points( points, domain, dim, omega, M, uniform_aux);
    
    if Sample.N > 0 
        %% compute the weights 
        d_start = 0;
        [ w_LS_dir, d_LS_dir, K_LS_dir] = compute_cubatureWeights( Sample, domain, weightFun, init_basis, 'LS_direct', d_start); 
        [ w_LS_opt, d_LS_opt, K_LS_opt] = compute_cubatureWeights( Sample, domain, weightFun, init_basis, 'LS_opt', d_start);
        %d_start = d_LS_opt;
        [ w_l1, d_l1, K_l1] = compute_cubatureWeights( Sample, domain, weightFun, init_basis, 'l1', d_start); 

        [dim, counter3, counter2, counter1, n, d_LS_dir, d_LS_opt, d_l1, min(w_LS_dir), min(w_LS_opt), min(w_l1)]
        e = -1e-14; % tollerance to allow small rounding errors
        if min(w_LS_dir) <= e || min(w_LS_opt) <= e || min(w_l1) <= e 
            error('negative weight!')
        end 
        
        % save points, weights, and d and K in a matrix
        CF_LS_dir = zeros(Sample.N,dim+2); % initiate 
        CF_LS_opt = zeros(Sample.N,dim+2); % initiate
        CF_l1 = zeros(Sample.N,dim+2); % initiate
        CF_LS_dir(:,1:dim) = Sample.coord; % store data points
        CF_LS_opt(:,1:dim) = Sample.coord; % store data points
        CF_l1(:,1:dim) = Sample.coord; % store data points 
        CF_LS_dir(:,dim+1) = w_LS_dir; % store cubature weights 
        CF_LS_opt(:,dim+1) = w_LS_opt; % store cubature weights
        CF_l1(:,dim+1) = w_l1; % store cubature weights 
        CF_LS_dir(1,dim+2) = d_LS_dir; % store d
        CF_LS_opt(1,dim+2) = d_LS_opt; % store d
        CF_l1(1,dim+2) = d_l1; % store d
        CF_LS_dir(2,dim+2) = d_LS_dir; % store K
        CF_LS_opt(2,dim+2) = d_LS_opt; % store K
        CF_l1(2,dim+2) = d_l1; % store K
        
        save( ['CFs/CF_LS_dir_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat'], 'CF_LS_dir' ); % safe matrix 
        save( ['CFs/CF_LS_opt_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat'], 'CF_LS_opt' ); % safe matrix
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