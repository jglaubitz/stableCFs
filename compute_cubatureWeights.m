%% compute_cubatureWeights
% Generates the vector of cubature weights. 
% 
% INPUT: 
%  Sample :     sample of data points
%  domain :     (integration) domain 
%  approach :   approach to compute the weights (LS, l1) 
%  d_start :    start value for d
%
% OUTPUT: 
%  w         : vector of cubature weights 
%  d       : highest possible degree of exactness 
%  K       : corresponding number of basis functions

function [w,d,K] = compute_cubatureWeights( Sample, domain, weightFun, approach, d_start)

    R = diag(Sample.r); % discrete weight matrix 
    w = zeros(Sample.N,1); % cubature weights 
    
    %% set up inital values 
    d = d_start; % degree of exactness (DoE)
    K = nchoosek(Sample.dim + d, Sample.dim); % number of basis elements 
    check_rank = K; % rank 
    w_min = min(w); % smallest cubature weight 
    
    %% loop in which d is increased until the CF becomes unstable 
    e = -1e-14; % tollerance to allow small rounding errors
    while check_rank == K && w_min >= e 
          
        d = d+1; % increase the DoE
        K = nchoosek(Sample.dim + d, Sample.dim); % number of basis elements 
        [ basis, m ] = generate_monomials( Sample.dim, domain, weightFun, d ); % monomials and their moments 
        P = basis(Sample.coord); % Vandermonde matrix
        check_rank = rank(P); % rank of P
        
        % check wheter the data points are d-unisolvent
        if check_rank == K 
        
            % different approaches 
            if strcmp( approach, 'LS') 
                A = P*sqrt(R); 
                v = lsqminnorm(A,m); % indirect computation using optimization tools 
                w = sqrt(R)*v;
            elseif strcmp( approach, 'l1') 
                options = optimoptions('linprog','Display','none');
                w = linprog( ones(Sample.N,1), [], [], P, m, zeros(Sample.N,1), [], options ); % l1 weights
            else 
                error('Unkwon approach!') 
            end
            
        end 
        
        if sum(w)==sum(w) && sum(abs(w))~=Inf && length(w) > 0 % w does not contain NaNs or Infs and is nonempty
            w_min = min(w); 
        else 
           w_min = -1; 
        end
        
    end
    
    % once the CF becomes unstable (is not unisolvent anymore) we go back
    % to the previous d 
    d = d-1; % increase the DoE
    K = nchoosek(Sample.dim + d, Sample.dim); % number of basis elements 
    [ basis, m ] = generate_monomials( Sample.dim, domain, weightFun, d ); % monomials and their moments 
   	P = basis(Sample.coord); % Vandermonde matrix
    % different approaches 
   	if strcmp( approach, 'LS') 
     	A = P*sqrt(R); 
     	v = lsqminnorm(A,m); % indirect computation using optimization tools 
       	w = sqrt(R)*v;
  	elseif strcmp( approach, 'l1') 
       	options = optimoptions('linprog','Display','none');
      	w = linprog( ones(Sample.N,1), [], [], P, m, zeros(Sample.N,1), [], options ); % l1 weights
    else 
      	error('Unkwon approach!') 
   	end
            
end