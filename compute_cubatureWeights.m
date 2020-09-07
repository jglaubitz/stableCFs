%% compute_cubatureWeights
% Generates the vector of cubature weights. 
% We use the modified Gram-Schmidt process! 
% As an intial basis, we use the Legendre polynomials L_k. 
% 
% INPUT: 
%  Sample :     sample of data points
%  domain :     (integration) domain 
%  init_basis : initial basis for the GS procedure 
%  approach :   approach to compute the weights (LS_direct, LS_opt, l1) 
%  d_start :    start value for d
%
% OUTPUT: 
%  w         : vector of cubature weights

function [w,d,K] = compute_cubatureWeights( Sample, domain, weightFun, init_basis, approach, d_start)

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
        
        check_rank = -1;
        if K <= Sample.N 
            P = dopBasis( Sample, d, Sample.coord, init_basis); % Vandermonde matrix 
            if P==P % P does not contain NaNs
                check_rank = rank(P); % rank of P 
            end       
        end
        
        % check wheter the data points are d-unisolvent
        if check_rank == K 
        
            m = generate_moments_GS( Sample, d, domain, weightFun ); % moments
            % different approaches 
            if strcmp( approach, 'LS_direct') 
                w = R*(P')*m; % direct computation using DOPs 
            elseif strcmp( approach, 'LS_opt') 
                A = P*sqrt(R); 
                v = lsqminnorm(A,m); % indirect computation using optimization tools 
                w = sqrt(R)*v; 
            elseif strcmp( approach, 'l1') 
                w = minL1lin( eye(Sample.N),zeros(Sample.N,1), [],[], P,m); % l1 weights 
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
    P = dopBasis( Sample, d, Sample.coord, init_basis); % Vandermonde matrix 
    m = generate_moments_GS( Sample, d, domain, weightFun ); % moments
    % different approaches 
    if strcmp( approach, 'LS_direct') 
        w = R*(P')*m; % direct computation using DOPs 
    elseif strcmp( approach, 'LS_opt') 
        A = P*sqrt(R); 
        v = lsqminnorm(A,m); % indirect computation using optimization tools 
        w = sqrt(R)*v; 
    elseif strcmp( approach, 'l1') 
        w = minL1lin( eye(Sample.N),zeros(Sample.N,1), [],[], P,m); % l1 weights
    else 
        error('Unkwon approach!') 
    end
    
end