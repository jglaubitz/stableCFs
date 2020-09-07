%% generate_discreteWeights
% Generates the vector of cubature weights. 
% We use the modified Gram-Schmidt process! 
% As an intial basis, we use the Legendre polynomials L_k. 
% 
% INPUT: 
%  Sample    : Sample of data points
%  weightFun : weight function 
%
% OUTPUT: 
%  w         : vector of cubature weights

function r = generate_discreteWeights( Sample, omega )

    % evaluate the weight function at the data points 
    omega_evaluated = zeros(Sample.N,1);
    for n = 1:Sample.N 
       if Sample.dim == 1  
            omega_evaluated(n) = omega( Sample.coord(n,1) ); 
       elseif Sample.dim == 2  
            omega_evaluated(n) = omega( Sample.coord(n,1), Sample.coord(n,2) );
       elseif Sample.dim == 3  
            omega_evaluated(n) = omega( Sample.coord(n,1), Sample.coord(n,2) , Sample.coord(n,3) );
       else 
            error('Desired dimension not yet implemented!') 
       end
    end

    % generate the discrete weights 
    r = omega_evaluated*Sample.volume/Sample.N;
    
end