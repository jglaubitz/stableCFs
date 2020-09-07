%% dopBasis 
% Generates a basis of DOPs \pi_k and evaluates them at some points. 
% We use the modified Gram-Schmidt process! 
% As an intial basis, we use the Legendre polynomials L_k. 
% 
% INPUT: 
%  Sample : Sample for the data points
%  d :      maximal degree 
%  Y :      array of evaluation points 
%  init_basis : initial basis from which the DOPs are constructed
%
% OUTPUT: 
%  A : matrix which contains the values of \pi_k at the points Y 

function A = dopBasis( Sample, d, Y, init_basis)

    N = Sample.N; % number of data points 
    X = Sample.coord; % data points
    [E, dim] = size(Y); % number of evaluation points 
    K = nchoosek(dim + d, dim); % binomial coefficient/ number of DOPs
    B = zeros(K,N); % auxiliary matrix 
    A = zeros(K,E); % initiate matrix A

    %% 1st DOP
    x_aux = ones(1,N); 
    y_aux = ones(1,E);
    x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
    B(1,:) = x_aux./x_aux_norm; % normalize
    A(1,:) = y_aux./x_aux_norm; % normalize

    %% dim = 1 
    if dim == 1
        for k=2:d+1 
            % select init basis
            if strcmp( init_basis, 'Legendre') 
                x_aux = legendreP(k-1,X)'; % values of L_k at the data points
                y_aux = legendreP(k-1,Y)'; % values of L_k at the evaluation points 
            elseif strcmp( init_basis, 'monomials') 
                x_aux = X.^(k-1)'; % values of e_k at the data points
                y_aux = Y.^(k-1)'; % values of e_k at the evaluation points
            else 
                error('Desired inital basis not yet implemented!')
            end
            % GS procedure
            for i=1:k 
                inner_prod = dot(x_aux.*B(i,:),Sample.r); % inner product
                y_aux = y_aux - inner_prod*A(i,:); % modified GM
                x_aux = x_aux - inner_prod*B(i,:); % modified GM
            end
            x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
            B(k,:) = x_aux./x_aux_norm; % normalize
            A(k,:) = y_aux./x_aux_norm; % normalize
        end

    %% dim = 2 
    elseif dim == 2
        k = 2; 
        for k1=0:d
        for k2=0:d
            if k1+k2>=1 && k1+k2<=d 
                % select init basis
                if strcmp( init_basis, 'Legendre') 
                    x_aux = ( legendreP(k1,X(:,1)).*legendreP(k2,X(:,2)) )'; % values of L_k at the data points
                    y_aux = ( legendreP(k1,Y(:,1)).*legendreP(k2,Y(:,2)) )'; % values of L_k at the evaluation points
                elseif strcmp( init_basis, 'monomials') 
                    x_aux = ( (X(:,1).^k1).*(X(:,2).^k2) )'; % values of e_k at the data points
                    y_aux = ( (Y(:,1).^k1).*(Y(:,2).^k2) )'; % values of e_k at the evaluation points
                else 
                    error('Desired inital basis not yet implemented!')
                end
                % GS procedure
                for i = 1:k % modified GS
                    inner_prod = dot(x_aux.*B(i,:),Sample.r); % inner product
                    y_aux = y_aux - inner_prod*A(i,:); % modified GM
                    x_aux = x_aux - inner_prod*B(i,:); % modified GM
                end
                x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
                B(k,:) = x_aux./x_aux_norm; % normalize
                A(k,:) = y_aux./x_aux_norm; % normalize 
                k = k+1;
            end
        end
        end

    %% dim = 3 
    elseif dim == 3
        k = 2; 
        for k1=0:d
        for k2=0:d 
        for k3=0:d
            if k1+k2+k3>=1 && k1+k2+k3<=d 
                % select init basis
                if strcmp( init_basis, 'Legendre') 
                    x_aux = ( legendreP(k1,X(:,1)).*legendreP(k2,X(:,2)).*legendreP(k3,X(:,3)) )'; % values of L_k at the data points
                    y_aux = ( legendreP(k1,Y(:,1)).*legendreP(k2,Y(:,2)).*legendreP(k3,Y(:,3)) )'; % values of L_k at the evaluation points
                elseif strcmp( init_basis, 'monomials') 
                    x_aux = ( (X(:,1).^k1).*(X(:,2).^k2).*(X(:,3).^k3) )'; % values of e_k at the data points
                    y_aux = ( (Y(:,1).^k1).*(Y(:,2).^k2).*(Y(:,3).^k3) )'; % values of e_k at the evaluation points
                else 
                    error('Desired inital basis not yet implemented!')
                end
                % GS procedure
                for i = 1:k % modified GS
                    inner_prod = dot(x_aux.*B(i,:),Sample.r); % inner product
                    y_aux = y_aux - inner_prod*A(i,:); % modified GM
                    x_aux = x_aux - inner_prod*B(i,:); % modified GM
                end
                x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
                B(k,:) = x_aux./x_aux_norm; % normalize
                A(k,:) = y_aux./x_aux_norm; % normalize 
                k = k+1;
            end
        end
        end
        end

    else 
        error('Desired dimension not yet implemented!') 
    end

end