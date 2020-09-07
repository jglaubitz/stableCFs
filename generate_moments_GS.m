%% generate_moments_GS
% Generates the moments corresponding to the basis of DOPs \pi_k 
% We follow the Gram-Schmidt procedure to compute the moments of the DOPs
% 
% INPUT: 
%  Sample :     Sample for the data points
%  d :          degree of exactness 
%  domain :     domain 
%  weightFun :  weight function 
%
% OUTPUT: 
%  A : matrix which contains the values of \pi_k at the points Y 

function m = generate_moments_GS( Sample, d, domain, weightFun )

    N = Sample.N; % number of data points 
    X = Sample.coord; % data points
    K = nchoosek(Sample.dim + d, Sample.dim); % binomial coefficient/ number of DOPs
    B = zeros(K,N); % auxiliary matrix 
    m = zeros(K,1); % moments 
    
    %% moments of the monomials 
    m_monom = zeros(K,1); 
    % cube 
    if strcmp( domain, 'cube') 
        % omega = 1
        if strcmp( weightFun, '1') 
            if Sample.dim == 1
                for k=0:2:d 
                   m_monom(k+1) = 2/(k+1); % moment 
                end
            elseif Sample.dim == 2 
                k = 1; 
                for k1=0:1:d
                for k2=0:1:d 
                    if k1+k2<=d 
                        % k1 and k2 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 
                            m_monom(k) = 4/( (k1+1)*(k2+1) ); % moment 
                        end
                        k = k+1;
                    end
                end
                end
            elseif Sample.dim == 3 
                k = 1; 
                for k1=0:d
                for k2=0:d 
                for k3=0:d
                    if k1+k2+k3<=d 
                        % k1, k2, and k3 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 && mod(k3,2) == 0
                            m_monom(k) = 8/( (k1+1)*(k2+1)*(k3+1) ); % moment 
                        end
                        k = k+1;
                    end
                end
                end
                end
            else 
                error('Desired dimension not yet implemented!') 
            end
            
        % omega = sqrt(1-x^2) (C2k)
        elseif strcmp( weightFun, 'C2k') 
            m_aux = zeros(d+1,1); 
            m_aux(1) = 0.5*pi; 
            for k=2:2:d
                m_aux(k+1) = ( (k-1)/(k+2) )*m_aux(k-1); 
            end 
            if Sample.dim == 1 
                m_monom = m_aux; % moments
            elseif Sample.dim == 2 
                k = 1; 
                for k1=0:1:d
                for k2=0:1:d 
                    if k1+k2<=d 
                        m_monom(k) = m_aux(k1+1)*m_aux(k2+1); % moment 
                        k = k+1;
                    end
                end
                end
            elseif Sample.dim == 3 
                k = 1; 
                for k1=0:d
                for k2=0:d 
                for k3=0:d
                    if k1+k2+k3<=d 
                        m_monom(k) = m_aux(k1+1)*m_aux(k2+1)*m_aux(k3+1); % moment 
                        k = k+1;
                    end
                end
                end
                end
            else 
                error('Desired dimension not yet implemented!') 
            end
            
        else 
            error('Desired weight function not yet implemented!')
        end
        
    % ball
    elseif strcmp( domain, 'ball') 
        % omega = 1
        if strcmp( weightFun, '1') 
            if Sample.dim == 1
                for k=0:2:d 
                   m_monom(k+1) = 2/(k+1); % moment 
                end
            elseif Sample.dim == 2 
                k = 1; 
                for k1=0:d
                for k2=0:d
                    if k1+k2<=d 
                        % k1 and k2 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 
                            b1 = 0.5*(k1+1); 
                            b2 = 0.5*(k2+1);
                            m_monom(k) = ( 2*gamma(b1)*gamma(b2)/gamma(b1+b2) )/(k1+k2+Sample.dim); % moment 
                        end
                        k = k+1;
                    end
                end
                end
            elseif Sample.dim == 3 
                k = 1; 
                for k1=0:d
                for k2=0:d 
                for k3=0:d
                    if k1+k2+k3<=d 
                        % k1, k2, and k3 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 && mod(k3,2) == 0
                            b1 = 0.5*(k1+1); 
                            b2 = 0.5*(k2+1);
                            b3 = 0.5*(k3+1);
                            m_monom(k) = ( 2*gamma(b1)*gamma(b2)*gamma(b3)/gamma(b1+b2+b3) )/(k1+k2+k3+Sample.dim); % moment 
                        end
                        k = k+1;
                    end
                end
                end
                end
            else 
                error('Desired dimension not yet implemented!') 
            end
        
        % omega = sqrt(r)
        elseif strcmp( weightFun, 'sqrt(r)') 
            if Sample.dim == 1
                for k=0:2:d 
                   m_monom(k+1) = 2/(k+Sample.dim+0.5); % moment 
                end
            elseif Sample.dim == 2 
                k = 1; 
                for k1=0:d
                for k2=0:d
                    if k1+k2<=d 
                        % k1 and k2 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 
                            b1 = 0.5*(k1+1); 
                            b2 = 0.5*(k2+1);
                            m_monom(k) = ( 2*gamma(b1)*gamma(b2)/gamma(b1+b2) )/(k1+k2+Sample.dim+0.5); % moment 
                        end
                        k = k+1;
                    end
                end
                end
            elseif Sample.dim == 3 
                k = 1; 
                for k1=0:d
                for k2=0:d 
                for k3=0:d
                    if k1+k2+k3<=d 
                        % k1, k2, and k3 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 && mod(k3,2) == 0
                            b1 = 0.5*(k1+1); 
                            b2 = 0.5*(k2+1);
                            b3 = 0.5*(k3+1);
                            m_monom(k) = ( 2*gamma(b1)*gamma(b2)*gamma(b3)/gamma(b1+b2+b3) )/(k1+k2+k3+Sample.dim+0.5); % moment 
                        end
                        k = k+1;
                    end
                end
                end
                end
            else 
                error('Desired dimension not yet implemented!') 
            end
            
        else
            error('Desired weight function not yet implemented!')
        end
        
    % else 
    else
        error('Desired domain not yet implemented!')
    end
    
    %% 1st DOP and moment 
    x_aux = ones(1,N); 
    x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
    B(1,:) = x_aux./x_aux_norm; % normalize
    m(1) = m_monom(1)/x_aux_norm; % moment
    
    %% dim = 1 
    if Sample.dim == 1
        for k=2:d+1 
            x_aux = X.^(k-1)'; % values of e_k at the data points 
            m(k) = m_monom(k); % moment
            % GS procedure
            for i=1:k 
                inner_prod = dot(x_aux.*B(i,:),Sample.r); % inner product
                x_aux = x_aux - inner_prod*B(i,:); % modified GM 
                m(k) = m(k) - inner_prod*m(i); % moment 
            end
            x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
            B(k,:) = x_aux./x_aux_norm; % normalize
            m(k) = m(k)/x_aux_norm; % moment 
        end

    %% dim = 2 
    elseif Sample.dim == 2
        k = 2; 
        for k1=0:d
        for k2=0:d
            if k1+k2>=1 && k1+k2<=d 
                x_aux = ( (X(:,1).^k1).*(X(:,2).^k2) )'; % values of e_k at the data points
                m(k) = m_monom(k); % moment
                % GS procedure
                for i = 1:k % modified GS
                    inner_prod = dot(x_aux.*B(i,:),Sample.r); % inner product
                    x_aux = x_aux - inner_prod*B(i,:); % modified GM
                    m(k) = m(k) - inner_prod*m(i); % moment 
                end
                x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
                B(k,:) = x_aux./x_aux_norm; % normalize
                m(k) = m(k)/x_aux_norm; % moment 
                k = k+1;
            end
        end
        end

    %% dim = 3 
    elseif Sample.dim == 3
        k = 2; 
        for k1=0:d
        for k2=0:d 
        for k3=0:d
            if k1+k2+k3>=1 && k1+k2+k3<=d 
                x_aux = ( (X(:,1).^k1).*(X(:,2).^k2).*(X(:,3).^k3) )'; % values of e_k at the data points
                m(k) = m_monom(k); % moment
                % GS procedure
                for i = 1:k % modified GS
                    inner_prod = dot(x_aux.*B(i,:),Sample.r); % inner product
                    x_aux = x_aux - inner_prod*B(i,:); % modified GM
                    m(k) = m(k) - inner_prod*m(i); % moment 
                end
                x_aux_norm = sqrt( dot(x_aux.^2,Sample.r) ); % norm
                B(k,:) = x_aux./x_aux_norm; % normalize
                m(k) = m(k)/x_aux_norm; % moment 
                k = k+1;
            end
        end
        end
        end

    else 
        error('Desired dimension not yet implemented!') 
    end
    
end