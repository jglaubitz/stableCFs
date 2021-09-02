%% generate_monomials
% Generates the basis of monomials and the corresponding moments
% 
% INPUT: 
%  dim :        dimension 
%  domain :     domain 
%  weightFun :  weight function 
%  d :          degree of exactness 
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 

function [ basis, m ] = generate_monomials( dim, domain, weightFun, d )

    K = nchoosek(dim + d, dim); % binomial coefficient/ dimension
    alpha1 = zeros(K,1); % vector of exponents for x
    alpha2 = zeros(K,1); % vector of exponents for y
    alpha3 = zeros(K,1); % vector of exponents for z 
    m = zeros(K,1); % moments 
    
    %% exponents and basis 
    if dim == 1 
        alpha1 = (0:d)'; 
        basis = @(x) x'.^alpha1; % basis 
    elseif dim == 2 
        k = 1; 
        for k1=0:1:d
        for k2=0:1:d 
            if k1+k2<=d 
                alpha1(k) = k1; 
                alpha2(k) = k2; 
                k = k+1;
            end
        end
        end 
        basis = @(x) x(:,1)'.^alpha1 .* x(:,2)'.^alpha2; % basis 
    elseif dim == 3 
    	k = 1; 
        for k1=0:d
        for k2=0:d 
        for k3=0:d
        	if k1+k2+k3<=d 
            	alpha1(k) = k1; 
                alpha2(k) = k2;
                alpha3(k) = k3; 
                k = k+1;
            end
        end
        end
        end 
        basis = @(x) x(:,1)'.^alpha1 .* x(:,2)'.^alpha2 .* x(:,3)'.^alpha3; % basis
    else 
    	error('Desired dimension not yet implemented!') 
    end   
    
    %% moments of the monomials 
    % cube 
    if strcmp( domain, 'cube') 
        % omega = 1
        if strcmp( weightFun, '1') 
            if dim == 1 
                for k=0:2:d 
                   m(k+1) = 2/(k+1); % moment 
                end
            elseif dim == 2 
                k = 1; 
                for k1=0:1:d
                for k2=0:1:d 
                    if k1+k2<=d 
                        % k1 and k2 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 
                            m(k) = 4/( (k1+1)*(k2+1) ); % moment 
                        end
                        k = k+1;
                    end
                end
                end
            elseif dim == 3 
                k = 1; 
                for k1=0:d
                for k2=0:d 
                for k3=0:d
                    if k1+k2+k3<=d 
                        % k1, k2, and k3 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 && mod(k3,2) == 0
                            m(k) = 8/( (k1+1)*(k2+1)*(k3+1) ); % moment 
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
            if dim == 1 
                m = m_aux; % moments
            elseif dim == 2 
                k = 1; 
                for k1=0:1:d
                for k2=0:1:d 
                    if k1+k2<=d 
                        m(k) = m_aux(k1+1)*m_aux(k2+1); % moment 
                        k = k+1;
                    end
                end
                end
            elseif dim == 3 
                k = 1; 
                for k1=0:d
                for k2=0:d 
                for k3=0:d 
                    if k1+k2+k3<=d 
                        m(k) = m_aux(k1+1)*m_aux(k2+1)*m_aux(k3+1); % moment 
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
            if dim == 1 
                for k=0:2:d 
                   m(k+1) = 2/(k+1); % moment 
                end
            elseif dim == 2 
                k = 1; 
                for k1=0:d
                for k2=0:d
                    if k1+k2<=d 
                        % k1 and k2 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 
                            b1 = 0.5*(k1+1); 
                            b2 = 0.5*(k2+1);
                            m(k) = ( 2*gamma(b1)*gamma(b2)/gamma(b1+b2) )/(k1+k2+dim); % moment 
                        end
                        k = k+1;
                    end
                end
                end
            elseif dim == 3 
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
                            m(k) = ( 2*gamma(b1)*gamma(b2)*gamma(b3)/gamma(b1+b2+b3) )/(k1+k2+k3+dim); % moment 
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
            if dim == 1 
                for k=0:2:d 
                   m(k+1) = 2/(k+dim+0.5); % moment 
                end
            elseif dim == 2 
                k = 1; 
                for k1=0:d
                for k2=0:d
                    if k1+k2<=d 
                        % k1 and k2 need to be even
                        if mod(k1,2) == 0 && mod(k2,2) == 0 
                            b1 = 0.5*(k1+1); 
                            b2 = 0.5*(k2+1);
                            m(k) = ( 2*gamma(b1)*gamma(b2)/gamma(b1+b2) )/(k1+k2+dim+0.5); % moment 
                        end
                        k = k+1;
                    end
                end
                end
            elseif dim == 3 
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
                            m(k) = ( 2*gamma(b1)*gamma(b2)*gamma(b3)/gamma(b1+b2+b3) )/(k1+k2+k3+dim+0.5); % moment 
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
    
    % combi
    elseif strcmp( domain, 'combi') 
        if dim ~= 2 
            error('Desired dimension not yet implemented!')
        end 
        % omega = 1
        if strcmp( weightFun, '1') 
            % m_ball and m_cube
            m_ball = zeros(K,1);
            m_cube = zeros(K,1);
            k = 1; 
            for k1=0:d
            for k2=0:d
            	if k1+k2<=d 
                    % m_ball 
                	% k1 and k2 need to be even
                    if mod(k1,2) == 0 && mod(k2,2) == 0 
                    	b1 = 0.5*(k1+1); 
                        b2 = 0.5*(k2+1);
                        m_ball(k) = ( 2*gamma(b1)*gamma(b2)/gamma(b1+b2) )/(k1+k2+dim); % moment 
                    end
                    % m_cube 
                    m_cube(k) = ( (2^(k1+1)-1)/(k1+1) )*( (2^(k2+1)-1)/(k2+1) );
                    k = k+1;
                end
            end
            end              
            % sum up 
            m = m_ball + m_cube; 
        else
            error('Desired weight function not yet implemented!')
        end
    
    % else 
    else
        error('Desired domain not yet implemented!')
    end
    
end