%% compute_LegendreRule
% Generates the vectors of Legender points and weights. 
% 
% INPUT: 
%  dim:    dimension
%  domain: integration domain 
%  n:      number of points in every direction 
%
% OUTPUT: 
%  w: vector of cubature weights 
%  d: degree of exavtness 
%  K: number of basis functions 

function [ X, W, d, K ] = compute_LegendreRule( dim, domain, n )

    %% cube 
    if strcmp( domain, 'cube') 
        M = n; % number of evaluation points in every direction 
        Y = zeros(M^dim,dim); % evaluation points 
        W = zeros(M^dim,1); % auxilary cubature weights
        [x,w]=lgwt(M,-1,1); 
        x = flip(x'); 
        w = flip(w');
        if dim == 1
            Y = x'; % evaluation points 
            W = w'; % weights
        elseif dim == 2 
            % evaluation points
            [PointsX, PointsY] = meshgrid(x, x); 
            Y = [ reshape(PointsX, numel(PointsX),1)'; 
                  reshape(PointsY, numel(PointsY),1)']'; 
            % weights 
            [weightsX, weightsY] = meshgrid(w, w); 
            W_aux = [ reshape(weightsX, numel(weightsX),1)'; 
                      reshape(weightsY, numel(weightsY),1)']'; 
            W = W_aux(:,1).*W_aux(:,2);        
        elseif dim == 3 
            % evaluation points
            [PointsX, PointsY, PointsZ] = meshgrid(x, x, x); 
            Y = [ reshape(PointsX, numel(PointsX),1)'; 
                  reshape(PointsY, numel(PointsY),1)';
                  reshape(PointsZ, numel(PointsZ),1)']'; 
            % weights 
            [weightsX, weightsY, weightsZ] = meshgrid(w, w, w); 
            W_aux = [ reshape(weightsX, numel(weightsX),1)'; 
                      reshape(weightsY, numel(weightsY),1)';
                      reshape(weightsZ, numel(weightsZ),1)']'; 
            W = W_aux(:,1).*W_aux(:,2).*W_aux(:,3);
        else 
        	error('Desired dimension not yet implemented!') 
        end
        
    %% ball 
    elseif strcmp( domain, 'ball')  
        M = n; % number of evaluation points in every direction 
        Y = zeros(M^dim,dim); % evaluation points 
        W = zeros(M^dim,1); % auxilary cubature weights
        [x,w]=lgwt(M,-1,1); 
        x = flip(x'); 
        w = flip(w');
    	if dim == 1 
            Y = x'; % evaluation points 
            W = w'; % weights
        elseif dim == 2 
            [r,w_r] = lgwt(M,0,1); % radius
            [phi,w_phi] = lgwt(M,0,2*pi); % angle
            r = flip(r'); 
            w_r = flip(w_r'); 
            phi = flip(phi'); 
            w_phi = flip(w_phi');
            % evaluation points - product rule
            [PointsR, PointsPhi] = meshgrid(r, phi); 
            RPhi = [ reshape(PointsR, numel(PointsR),1)'; 
            reshape(PointsPhi, numel(PointsPhi),1)']'; 
            % weights - product rule 
            [weightsR, weightsPhi] = meshgrid(w_r, w_phi); 
            W_RPhi_aux = [ reshape(weightsR, numel(weightsR),1)'; 
                           reshape(weightsPhi, numel(weightsPhi),1)']'; 
            W_RPhi = W_RPhi_aux(:,1).*W_RPhi_aux(:,2); 
            % transform to the ball 
            Y(:,1) = RPhi(:,1).*cos( RPhi(:,2) ); % x coordinate 
            Y(:,2) = RPhi(:,1).*sin( RPhi(:,2) ); % y coordinate 
            W = W_RPhi.*RPhi(:,1); % weights
        elseif dim == 3 
            [r,w_r] = lgwt(M,0,1); % radius
            [phi1,w_phi1] = lgwt(M,0,pi); % angle 1
            [phi2,w_phi2] = lgwt(M,0,2*pi); % angle 2
            r = flip(r'); 
            w_r = flip(w_r'); 
            phi1 = flip(phi1'); 
            w_phi1 = flip(w_phi1');
            phi2 = flip(phi2'); 
            w_phi2 = flip(w_phi2');
            % evaluation points - product rule
            [PointsR, PointsPhi1, PointsPhi2] = meshgrid(r, phi1, phi2); 
            RPhi = [ reshape(PointsR, numel(PointsR),1)'; 
                     reshape(PointsPhi1, numel(PointsPhi1),1)'; 
                     reshape(PointsPhi2, numel(PointsPhi2),1)']'; 
            % weights - product rule 
            [weightsR, weightsPhi1, weightsPhi2] = meshgrid(w_r, w_phi1, w_phi2); 
            W_RPhi_aux = [ reshape(weightsR, numel(weightsR),1)'; 
                           reshape(weightsPhi1, numel(weightsPhi1),1)'; 
                           reshape(weightsPhi2, numel(weightsPhi2),1)']'; 
            W_RPhi = W_RPhi_aux(:,1).*W_RPhi_aux(:,2).*W_RPhi_aux(:,3); 
            % transform to the ball 
            Y(:,1) = RPhi(:,1).*sin( RPhi(:,2) ).*cos( RPhi(:,3) ); % x coordinate 
            Y(:,1) = RPhi(:,1).*sin( RPhi(:,2) ).*sin( RPhi(:,3) ); % y coordinate 
            Y(:,1) = RPhi(:,1).*cos( RPhi(:,2) ); % z coordinate
            W = W_RPhi.*(RPhi(:,1).^2).*sin( RPhi(:,2) ); % weights
        else 
            error('Desired dimension not yet implemented!') 
        end
        
    else 
        error('Desired domain not yet implemented!')
end
    
    X = Y;
    d = 2*n-1; % degree if exactness 
    K = (d+1)^dim; % number of basis functions 
    
end