function Sample = generate_points( points, domain, dim, omega, M, uniform_aux )

% This function generates the desired sample of data points. 
%   
% Input: 
%   points: String, type of data points 
%   domain: String, type of domain 
%   dim:    Integer, dimension 
%   omega:  weight function 
%   M:      number of initial points
%   uniform_aux: auxiliary vector of random points
%
% Output: 
%   Sample: Structure containing
%       N:      number of data points in the domain 
%       coord:  the (X, XY, or XYZ) coordinates of the data points 
%       

Sample.N = M; % number of points N
Sample.dim = dim;

%% Volume 
if strcmp( domain, 'cube') 
    Sample.volume = 2^dim; % volume
elseif strcmp( domain, 'ball') 
    Sample.volume = pi^(0.5*dim)/gamma(0.5*dim + 1); % volume
else
    error('Desired domain not yet implemented!')
end


%% Let us assume a cube first
%% equidistant point
if strcmp( points, 'equid')  
    n = ceil( Sample.N^(1/dim) ); % number of points in every direction 
    if n == 1 
        coord_aux = 0; 
    else
        coord_aux = linspace(-1, 1, n); 
    end
    % different dimensions 
    if dim == 1 
        Sample.coord = coord_aux';
    elseif dim == 2 
        [PointsX, PointsY] = meshgrid(coord_aux, coord_aux); 
        Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                         reshape(PointsY, numel(PointsY),1)']';
    elseif dim == 3 
        [PointsX, PointsY, PointsZ] = meshgrid(coord_aux, coord_aux, coord_aux); 
        Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                         reshape(PointsY, numel(PointsY),1)'; 
                         reshape(PointsZ, numel(PointsZ),1)' ]';
    else 
        error('Desired dimension not yet implemented!') 
    end   
    
%% (product rule) Legendre points
elseif strcmp( points, 'Legendre') 
    n = ceil( Sample.N^(1/dim) ); % number of points in every direction 
    [x,w]=lgwt(n,-1,1);
    coord_aux = flip(x'); 
    if dim == 1
        Sample.coord = coord_aux';
    elseif dim == 2 
        [PointsX, PointsY] = meshgrid(coord_aux, coord_aux); 
        Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                         reshape(PointsY, numel(PointsY),1)']';
    elseif dim == 3 
        [PointsX, PointsY, PointsZ] = meshgrid(coord_aux, coord_aux, coord_aux); 
        Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                         reshape(PointsY, numel(PointsY),1)'; 
                         reshape(PointsZ, numel(PointsZ),1)' ]';
    else 
        error('Desired dimension not yet implemented!') 
    end
        
%% uniformly distributed points 
elseif strcmp( points, 'uniform') 
    Sample.coord = uniform_aux(1:Sample.N,:);

%% semi-uniformly distributed points (add noise to equidistant points)
elseif strcmp( points, 'semi-uniform') 
    n = ceil( Sample.N^(1/dim) ); % number of points in every direction 
    if n == 1 
        coord_aux = 0; 
    else
        coord_aux = linspace(-1, 1, n); 
    end
    % different dimensions 
    if dim == 1 
        Sample.coord = coord_aux';
    elseif dim == 2 
        [PointsX, PointsY] = meshgrid(coord_aux, coord_aux); 
        Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                         reshape(PointsY, numel(PointsY),1)']';
    elseif dim == 3 
        [PointsX, PointsY, PointsZ] = meshgrid(coord_aux, coord_aux, coord_aux); 
        Sample.coord = [ reshape(PointsX, numel(PointsX),1)'; 
                         reshape(PointsY, numel(PointsY),1)'; 
                         reshape(PointsZ, numel(PointsZ),1)' ]';
    else 
        error('Desired dimension not yet implemented!') 
    end  
    % add noise 
    Sample.coord = Sample.coord + uniform_aux(1:Sample.N,:)/(4*n);
    % remove points which are outside of the ball 
    n = 1; 
    while n <= Sample.N 
        if max(abs(Sample.coord(n,:))) > 1 % point lies outside of the cube
            Sample.coord(n,:) = []; % remove point 
            Sample.N = Sample.N - 1; % decrease number of points by 1
        else
            n = n+1; % check next point
        end 
    end
    
%% Halton points 
elseif strcmp( points, 'Halton')    
    p = haltonset(dim); % generate Halton point set 
    p = scramble(p,'RR2'); % scramble point set 
    Sample.coord = 2*net(p,Sample.N)-1; % generate the first N points
        
%% otherwise 
else
    error('Desired points not yet implemented!')    
end

%% initiate discrete weights 
Sample.r = zeros( Sample.N, 1); % discrete weights
Sample.r = generate_discreteWeights( Sample, omega);

%% Cube 
if strcmp( domain, 'cube') 
    n = 1; 
    while n <= Sample.N 
        if Sample.r(n) == 0 % weight function is zero at the points
            Sample.coord(n,:) = []; % remove point 
            Sample.r(n) = [];
            Sample.N = Sample.N - 1; % decrease number of points by 1
        else
            n = n+1; % check next point
        end 
    end
    
%% Ball (only use points with radius less or equal to 1)
elseif strcmp( domain, 'ball') 
    n = 1; 
    while n <= Sample.N 
        if norm(Sample.coord(n,:)) > 1 || Sample.r(n) == 0 % point lies outside of the ball or weight function is zero there
            Sample.coord(n,:) = []; % remove point 
            Sample.r(n) = [];
            Sample.N = Sample.N - 1; % decrease number of points by 1
        else
            n = n+1; % check next point
        end 
    end
    
%% Otherwise
else
    error('Desired domain not yet implemented!')
end