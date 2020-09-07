%% generate_weights
% Generates the weight function .  
% 
% INPUT: 
%  weightFun : weight function 
%  dim :       weight function  
%
% OUTPUT: 
%  omega : weight function

function omega = generate_weightFun( weightFun, dim)

    %% 1 
    if strcmp( weightFun, '1') 
        if dim == 1
            omega = @(x) 1 + 0*x; 
        elseif dim == 2 
            omega = @(x,y) 1 + 0*x.*y; 
        elseif dim == 3 
            omega = @(x,y,z) 1 + 0*x.*y.*z; 
        else 
            error('Desired dimension not yet implemented!')
        end 
    
    %% C2k
    elseif strcmp( weightFun, 'C2k') 
        if dim == 1
            omega = @(x) sqrt(1 - x.^2); 
        elseif dim == 2 
            omega = @(x,y) sqrt(1 - x.^2).*sqrt(1 - y.^2); 
        elseif dim == 3 
            omega = @(x,y,z) sqrt(1 - x.^2).*sqrt(1 - y.^2).*sqrt(1 - z.^2); 
        else 
            error('Desired dimension not yet implemented!')
        end 
        
    %% sqrt(r)
    elseif strcmp( weightFun, 'sqrt(r)') 
        if dim == 1
            omega = @(x) sqrt( sqrt( x.^2 ) ); 
        elseif dim == 2 
            omega = @(x,y) sqrt( sqrt( x.^2 + y.^2 ) ); 
        elseif dim == 3 
            omega = @(x,y,z) sqrt( sqrt( x.^2 + y.^2 + z.^2 ) ); 
        else 
            error('Desired dimension not yet implemented!')
        end 
        
    %% else
    else 
        error('Desired weight function not yet implemented!') 
    end
end