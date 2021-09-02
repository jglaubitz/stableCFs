%% Script to investigate accuracy 

%% Setting up the script 
clc, clear 

% free parameters
dim = 3; % dimension (1,2,3)
points = 'Halton'; % points (equid, uniform, Halton)
noise_level = 0; % 0, 10^(-6)

% fixed parameters
domain = 'ball'; % domain (cube, ball) 
weightFun = 'sqrt(r)'; % weight function - 1, C2k, sqrt(r)

omega = generate_weightFun( weightFun, dim); 

if dim == 1 
    volume = 2;
    f = @(x) 1./(1 + x.^2 ) + sin(x); % test function 
    I_aux = -2^(-3/2)*( log(2+sqrt(2)) - log(2-sqrt(2)) - 2*atan(sqrt(2)+1) - 2*atan(sqrt(2)-1) ); 
    n = 20;
    n_max = 400; 
elseif dim == 2 
    volume = pi;
    f = @(x,y) 1./(1 + x.^2 + y.^2 ) + sin(x); % test function 
    I_aux = 2^(-3/2)*( -log(2+sqrt(2)) + log(2-sqrt(2)) + 2^(5/2) - 2*atan(sqrt(2)+1) - 2*atan(sqrt(2)-1) ); 
    n = 4;
    n_max = 40;  
else 
    volume = 4*pi/3;
    f = @(x,y,z) 1./(1 + x.^2 + y.^2 + z.^2 ) + sin(x); % test function 
    I_aux = 2^(-3/2)*( log(2+sqrt(2)) - log(2-sqrt(2)) ) + 2/3 -2^(-1/2)*( atan(sqrt(2)+1) + atan(sqrt(2)-1) ); 
    n = 4;
    n_max = 16; 
end 
I = I_aux*( 2*pi^(dim/2)/gamma(dim/2) ); % exact integral

NN_Leg = []; NN_MC = []; NN_LS = []; NN_l1 = []; % number of data points 
err_Leg = []; err_MC = []; err_LS = []; err_l1 = []; % errors 

while n <= n_max 
    
    % Legendre rule 
    example = matfile(['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_n=',num2str(n),'.mat']);
    C = example.CF_Leg; 
    [ N, aux] = size(C); 
    NN_Leg = [NN_Leg; N];
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(N,1); 
    for m = 1:N 
       if dim == 1  
            f_values(m) = f( X(m,1) ).*omega( X(m,1) ); 
       elseif dim == 2  
            f_values(m) = f( X(m,1), X(m,2) ).*omega( X(m,1), X(m,2) );
       elseif dim == 3  
            f_values(m) = f( X(m,1), X(m,2) , X(m,3) ).*omega( X(m,1), X(m,2) , X(m,3) );
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    % generate and add uniform noise 
    noise = noise_level*(2*rand(N,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_Leg = [err_Leg; abs( I - CF )]; % absolute error
    
    % LS rule
    example = matfile(['CFs/CF_LS_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS;
    [ N, aux] = size(C); 
    NN_LS = [NN_LS; N];
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(N,1); 
    for m = 1:N 
       if dim == 1  
            f_values(m) = f( X(m,1) ); 
       elseif dim == 2  
            f_values(m) = f( X(m,1), X(m,2) );
       elseif dim == 3  
            f_values(m) = f( X(m,1), X(m,2) , X(m,3) );
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    % generate and add uniform noise 
    noise = noise_level*(2*rand(N,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_LS = [err_LS; abs( I - CF )]; % absolute error
    
    % l1 rule 
    example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_l1; 
    [ N, aux] = size(C); 
    NN_l1 = [NN_l1; N]; 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(N,1); 
    for m = 1:N 
       if dim == 1  
            f_values(m) = f( X(m,1) ); 
       elseif dim == 2  
            f_values(m) = f( X(m,1), X(m,2) );
       elseif dim == 3  
            f_values(m) = f( X(m,1), X(m,2) , X(m,3) );
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    % generate and add uniform noise 
    noise = noise_level*(2*rand(N,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_l1 = [err_l1; abs( I - CF )]; % absolute error
    
    % MC integration
    NN_MC = [NN_MC; N]; 
    % MC weights 
    for m = 1:N 
       if dim == 1  
            w(m) = volume*omega( X(m,1) )/N; 
       elseif dim == 2  
            w(m) = volume*omega( X(m,1), X(m,2) )/N;
       elseif dim == 3  
            w(m) = volume*omega( X(m,1), X(m,2) , X(m,3) )/N;
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    CF = dot( w, f_values ); % value of the CF 
    err_MC = [err_MC; abs( I - CF )]; % absolute error
    
    % increase n
    if dim == 1 
        n = n + 20;
    elseif dim == 2 
        n = n + 2;
    else 
        n = n + 1;
    end
    
end

% Plot the results
figure(1) 
p = plot( NN_MC,err_MC,'ms', NN_LS,err_LS,'r+', NN_l1,err_l1,'b^', NN_Leg,err_Leg,'ko');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
ylim([ min([err_MC;err_LS;err_l1;err_Leg])/10, max([err_MC;err_LS;err_l1;err_Leg])*10 ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C[f] - I[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_accuracy/accuracy_test2_dim=',num2str(dim),'_',points,'.fig'] );
%savefig(str);