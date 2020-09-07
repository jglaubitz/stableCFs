%% Script to investigate accuracy 

%% Setting up the script 
clc, clear 

% free parameters
dim = 2; % dimension (1,2,3)
points = 'equid'; % points (equid, uniform, Halton)
noise_level = 0; % 0, 10^(-6)

% fixed parameters
domain = 'ball'; % domain (cube, ball) 
weightFun = '1'; % weight function - 1, C2k, sqrt(r)

omega = generate_weightFun( weightFun, dim); 

if dim == 1 
    f = @(x) 1./(1 + x.^2 ) + sin(x); % test function 
    I_aux = atan(1); 
    n = 20;
    n_max = 400; 
    n_lenght = 20; 
elseif dim == 2 
    f = @(x,y) 1./(1 + x.^2 + y.^2 ) + sin(x); % test function 
    I_aux = 0.5*log(2);
    n = 4;
    n_max = 40; 
    n_lenght = 19; 
else 
    f = @(x,y,z) 1./(1 + x.^2 + y.^2 + z.^2 ) + sin(x); % test function 
    I_aux = 1 - atan(1);
    n = 4;
    n_max = 16; 
    n_lenght = 13;
end 
I = I_aux*( 2*pi^(dim/2)/gamma(dim/2) ); % exact integral

NN_Leg = zeros(n_lenght,1); 
err_Leg = zeros(n_lenght,1);
NN_MC = zeros(n_lenght,1); 
err_MC = zeros(n_lenght,1);
NN_LS = zeros(n_lenght,1); 
err_LS = zeros(n_lenght,1);
NN_LS2 = zeros(n_lenght,1); 
err_LS2 = zeros(n_lenght,1);
NN_l1 = zeros(n_lenght,1); 
err_l1 = zeros(n_lenght,1);
i = 1;

while n <= n_max 
    
    % Legendre rule 
    example = matfile(['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_n=',num2str(n),'.mat']);
    C = example.CF_Leg; 
    [ M, aux] = size(C); % number of data points 
    NN_Leg(i) = M;
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(M,1); 
    for m = 1:M 
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
    noise = noise_level*(2*rand(M,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_Leg(i) = abs( I - CF ); % absolute error
    
    % LS rule - dir
    example = matfile(['CFs/CF_LS_dir_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS_dir;
    [ M, aux] = size(C); % number of data points 
    NN_LS(i) = M;
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(M,1); 
    for m = 1:M 
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
    noise = noise_level*(2*rand(M,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_LS(i) = abs( I - CF ); % absolute error
    
    % LS rule - opt
    example = matfile(['CFs/CF_LS_opt_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS_opt;
    [ M, aux] = size(C); % number of data points 
    NN_LS2(i) = M; 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(M,1); 
    for m = 1:M 
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
    noise = noise_level*(2*rand(M,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_LS2(i) = abs( I - CF ); % absolute error
    
    % l1 rule 
    example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_l1; 
    [ M, aux] = size(C); 
    NN_l1(i) = M; 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(M,1); 
    for m = 1:NN_l1(i) 
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
    noise = noise_level*(2*rand(M,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_l1(i) = abs( I - CF ); % absolute error
    
    % MC integration
    NN_MC(i) = M; 
    % MC weights 
    for m = 1:M 
       if dim == 1  
            w(m) = 2*omega( X(m,1) )/M; 
       elseif dim == 2  
            w(m) = 4*omega( X(m,1), X(m,2) )/M;
       elseif dim == 3  
            w(m) = 8*omega( X(m,1), X(m,2) , X(m,3) )/M;
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    CF = dot( w, f_values ); % value of the CF 
    err_MC(i) = abs( I - CF ); % absolute error
    
    % increase n
    if dim == 1 
        n = n + 20;
    elseif dim == 2 
        n = n + 2;
    else 
        n = n + 1;
    end
    i = i+1;
    
end 

% plot the results
figure(1) 
p = plot( NN_MC,err_MC,'m-.', NN_LS2,err_LS2,'r-', NN_l1,err_l1,'b--', NN_Leg,err_Leg,'k:' ); 
set(p, 'LineWidth',2.5)
set(gca, 'FontSize', 18)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS2(1);NN_l1(1)]), min([NN_Leg(end);NN_LS2(end);NN_l1(end)]) ]) 
ylim([ min([err_LS2;err_l1;err_MC;err_Leg]), max([err_LS2;err_l1;err_MC;err_Leg]) ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C[f] - I[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','best');
set(id, 'Interpreter','latex', 'FontSize',22)
str = sprintf( ['accuracy_test3_dim=',num2str(dim),'_',points,'.fig'] );
%savefig(str); 