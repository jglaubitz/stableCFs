%% Script to investigate accuracy 

%% Setting up the script 
clc, clear 

% free parameters
points = 'Halton'; % points (equid, uniform, Halton)
n_start = 10; 
n_step = 10;
n_max = 80; 

% fixed parameters
domain = 'combi'; % domain  
dim = 2; % dimension 
volume = pi + 1;
weightFun = '1'; % weight function - 1, C2k, sqrt(r)

% Test and weight function as well as the referecne value 
omega = generate_weightFun( weightFun, dim); 
f = @(x,y) exp( - x.^2 - y.^2 ); % test function 
I1 = pi*( 1 - exp(-1) );
I2 = integral2(f,1,2,1,2,'AbsTol',1e-14);
I = I1 + I2; % exact integral 


NN_MC = []; NN_LS = []; NN_l1 = []; % number of data points 
err_MC = []; err_LS = []; err_l1 = []; % errors 

for n=n_start:n_step:n_max     
    
    % LS rule
    example = matfile(['CFs/CF_LS_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS;
    [ N, aux] = size(C); 
    NN_LS = [NN_LS; N];
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = f(X(:,1),X(:,2));
    err_LS = [err_LS; abs( I - dot(w,f_values) )]; % absolute error
    
    % l1 rule 
    example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_l1; 
    [ N, aux] = size(C); 
    NN_l1 = [NN_l1; N]; 
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = f(X(:,1),X(:,2));
    err_l1 = [err_l1; abs( I - dot(w,f_values) )]; % absolute error
    
    % MC integration
    NN_MC = [NN_MC; N]; 
    % MC weights 
    w = volume*omega(X(:,1),X(:,2))/N;
    CF = dot( w, f_values ); % value of the CF 
    err_MC = [err_MC; abs( I - dot(w,f_values) )]; % absolute error
    
end

% Plot the results
figure(1) 
p = plot( NN_MC,err_MC,'ms', NN_LS,err_LS,'r+', NN_l1,err_l1,'b^');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ NN_LS(1), NN_LS(end) ]) 
ylim([ min([err_MC;err_LS;err_l1])/10, max([err_MC;err_LS;err_l1])*10 ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C[f] - I[f]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_nonstandard_domain/nonstandard_',points,'.fig'] );
%savefig(str);