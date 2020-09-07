%% Script to plot and compute the ratio between N, d, and K  

%% Setting up the script 
clc, clear 

dim = 2; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
weightFun = 'C2k'; % weight function - 1, C2k, sqrt(r)
points = 'equid'; % points (equid, semi-uniform, uniform, Halton) 

if dim == 1 
    n = 20;
    n_max = 400; 
    n_lenght = 20;
elseif dim == 2 
    n = 4;
    n_max = 40; 
    n_lenght = 19;
else 
    n = 4;
    n_max = 16; 
    n_lenght = 13;
end

NN_Leg = zeros(n_lenght,1); 
dd_Leg = zeros(n_lenght,1);
NN_LS = zeros(n_lenght,1); 
dd_LS = zeros(n_lenght,1);
NN_LS2 = zeros(n_lenght,1); 
dd_LS2 = zeros(n_lenght,1);
NN_l1 = zeros(n_lenght,1); 
dd_l1 = zeros(n_lenght,1);
i = 1;

while n <= n_max 
    
    % Legendre rule 
    example = matfile(['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_n=',num2str(n),'.mat']);
    C = example.CF_Leg; 
    [ NN_Leg(i), aux] = size(C); 
    dd_Leg(i) = C(1,dim+2);
    
    % LS rule - dir
    example = matfile(['CFs/CF_LS_dir_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS_dir; 
    [ NN_LS(i), aux] = size(C); 
    dd_LS(i) = C(1,dim+2);
    
    % LS rule - opt
    example = matfile(['CFs/CF_LS_opt_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS_opt; 
    [ NN_LS2(i), aux] = size(C); 
    dd_LS2(i) = C(1,dim+2);
    
    % l1 rule 
    example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_l1; 
    [ NN_l1(i), aux] = size(C); 
    dd_l1(i) = C(1,dim+2);
    
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

if strcmp( domain, 'cube') && strcmp( weightFun, '1')
    figure(1) 
    p = plot( NN_Leg,dd_Leg,'k:', NN_LS2,dd_LS2,'r-', NN_l1,dd_l1,'b--' ); 
    set(p, 'LineWidth',2.5)
    set(gca, 'FontSize', 18)  % Increasing ticks fontsize
    xlim([ max([NN_Leg(1);NN_LS2(1);NN_l1(1)]), min([NN_Leg(end);NN_LS2(end);NN_l1(end)]) ]) 
    ylim([ min([dd_Leg;dd_LS2;dd_l1])-1, max([dd_Leg;dd_LS2;dd_l1])+1 ])
    xlabel('$N$','Interpreter','latex') 
    ylabel('$d$','Interpreter','latex')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    id = legend('Legendre','LS','$\ell^1$','Interpreter','latex','Location','northwest');
    set(id, 'Interpreter','latex', 'FontSize',26)
    str = sprintf( ['ratio_plots/ratio_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'.fig'] );
    savefig(str);
else 
    figure(1) 
    p = plot( NN_LS2,dd_LS2,'r-', NN_l1,dd_l1,'b--' ); 
    set(p, 'LineWidth',2.5)
    set(gca, 'FontSize', 18)  % Increasing ticks fontsize
    xlim([ max([NN_LS2(1);NN_l1(1)]), min([NN_LS2(end);NN_l1(end)]) ]) 
    ylim([ min([dd_LS2;dd_l1])-1, max([dd_LS2;dd_l1])+1 ]) 
    xlabel('$N$','Interpreter','latex') 
    ylabel('$d$','Interpreter','latex')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    id = legend('LS','$\ell^1$','Interpreter','latex','Location','northwest');
    set(id, 'Interpreter','latex', 'FontSize',26)
    str = sprintf( ['ratio_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'.fig'] );
    %savefig(str);
end

if strcmp( domain, 'cube') && strcmp( weightFun, '1')
    % LS fit for the parameters s and C - Legendre 
    f_Leg = fittype('a*x.^b'); % set up model for parameters 
    xdata = NN_Leg; 
    ydata = dd_Leg; 
    fit( ydata, xdata, f_Leg) 
end

% LS fit for the parameters s and C - LS 
f_LS = fittype('a*x.^b'); % set up model for parameters 
xdata = NN_LS2; 
ydata = dd_LS2; 
fit( ydata, xdata, f_LS, 'Lower', [0,0] )

% LS fit for the parameters s and C - l1
f_l1 = fittype('a*x.^b'); % set up model for parameters 
xdata = NN_l1; 
ydata = dd_l1; 
fit( ydata, xdata, f_l1, 'Lower', [0,0] )