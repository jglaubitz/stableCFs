%% Script to plot and compute the ratio between N, d, and K  

%% Setting up the script 
clc, clear 

dim = 3; % dimension (1,2,3)
domain = 'ball'; % domain (cube, ball) 
weightFun = 'sqrt(r)'; % weight function - 1, C2k, sqrt(r)
points = 'Halton'; % points (equid, semi-uniform, uniform, Halton) 

if dim == 1 
    n = 20;
    n_max = 400; 
elseif dim == 2 
    n = 4;
    n_max = 40; 
else 
    n = 4;
    n_max = 16; 
end

NN_Leg = []; NN_LS = []; NN_l1 = []; % number of data points 
dd_Leg = []; dd_LS = []; dd_l1 = []; % degree of exactness 
KK_Leg = []; KK_LS = []; KK_l1 = []; % corresponding number of basis functions 

while n <= n_max 
    
    % Legendre rule 
    example = matfile(['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_n=',num2str(n),'.mat']);
    C = example.CF_Leg; 
    [ N, aux] = size(C); 
    NN_Leg = [NN_Leg; N]; 
    dd_Leg = [dd_Leg; C(1,dim+2); ]; 
    KK_Leg = [KK_Leg; C(2,dim+2); ];
    
    % LS rule
    example = matfile(['CFs/CF_LS_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_LS; 
    [ N, aux] = size(C); 
    NN_LS = [NN_LS; N]; 
    dd_LS = [dd_LS; C(1,dim+2); ]; 
    KK_LS = [KK_LS; C(2,dim+2); ];
    
    % l1 rule 
    example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_l1; 
    [ N, aux] = size(C); 
    NN_l1 = [NN_l1; N]; 
    dd_l1 = [dd_l1; C(1,dim+2); ]; 
    KK_l1 = [KK_l1; C(2,dim+2); ];
    
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
    p = plot( NN_Leg,KK_Leg,'ko', NN_LS,KK_LS,'r+', NN_l1,KK_l1,'b^');
    set(p, 'LineWidth',1.5)
    set(p, 'markersize',8)
    set(gca, 'FontSize', 20)  % Increasing ticks fontsize
    xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
    ylim([ min([KK_Leg;KK_LS;KK_l1])-1, max([KK_Leg;KK_LS;KK_l1])+20 ])
    xlabel('$N$','Interpreter','latex') 
    ylabel('$K$','Interpreter','latex')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    id = legend('Legendre','LS','$\ell^1$','Interpreter','latex','Location','northwest');
    set(id, 'Interpreter','latex', 'FontSize',26)
    grid on
    str = sprintf( ['ratio_plots/ratio_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'.fig'] );
    %savefig(str);
else 
    figure(1) 
    p = plot( NN_LS,KK_LS,'r+', NN_l1,KK_l1,'b^');
    set(p, 'LineWidth',1.5)
    set(p, 'markersize',8)
    set(gca, 'FontSize', 20)  % Increasing ticks fontsize
    xlim([ NN_LS(1), NN_LS(end) ]) 
    ylim([ min([KK_LS;KK_l1])-1, max([KK_LS;KK_l1])+20 ])
    xlabel('$N$','Interpreter','latex') 
    ylabel('$K$','Interpreter','latex')
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    id = legend('LS','$\ell^1$','Interpreter','latex','Location','northwest');
    set(id, 'Interpreter','latex', 'FontSize',26)
    grid on
    str = sprintf( ['ratio_plots/ratio_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'.fig'] );
    %savefig(str);
end

if strcmp( domain, 'cube') && strcmp( weightFun, '1')
    % LS fit for the parameters s and C - Legendre 
    f_Leg = fittype('a*x.^b'); % set up model for parameters 
    xdata = NN_Leg; 
    ydata = KK_Leg; 
    fit( ydata, xdata, f_Leg) 
end

% LS fit for the parameters s and C - LS 
f_LS = fittype('a*x.^b'); % set up model for parameters 
xdata = NN_LS; 
ydata = KK_LS; 
fit( ydata, xdata, f_LS, 'Lower', [0,0] )

% LS fit for the parameters s and C - l1
f_l1 = fittype('a*x.^b'); % set up model for parameters 
xdata = NN_l1; 
ydata = KK_l1; 
fit( ydata, xdata, f_l1, 'Lower', [0,0] )