%% Script to investigate accuracy for the Genz test functions
% The unit cube [0,1]^q is considered! 

%% Setting up the script 
clc, clear 

% free parameters
dim = 2; % dimension (1,2,3)
points = 'Halton'; % points (equid, uniform, Halton) 
CC = 50; % number of tests 

% fixed parameters
domain = 'cube'; % domain (cube, ball) 
volume = 2^dim; 
weightFun = '1'; % weight function (1)

% Set up weight functions 
omega = generate_weightFun( weightFun, dim); 

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

NN_Leg = zeros(n_lenght,1); NN_MC = zeros(n_lenght,1); 
NN_LS = zeros(n_lenght,1); NN_l1 = zeros(n_lenght,1); 
err_Leg = zeros(n_lenght,5); err_MC = zeros(n_lenght,5); 
err_LS = zeros(n_lenght,5); err_l1 = zeros(n_lenght,5);
i = 1;

for c=1:CC 
    
    % Set up Genz test functions 
	a = rand(1,dim); a = a*5/(2*norm(a)); b = rand(1,dim); % random parameters
	if dim == 1  
      	genz1 = @(x) cos( 2*pi*b(1) + a(1)*x ); % osciallatory 
       	genz2 = @(x) ( a(1)^(-2) + (x-b(1)).^2 ).^(-1); % product peak 
       	genz3 = @(x) ( 1 + a(1)*x ).^(-2); % corner peak 
       	genz4 = @(x) exp( -a(1)^2*(x-b(1)).^2 ); % Gaussian
        genz5 = @(x) exp( -a(1)*abs(x-b(1)) ); % C0 function 
        I1 = integral(genz1,0,1,'AbsTol',1e-14); I2 = integral(genz2,0,1,'AbsTol',1e-14); 
        I3 = integral(genz3,0,1,'AbsTol',1e-14); I4 = integral(genz4,0,1,'AbsTol',1e-14); 
        I5 = integral(genz5,0,1,'AbsTol',1e-14); 
  	elseif dim == 2  
      	genz1 = @(x,y) cos( 2*pi*b(1) + a(1)*x + a(2)*y ); % osciallatory 
       	genz2 = @(x,y) ( a(1)^(-2) + (x-b(1)).^2 ).^(-1).*( a(2)^(-2) + (y-b(2)).^2 ).^(-1); % product peak 
       	genz3 = @(x,y) ( 1 + a(1)*x + a(2)*y ).^(-3); % corner peak 
       	genz4 = @(x,y) exp( -a(1)^2*(x-b(1)).^2 - a(2)^2*(y-b(2)).^2 ); % Gaussian 
        genz5 = @(x,y) exp( -a(1)*abs(x-b(1)) - a(2)*abs(y-b(2)) ); % C0 function 
        I1 = integral2(genz1,0,1,0,1,'AbsTol',1e-14); I2 = integral2(genz2,0,1,0,1,'AbsTol',1e-14); 
        I3 = integral2(genz3,0,1,0,1,'AbsTol',1e-14); I4 = integral2(genz4,0,1,0,1,'AbsTol',1e-14); 
        I5 = integral2(genz5,0,1,0,1,'AbsTol',1e-14); 
  	elseif dim == 3  
       	genz1 = @(x,y,z) cos( 2*pi*b(1) + a(1)*x + a(2)*y + a(3)*z ); % osciallatory 
       	genz2 = @(x,y,z) ( a(1)^(-2) + (x-b(1)).^2 ).^(-1).*( a(2)^(-2) + (y-b(2)).^2 ).^(-1).*( a(3)^(-2) + (z-b(3)).^2 ).^(-1); % product peak 
       	genz3 = @(x,y,z) ( 1 + a(1)*x + a(2)*y + a(3)*z).^(-4); % corner peak 
       	genz4 = @(x,y,z) exp( -a(1)^2*(x-b(1)).^2 - a(2)^2*(y-b(2)).^2 - a(3)^2*(z-b(3)).^2); % Gaussian 
        genz5 = @(x,y,z) exp( -a(1)*abs(x-b(1)) - a(2)*abs(y-b(2)) - a(3)*abs(z-b(3))); % C0 function 
        I1 = integral3(genz1,0,1,0,1,0,1,'AbsTol',1e-14); I2 = integral3(genz2,0,1,0,1,0,1,'AbsTol',1e-14); 
        I3 = integral3(genz3,0,1,0,1,0,1,'AbsTol',1e-14); I4 = integral3(genz4,0,1,0,1,0,1,'AbsTol',1e-14); 
        I5 = integral3(genz5,0,1,0,1,0,1,'AbsTol',1e-14); 
    else 
      	error('Desired dimension not yet implemented!') 
    end

    while n <= n_max 
    
        % Legendre rule 
        example = matfile(['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_n=',num2str(n),'.mat']);
        C = example.CF_Leg; 
        [ N, aux] = size(C); 
        NN_Leg(i) = N;
        X = (C(:,1:dim)+1)/2; % data points 
        w = C(:,dim+1)/volume; % weights 
        % Function values 
        if dim == 1  
            genz1_values = genz1(X);
            genz2_values = genz2(X);
            genz3_values = genz3(X);
            genz4_values = genz4(X);
            genz5_values = genz5(X);
        elseif dim == 2  
            genz1_values = genz1(X(:,1),X(:,2)); 
            genz2_values = genz2(X(:,1),X(:,2));
            genz3_values = genz3(X(:,1),X(:,2));
            genz4_values = genz4(X(:,1),X(:,2));
            genz5_values = genz5(X(:,1),X(:,2));
        elseif dim == 3  
            genz1_values = genz1(X(:,1),X(:,2),X(:,3));
            genz2_values = genz2(X(:,1),X(:,2),X(:,3));
            genz3_values = genz3(X(:,1),X(:,2),X(:,3));
            genz4_values = genz4(X(:,1),X(:,2),X(:,3));
            genz5_values = genz5(X(:,1),X(:,2),X(:,3));
        else 
            error('Desired dimension not yet implemented!') 
        end
        % Errors       
        err_Leg(i,1) = err_Leg(i,1) + abs(I1-dot(w,genz1_values)); % absolute error
        err_Leg(i,2) = err_Leg(i,2) + abs(I2-dot(w,genz2_values)); % absolute error
        err_Leg(i,3) = err_Leg(i,3) + abs(I3-dot(w,genz3_values)); % absolute error
        err_Leg(i,4) = err_Leg(i,4) + abs(I4-dot(w,genz4_values)); % absolute error
        err_Leg(i,5) = err_Leg(i,5) + abs(I5-dot(w,genz5_values)); % absolute error

        % LS rule
        example = matfile(['CFs/CF_LS_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
        C = example.CF_LS;       
        NN_LS(i) = N;
        X = (C(:,1:dim)+1)/2; % data points 
        w = C(:,dim+1)/volume; % weights
        % Function values 
        if dim == 1  
            genz1_values = genz1(X);
            genz2_values = genz2(X);
            genz3_values = genz3(X);
            genz4_values = genz4(X);
            genz5_values = genz5(X);
        elseif dim == 2  
            genz1_values = genz1(X(:,1),X(:,2)); 
            genz2_values = genz2(X(:,1),X(:,2));
            genz3_values = genz3(X(:,1),X(:,2));
            genz4_values = genz4(X(:,1),X(:,2));
            genz5_values = genz5(X(:,1),X(:,2));
        elseif dim == 3  
            genz1_values = genz1(X(:,1),X(:,2),X(:,3));
            genz2_values = genz2(X(:,1),X(:,2),X(:,3));
            genz3_values = genz3(X(:,1),X(:,2),X(:,3));
            genz4_values = genz4(X(:,1),X(:,2),X(:,3));
            genz5_values = genz5(X(:,1),X(:,2),X(:,3));
        else 
            error('Desired dimension not yet implemented!') 
        end
        % Errors       
        err_LS(i,1) = err_LS(i,1) + abs(I1-dot(w,genz1_values)); % absolute error
        err_LS(i,2) = err_LS(i,2) + abs(I2-dot(w,genz2_values)); % absolute error
        err_LS(i,3) = err_LS(i,3) + abs(I3-dot(w,genz3_values)); % absolute error
        err_LS(i,4) = err_LS(i,4) + abs(I4-dot(w,genz4_values)); % absolute error
        err_LS(i,5) = err_LS(i,5) + abs(I5-dot(w,genz5_values)); % absolute error

        % l1 rule 
        example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
        C = example.CF_l1; 
        [ N, aux] = size(C); 
        NN_l1(i) = N;
        X = (C(:,1:dim)+1)/2; % data points 
        w = C(:,dim+1)/volume; % weights
        % Function values 
        if dim == 1  
            genz1_values = genz1(X);
            genz2_values = genz2(X);
            genz3_values = genz3(X);
            genz4_values = genz4(X);
            genz5_values = genz5(X);
        elseif dim == 2  
            genz1_values = genz1(X(:,1),X(:,2)); 
            genz2_values = genz2(X(:,1),X(:,2));
            genz3_values = genz3(X(:,1),X(:,2));
            genz4_values = genz4(X(:,1),X(:,2));
            genz5_values = genz5(X(:,1),X(:,2));
        elseif dim == 3  
            genz1_values = genz1(X(:,1),X(:,2),X(:,3));
            genz2_values = genz2(X(:,1),X(:,2),X(:,3));
            genz3_values = genz3(X(:,1),X(:,2),X(:,3));
            genz4_values = genz4(X(:,1),X(:,2),X(:,3));
            genz5_values = genz5(X(:,1),X(:,2),X(:,3));
        else 
            error('Desired dimension not yet implemented!') 
        end
        % Errors       
        err_l1(i,1) = err_l1(i,1) + abs(I1-dot(w,genz1_values)); % absolute error
        err_l1(i,2) = err_l1(i,2) + abs(I2-dot(w,genz2_values)); % absolute error
        err_l1(i,3) = err_l1(i,3) + abs(I3-dot(w,genz3_values)); % absolute error
        err_l1(i,4) = err_l1(i,4) + abs(I4-dot(w,genz4_values)); % absolute error
        err_l1(i,5) = err_l1(i,5) + abs(I5-dot(w,genz5_values)); % absolute error

        % MC integration
        NN_MC(i) = N; 
        % MC weights 
        if dim == 1  
            w = omega( X(:,1) )/N;
       	elseif dim == 2  
           	w = omega( X(:,1), X(:,2) )/N;
      	elseif dim == 3  
          	w = omega( X(:,1), X(:,2), X(:,3) )/N;
        else 
          	error('Desired dimension not yet implemented!') 
        end
        % Function values 
        if dim == 1  
            genz1_values = genz1(X);
            genz2_values = genz2(X);
            genz3_values = genz3(X);
            genz4_values = genz4(X);
            genz5_values = genz5(X);
        elseif dim == 2  
            genz1_values = genz1(X(:,1),X(:,2)); 
            genz2_values = genz2(X(:,1),X(:,2));
            genz3_values = genz3(X(:,1),X(:,2));
            genz4_values = genz4(X(:,1),X(:,2));
            genz5_values = genz5(X(:,1),X(:,2));
        elseif dim == 3  
            genz1_values = genz1(X(:,1),X(:,2),X(:,3));
            genz2_values = genz2(X(:,1),X(:,2),X(:,3));
            genz3_values = genz3(X(:,1),X(:,2),X(:,3));
            genz4_values = genz4(X(:,1),X(:,2),X(:,3));
            genz5_values = genz5(X(:,1),X(:,2),X(:,3));
        else 
            error('Desired dimension not yet implemented!') 
        end
        % Errors       
        err_MC(i,1) = err_MC(i,1) + abs(I1-dot(w,genz1_values)); % absolute error
        err_MC(i,2) = err_MC(i,2) + abs(I2-dot(w,genz2_values)); % absolute error
        err_MC(i,3) = err_MC(i,3) + abs(I3-dot(w,genz3_values)); % absolute error
        err_MC(i,4) = err_MC(i,4) + abs(I4-dot(w,genz4_values)); % absolute error
        err_MC(i,5) = err_MC(i,5) + abs(I5-dot(w,genz5_values)); % absolute error

        % increase i and n 
        i = i+1; 
        if dim == 1 
            n = n + 20;
        elseif dim == 2 
            n = n + 2;
        else 
            n = n + 1;
        end
        
    end
    
end 

% Compute average errors 
err_Leg = err_Leg/CC; err_LS = err_LS/CC; err_l1 = err_l1/CC; err_MC = err_MC/CC;

% Plot - Genz 1 
g = 1;
figure(g) 
p = plot( NN_MC,err_MC(:,g),'ms', NN_LS,err_LS(:,g),'r+', NN_l1,err_l1(:,g),'b^', NN_Leg,err_Leg(:,g),'ko');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ min([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]), max([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]) ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C_N[g_1] - I[g_1]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_Genz/accuracy_Genz1_dim=',num2str(dim),'_',points,'.fig'] );
savefig(str);
str = sprintf( ['plots_Genz/accuracy_Genz1_dim=',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc')

% Plot - Genz 2 
g = 2;
figure(g) 
p = plot( NN_MC,err_MC(:,g),'ms', NN_LS,err_LS(:,g),'r+', NN_l1,err_l1(:,g),'b^', NN_Leg,err_Leg(:,g),'ko');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ min([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]), max([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]) ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C_N[g_2] - I[g_2]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_Genz/accuracy_Genz2_dim=',num2str(dim),'_',points,'.fig'] );
savefig(str);
str = sprintf( ['plots_Genz/accuracy_Genz2_dim=',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc')

% Plot - Genz 3 
g = 3;
figure(g) 
p = plot( NN_MC,err_MC(:,g),'ms', NN_LS,err_LS(:,g),'r+', NN_l1,err_l1(:,g),'b^', NN_Leg,err_Leg(:,g),'ko');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ min([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]), max([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]) ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C_N[g_3] - I[g_3]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_Genz/accuracy_Genz3_dim=',num2str(dim),'_',points,'.fig'] );
savefig(str);
str = sprintf( ['plots_Genz/accuracy_Genz3_dim=',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc')

% Plot - Genz 4 
g = 4;
figure(g) 
p = plot( NN_MC,err_MC(:,g),'ms', NN_LS,err_LS(:,g),'r+', NN_l1,err_l1(:,g),'b^', NN_Leg,err_Leg(:,g),'ko');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
ylim([ min([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]), max([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]) ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C_N[g_4] - I[g_4]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_Genz/accuracy_Genz4_dim=',num2str(dim),'_',points,'.fig'] );
savefig(str);
str = sprintf( ['plots_Genz/accuracy_Genz4_dim=',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc')

% Plot - Genz 5 
g = 5;
figure(g) 
p = plot( NN_MC,err_MC(:,g),'ms', NN_LS,err_LS(:,g),'r+', NN_l1,err_l1(:,g),'b^', NN_Leg,err_Leg(:,g),'ko');
set(p, 'LineWidth',1.5)
set(p, 'markersize',8)
set(gca, 'FontSize', 20)  % Increasing ticks fontsize
xlim([ max([NN_Leg(1);NN_LS(1)]), min([NN_Leg(end);NN_LS(end)]) ]) 
%ylim([ min([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]), max([err_MC(:,g);err_LS(:,g);err_l1(:,g);err_Leg(:,g)]) ])
xlabel('$N$','Interpreter','latex') 
ylabel('$|C_N[g_5] - I[g_5]|$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if strcmp( points, 'Halton')
    id = legend('QMC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
else
    id = legend('MC','LS','$\ell^1$','Legendre','Interpreter','latex','Location','southwest');
end
set(id, 'Interpreter','latex', 'FontSize',26)
grid on
str = sprintf( ['plots_Genz/accuracy_Genz5_dim=',num2str(dim),'_',points,'.fig'] );
savefig(str);
str = sprintf( ['plots_Genz/accuracy_Genz5_dim=',num2str(dim),'_',points,'.eps'] );
%saveas(gcf,str,'epsc')