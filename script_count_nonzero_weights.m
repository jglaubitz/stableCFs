%% Script to count the number of nonzero weights in the l1-CF  

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

NN_l1 = zeros(n_lenght,1); 
dd_l1 = zeros(n_lenght,1);
KK_l1 = zeros(n_lenght,1); 
nonzero_l1 = zeros(n_lenght,1); 
i = 1;

while n <= n_max 
    
    % l1 rule 
    example = matfile(['CFs/CF_l1_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'_n=',num2str(n),'.mat']);
    C = example.CF_l1; 
    [ NN_l1(i), aux] = size(C); 
    dd_l1(i) = C(1,dim+2);
    KK_l1(i) = nchoosek(dim + dd_l1(i), dim);
    % count nonzero weights 
    w = C(:,dim+1); % weights 
    j = 0; 
    for m = 1:NN_l1(i) 
        if w(m) ~= 0 
            j = j+1; 
        end
    end 
    nonzero_l1(i) = j; 
    
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

mmin = min([KK_l1;nonzero_l1]); 
mmax = max([KK_l1;nonzero_l1]);

figure(1) 
p = plot( NN_l1,KK_l1,'k:', NN_l1,nonzero_l1,'r--' ); 
set(p, 'LineWidth',2.5)
set(gca, 'FontSize', 18)  % Increasing ticks fontsize
xlim([ NN_l1(1), NN_l1(end) ]) 
ylim([ mmin-0.1*abs(mmin), mmax+0.1*abs(mmax) ]) 
xlabel('$N$','Interpreter','latex') 
ylabel('$K$ / $\| \mathbf{w}^{\ell^1} \|_0$','Interpreter','latex')
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
id = legend('K','$\| \mathbf{w}^{\ell^1} \|_0$','Interpreter','latex','Location','northwest');
set(id, 'Interpreter','latex', 'FontSize',30)
str = sprintf( ['count_nonzero_weights_dim=',num2str(dim),'_',domain,'_',weightFun,'_',points,'.fig'] );
%savefig(str);