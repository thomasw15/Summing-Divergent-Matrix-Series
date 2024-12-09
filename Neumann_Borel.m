close all
clear
clc

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

n=1000;
q=1000;
trials = 50;
m_max = 10^8;
logerrors_StrongBorel = zeros(trials,1);
logerrors_eq = zeros(trials,1);
logerrors_inv = zeros(trials,1);
forwarderrors_StrongBorel = zeros(trials,1);

parfor i = 1:trials
    disp(i);
    
    Y = randn(n) + 1i * randn(n);
    [Q, ~] = qr(Y);
    T2 = 100*rand(n)-50;
    T3 = 1i*(100*rand(n)-50);
    T2 = T2 + T3;
    v1 = 10^3*rand(n,1)-(10^3+100);
    v2 = 1i*(10^3*rand(n,1)-500);
    v = v1+v2;
    Z = triu(T2);
    Z = Z- diag(diag(Z)) + diag(v);
    Z = Q * Z * ctranspose(Q);
    II = eye(n);
        
    %calculate the Strong Borel summation
    disp("Borel");
    tic
    M_StrongBorel = StrongBorel(Z);
    Borel_error = norm((double(M_StrongBorel)*double(II-Z)-double(II)),2);
    logerrors_StrongBorel(i) = log(Borel_error);
    toc

    inverse = double(inv(double(II-Z)));
    logerrors_inv(i) = log(norm((inverse*double(II-Z)-double(II)),2));
    
    forwarderrors_StrongBorel(i) = norm(double(M_StrongBorel)-inverse,2);
end

filename = 'Neumann_Borel_1114.mat';
save(filename)

hold on
fig = figure(1);
plot(logerrors_StrongBorel,'LineWidth',1)
plot(logerrors_inv,'LineWidth',1)
title('Error from the canonical matrix inversion','Interpreter','latex','Fontsize',16)
leg1 = legend('Strong Borel summation', 'inverse function');
set(leg1,'Interpreter','latex');
ylabel('log error','Fontsize',16,'Interpreter','latex')
xlabel('trials','Fontsize',16,'Interpreter','latex')
hold off

fig2 = figure(2);
hold on
plot(log10(forwarderrors_StrongBorel),'LineWidth',1)
%title('Error from the canonical matrix inversion','Interpreter','latex','Fontsize',16)
leg1 = legend('Strong Borel','Fontsize',16,'Interpreter','latex','location','northwest');
set(leg1,'Interpreter','latex');
ylabel('log forward error','Fontsize',16,'Interpreter','latex')
xlabel('trials','Fontsize',16,'Interpreter','latex')
hold off

print(fig,'Neumann_Borel_1114','-dpdf')
print(fig2,'Neumann_Borel_forward_1114','-dpdf')

