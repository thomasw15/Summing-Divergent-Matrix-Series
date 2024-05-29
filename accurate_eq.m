close all
clear
clc

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

n=1000; I = eye(n); m_max = 100; c = zeros(n); count = 20;

original =zeros(100,1);
oneerror =zeros(100,1);
halferror =zeros(100,1);
threequartererror =zeros(100,1);
oneerrorsp =zeros(100,1);
halferrorsp =zeros(100,1);
threequartererrorsp =zeros(100,1);
timeeq = zeros(100,1);
timesp = zeros(100,1);

parfor k = 1:100
    disp(k);
    olist = zeros(count,1);
    onelist = zeros(count,1);
    halflist = zeros(count,1);
    threequarterlist = zeros(count,1);
    onelistsp = zeros(count,1);
    halflistsp = zeros(count,1);
    threequarterlistsp = zeros(count,1);

    teq = zeros(count,1);
    tsp = zeros(count,1);
    for j = 1:count
        c = zeros(n);
        [V,W] = bidiagonal(n,1-k/100);
        sum = I; prod = I;
        for i = 1: m_max
            prod = double(prod * V);
            y = double(prod - c);
            t = double(sum + y);
            c = double((t - sum) - y); 
            sum = double(t);
        end
        
        originalerror = norm(W-sum);
        olist(j) = log(originalerror);
        
        %compensated
        tic;
        approx = eqsum_simplified(m_max, 1, V);
        error = norm(W-approx);
        onelist(j) = log(error);

        approx = eqsum_simplified(m_max, 0.75, V);
        error = norm(W-approx);
        threequarterlist(j) = log(error);

        approx = eqsum_simplified(m_max, 0.5, V);
        error = norm(W-approx);
        halflist(j) = log(error);
        teq(j) = toc;
        
        %Schur-Parlett
        tic;
        approx = SchurParlett(V, 0.1, m_max, 0.5);
        error = norm(W-approx);
        halflistsp(j) = log(error);

        approx = SchurParlett(V, 0.1, m_max, 1);
        error = norm(W-approx);
        onelistsp(j) = log(error);

        approx = SchurParlett(V, 0.1, m_max, 0.75);
        error = norm(W-approx);
        threequarterlistsp(j) = log(error);
        tsp(j) = toc;
    
    end
    original(k) = mean(olist);
    oneerror(k) = mean(onelist);
    halferror(k) = mean(halflist);
    threequartererror(k) = mean(threequarterlist);
    oneerrorsp(k) = mean(onelistsp);
    halferrorsp(k) = mean(halflistsp);
    threequartererrorsp(k) = mean(threequarterlistsp);

    timeeq(k) = mean(teq);
    timesp(k) = mean(tsp);
end

disp('end');

filename = 'accurate_eq.mat';
save(filename)

fig1 = figure(1);
hold on
plot(original,'LineWidth',1)
plot(oneerror,'LineWidth',1)
plot(halferror,'LineWidth',1)
plot(threequartererror,'LineWidth',1)
%title('Error from Eq summation','Interpreter','latex','Fontsize',16)
leg1 = legend('ordinary', '$\rho=1$','$\rho=1/2$','$\rho=3/4$','Fontsize',16,'Interpreter','latex');
set(leg1,'Interpreter','latex');
ylabel('average log error','Fontsize',16,'Interpreter','latex')
xlabel('Percentage','Fontsize',16,'Interpreter','latex')
hold off


fig2 = figure(2);
hold on
plot(original,'LineWidth',1)
plot(oneerrorsp,'LineWidth',1)
plot(halferrorsp,'LineWidth',1)
plot(threequartererrorsp,'LineWidth',1)
%title('Error from Eq summation with Schur Parlett','Interpreter','latex','Fontsize',16)
leg1 = legend('ordinary', '$\rho=1$','$\rho=1/2$','$\rho=3/4$','Fontsize',16,'Interpreter','latex');
set(leg1,'Interpreter','latex');
ylabel('average log error','Fontsize',16,'Interpreter','latex')
xlabel('Percentage','Fontsize',16,'Interpreter','latex')
hold off

fig3 = figure(3);
hold on
plot(timeeq,'LineWidth',1)
plot(timesp,'LineWidth',1)
title('Time of normal vs Schur Parlett','Interpreter','latex','Fontsize',16)
leg1 = legend('normal', 'Schur-Parlett','Fontsize',16,'Interpreter','latex');
set(leg1,'Interpreter','latex');
ylabel('average log error','Fontsize',16,'Interpreter','latex')
xlabel('Percentage','Fontsize',16,'Interpreter','latex')
hold off

print(fig1,'accurate','-dpdf')
print(fig2,'accurate_sp','-dpdf')
print(fig3,'accurate_Time','-dpdf')

