close all
clear
clc

pc = parcluster('local')
parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

n=1000;
q=double(10^4);
trials = 50;
m_max = 10^4;
logerrors_eq = zeros(trials,1);
logerrors_inv = zeros(trials,1);
logerrors_sp = zeros(trials,1);
timeeq = zeros(trials,1);
timesp = zeros(trials,1);
forwarderrors_eq = zeros(trials,1);
forwarderrors_sp = zeros(trials,1);

parfor i = 1:trials
    disp(i);
    
    Y = randn(n) + 1i * randn(n);
    [Q, ~] = qr(Y);
    T2 = double(100*rand(n)-50);
    T3 = double(1i*(100*rand(n)-50));
    T2 = T2 + T3;
    v1 = double(10^3*rand(n,1)-(10^3+100));
    v2 = double(1i*(10^3*rand(n,1)-500));
    v = v1+v2;
    Z = triu(T2);
    Z = double(Z- diag(diag(Z)) + diag(v));
    Z = Q * Z * ctranspose(Q);
    II = double(eye(n));
    tic
    M_eq = eqsum_simplified(m_max, q, Z);
    M_eq = double(M_eq);
    eq_error = norm((M_eq*(II-Z)-II),2);
    logerrors_eq(i) = log(eq_error);
    timeeq(i) = toc;

    inverse = inv((II-Z));
    logerrors_inv(i) = log(norm((inverse*(II-Z)-II),2));
    
    tic
    F = SchurParlett(Z, 0.1, m_max, q);
    logerrors_sp(i) = log(norm((F*(II-Z)-II),2));
    timesp(i) = toc;

    forwarderrors_eq(i) = norm(M_eq-inverse,2);
    forwarderrors_sp(i) = norm(F-inverse,2);
end

disp("end");

filename = 'Neumann_eq_0527.mat';
save(filename)

fig = figure(1);
hold on
plot(logerrors_eq,'LineWidth',1)
plot(logerrors_sp,'LineWidth',1)
plot(logerrors_inv,'LineWidth',1)
leg1 = legend('recursive Euler', 'Schur-Parlett Euler','inverse function');
set(leg1,'Interpreter','latex');
ylabel('log error','Fontsize',16,'Interpreter','latex')
xlabel('trials','Fontsize',16,'Interpreter','latex')
hold off

fig2 = figure(2);
hold on
plot(timeeq,'LineWidth',1)
plot(timesp,'LineWidth',1)
leg1 = legend('recursive Euler', 'Schur-Parlett Euler','Fontsize',16,'Interpreter','latex');
set(leg1,'Interpreter','latex');
ylabel('time','Fontsize',16,'Interpreter','latex')
xlabel('trials','Fontsize',16,'Interpreter','latex')
hold off

fig3 = figure(3);
hold on
plot(log10(forwarderrors_eq),'LineWidth',1)
plot(log10(forwarderrors_sp),'LineWidth',1)
leg1 = legend('recursive Euler','Schur-Parlett Euler','Fontsize',16,'Interpreter','latex','location','northwest');
set(leg1,'Interpreter','latex');
ylabel('log forward error','Fontsize',16,'Interpreter','latex')
xlabel('trials','Fontsize',16,'Interpreter','latex')
hold off

print(fig,'Neumann_eq_sp_0527','-dpdf')
print(fig2,'SchurParlett_Time_0527','-dpdf')
print(fig3,'Neumann_eq_sp_forward','-dpdf')
