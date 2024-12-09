close all
clear
clc

%pc = parcluster('local')
%parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

%Pade approximation of the binomial series


d=10;
converror = zeros(45,1);
eulererror = zeros(45,1);
convcond = zeros(45,1);
eulercond = zeros(45,1);

for k= 5:45
    error1 = zeros(10,1);
    error2 = zeros(10,1);
    cond1 = zeros(10,1);
    cond2 = zeros(10,1);
    for j = 1:10
        m=k;
        n=k;
        
        Y = randn(d);
        [Q, ~] = qr(Y);
        v = 0.5*rand(d,1)-0.25;
        Z = triu(0.25*randn(d));
        Z = Z- diag(diag(Z)) + diag(v);
        X = double(Q * Z * ctranspose(Q));
        II = eye(d);
        row = zeros(n+1,1);
        a = zeros(m+n+1,1);
        a(1) = 1;
        alpha = -3/4;
        %alpha = 1/4;
        %alpha = -1/4;
        %alpha = 3/4;
        %alpha = -3/5; (31,31) approximant does not exist
        %alpha = 3/5;
        %alpha = -1/2;
        %alpha = 1/2;
        %alpha = 11/7;
        %alpha = -11/7;
        coef = alpha;
        
        for i = 2:m+n+1
            a(i) = double(coef); %This is binomial series (1+x)^alpha
            coef = double(coef * (alpha-i+1)/i);
        end
        row(1) = a(1);
        
        %regular sum
        R = double(toeplitz(a(m+1:m+n),a(m+1:-1:m-n+2)));
        B = double(a(m+2:m+n+1));
        if rank(R) ~= rank([R B])
            disp("standard Pade is ill-conditioned");
            %disp(k);
            %disp(R);
            break
        end
        gamma = double(R\(-B));
        gamma = double([1, gamma']');
        beta = double(toeplitz(a(1:m+1),row) * gamma);
        
        P = double(beta(1)*eye(d) + beta(2)*X); 
        Q = double(gamma(1)*eye(d) + gamma(2)*X);  
        power = double(X);
        C1 = double(zeros(d));
        C2 = double(zeros(d));
        for i = 2:max([m n])
            power = double(power*X);
            if i <= m+1
                Y1 = double(a(i+1)*power - C1);
                T1 = double(P + Y1);
                C1 = double(double(T1-P)-Y1);
                P = T1;
            end    
            if i <= n+1
                Y2 = double(a(i+1)*power - C2);
                T2 = double(Q + Y2);
                C2 = double(double(T2-Q)-Y2);
                Q = T2;
            end   
        end
        conventional = double(P*double(Q\eye(d)));
        actual = double((II+X)^alpha);
        error1(j) = norm(double(conventional-actual),2);
        cond1(j) = cond(Q);

        %Euler sum
        p = 10;
        dE = zeros(m+n+1,1);
        cE = zeros(m+n+1,1);
        for i=1:m+n+1
            cE(i) = double(nchoosek(m+n+1,i)*p^(m+n-i+1)/(1+p)^(m+n+1));
        end
        for i=1:m+n+1
            dE(i) = double(a(i) * sum(cE(i:m+n+1)));
        end
        row = zeros(n+1,1);
        row(1) = dE(1);
        
        R = toeplitz(dE(m+1:m+n),dE(m+1:-1:m-n+2));
        B = dE(m+2:m+n+1);
        if rank(R) ~= rank([R B])
            disp("Euler--Pade is ill-conditioned");
            break
            %disp(R);
        end
        gamma = R\(-B);
        gamma = [1, gamma']';
        beta = toeplitz(dE(1:m+1),row) * gamma;
        
        P = double(beta(1)*eye(d) + beta(2)*X); 
        Q = double(gamma(1)*eye(d) + gamma(2)*X);  
        power = double(X);
        C1 = double(zeros(d));
        C2 = double(zeros(d));
        for i = 2:max([m n])
            power = double(power*X);
            if i <= m+1
                Y1 = double(a(i+1)*power - C1);
                T1 = double(P + Y1);
                C1 = double(double(T1-P)-Y1);
                P = T1;
            end    
            if i <= n+1
                Y2 = double(a(i+1)*power - C2);
                T2 = double(Q + Y2);
                C2 = double(double(T2-Q)-Y2);
                Q = T2;
            end   
        end
        euler = double(P*double(Q\eye(d)));
        error2(j) = norm(euler-actual,2);
        cond2(j) = cond(Q);
    end
    converror(k) = mean(error1);
    eulererror(k) = mean(error2);
    convcond(k) = mean(cond1);
    eulercond(k) = mean(cond2);
end

fig = figure(1);
hold on
plot(log10(converror),'LineWidth',1)
plot(log10(eulererror),'LineWidth',1)
ylim([-1 8])
leg1 = legend("Pad\'e","Pad\'e--Euler",'Fontsize',16,'Interpreter','latex','location','southeast');
set(leg1,'Interpreter','latex');
ylabel('log forward error','Fontsize',16,'Interpreter','latex')
xlabel('$n$','Fontsize',16,'Interpreter','latex')
hold off

fig2 = figure(2);
hold on
plot(log10(convcond),'LineWidth',1)
plot(log10(eulercond),'LineWidth',1)
ylim([-1 10])
leg1 = legend("Pad\'e","Pad\'e--Euler",'Fontsize',16,'Interpreter','latex','location','southeast');
set(leg1,'Interpreter','latex');
ylabel('log condition number of $Q$','Fontsize',16,'Interpreter','latex')
xlabel('$n$','Fontsize',16,'Interpreter','latex')
hold off

%print(fig,'pade_binomial_k','-dpdf')
