close all
clear
clc

%Pade approximation of the binomial series

d=10;
eulererror = zeros(20,1);
Cerror = zeros(20,1);
p = 100;
alpha = -3/4;
%alpha = 1/4;
%alpha = -1/4;
%alpha = 3/4;
%alpha = -3/5; 
%alpha = 3/5;
%alpha = -1/2;
%alpha = 1/2;
%alpha = 4/7;
%alpha = -4/7;

r = 150;
for k= 1:20
    disp(k);
    error1 = zeros(10,1);
    error2 = zeros(10,1);
    cond1 = zeros(10,1);
    cond2 = zeros(10,1);
    errorC = zeros(10,1);
    condC = zeros(10,1);
    for j = 1:10
        m=k;
        n=k;
        
        Y = randn(d);
        [Q, ~] = qr(Y);
        v = 3*r/4+r/4*rand(d,1);
        Z = triu(r/4*randn(d));
        Z = Z- diag(diag(Z)) + diag(v);
        X = double(Q * Z * ctranspose(Q));
        II = eye(d);
        row = zeros(n+1,1);
        a = zeros(m+n+1,1);
        a(1) = 1;
        
        coef = alpha;
        
        for i = 2:m+n+1
            a(i) = double(coef); %This is binomial series (1+x)^alpha
            coef = double(coef * (alpha-i+1)/i);
        end
        actual = double((II+X)^alpha);

        %Cesaro sum
        dC = zeros(m+n+1,1);
        for i=1:m+n+1
            dC(i) = a(i) * (m+n+1-i) /(m+n+1);
        end
        row = zeros(n+1,1);
        row(1) = dC(1);
        
        R = toeplitz(dC(m+1:m+n),dC(m+1:-1:m-n+2));
        B = dC(m+2:m+n+1);
        if rank(R) ~= rank([R -B])
            disp("Cesaro Pade is ill-conditioned");
            break
        end
        gamma = R\(-B);
        gamma = [1, gamma']';
        beta = toeplitz(dC(1:m+1),row) * gamma;
        
        P = double(beta(1)*eye(d) + beta(2)*X); 
        Q = double(gamma(1)*eye(d) + gamma(2)*X);
        C1 = double(zeros(d));
        C2 = double(zeros(d));
        power = double(X);
        for i = 2:max([m n])
            power = double(power*X);
            if i <= m+1
                Y1 = double(beta(i+1)*power - C1);
                T1 = double(P + Y1);
                C1 = double(double(T1-P)-Y1);
                P = T1;
            end    
            if i <= n+1
                Y2 = double(gamma(i+1)*power - C2);
                T2 = double(Q + Y2);
                C2 = double(double(T2-Q)-Y2);
                Q = T2;
            end   
        end
        cesaro = double((Q'\P')');
        errorC(j) = norm(reshape(cesaro-actual,[],1),Inf);

        %Euler sum
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
                Y1 = double(beta(i+1)*power - C1);
                T1 = double(P + Y1);
                C1 = double(double(T1-P)-Y1);
                P = T1;
            end    
            if i <= n+1
                Y2 = double(gamma(i+1)*power - C2);
                T2 = double(Q + Y2);
                C2 = double(double(T2-Q)-Y2);
                Q = T2;
            end   
        end
        euler = double((Q'\P')');
        error2(j) = norm(reshape(euler-actual,[],1),Inf);
    end
    eulererror(k) = mean(error2);
    Cerror(k) = mean(errorC);
end
disp(eulererror);

fig = figure(1);
hold on
plot(log10(eulererror),'LineWidth',1)
plot(log10(Cerror),'LineWidth',1)
leg1 = legend("Pad\'e--Euler", "Pad\'e--Ces\`aro",'Fontsize',16,'Interpreter','latex','location','southeast');
set(leg1,'Interpreter','latex');
ylabel('log forward error','Fontsize',16,'Interpreter','latex')
xlabel('$n$','Fontsize',16,'Interpreter','latex')
hold off

