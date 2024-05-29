close all
clear
clc

Freq = 100;
X = linspace(-1.5*pi,1.5*pi,6000);
sums = 0;
Sigma = 0;

%simply taking the truncated Fourier series
for k = 1:Freq
    
    sums = sums + 2/(pi*k)*(1-(-1)^k)*sin(k*X);
    Sigma = Sigma +sums;

end

hold on
figure(1)
plot(X,sums,'LineWidth',1)
plot(X,Sigma/(Freq+1),'LineWidth',1)
leg1 = legend('Fourier',"Ces\'{a}ro",'Fontsize',16,'Interpreter','latex','Location','best');
ylabel('$y$','Fontsize',16,'Interpreter','latex')
ylim([-2 2])
xlabel('$x$','Fontsize',16,'Interpreter','latex')
set(leg1,'Interpreter','latex');

figure(2)
hold on
plot(X,abs(sums),'LineWidth',1)
plot(X,abs(Sigma/(Freq+1)),'LineWidth',1)
leg2 = legend('Fourier',"Ces\'{a}ro",'Fontsize',16,'Interpreter','latex');
set(leg2,'Interpreter','latex');
ylabel('absolute values','Fontsize',16,'Interpreter','latex')
xlabel('$x$','Fontsize',16,'Interpreter','latex')
hold off