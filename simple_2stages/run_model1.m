close all
clear
clc

tRange = [0 200]; %t_start and t_end
Y0 = [0.15; 0.1; 0.1; 0.1]; %initial conditions for P, R, J, A

[tSol, YSol] = ode45(@model1, tRange, Y0);

P = YSol(:,1);
R = YSol(:,2);
J = YSol(:,3);
A = YSol(:,4);

figure
subplot(3,1,1)
plot(tSol, P, '-k')
ylabel('P concentration (mg P/L)') 
subplot(3,1,2)
plot(tSol, R, '-g')
ylabel('Algal density (mg C/L)')
subplot(3,1,3)
plot(tSol, J, '-b')
hold on
plot(tSol, A, '-r')
xlabel('Time (d)') 
ylabel('Grazer density (mg C/L)')
legend('Juvenile','Adult')
hold off