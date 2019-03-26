%Initial conditions
t_span = 0:1:420; 
[m,n] = size(t_span); 
I = zeros(n+1,1);
I(120:420,1) = 10;
x0 = [0;0;0;0;0;0];

[t,X] = ode45(@(t,x) prelim1problem2system(t,x,I),t_span,x0);

figure(1)
q = plot(t_span,X(:,4),t_span,X(:,5),t_span,X(:,6));
xlabel('Time (min)')
ylabel('Protein concentration (umol/gDW)')
legend('Protein 1','Protein 2','Protein 3')
title('Incoherent Feed Forward Loop')