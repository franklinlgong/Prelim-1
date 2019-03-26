%Unperturbed Concentration Averages 
phase1p1 = 0.0018; %mean(X([60:80],4))
phase2p1 = 2.6145e+03; %mean(X([120:140],4))
phase3p1 = 1.9537e+04; %mean(X([400:420],4))
phase1p2 = 6.3135;
phase2p2 = 1.7808e+03;
phase3p2 = 1.3239e+04;
phase1p3 = 0.5429;
phase2p3 = 620.4449;
phase3p3 = 4.8475e+03;

t_span = 0:1:420; 
[m,n] = size(t_span); 
I = zeros(n+1,1);
I(120:420,1) = 10;
x0 = [0;0;0;0;0;0];

[t,X] = ode45(@(t,x) prelim1problem2system(t,x,I),t_span,x0);

Unperturbed = [phase1p1, phase2p1, phase3p1;
                phase1p2, phase2p2, phase3p2;
                phase1p3, phase2p3, phase3p3;];

Array1 = [mean(X([60:80],4)), mean(X([120:140],4)), mean(X([400:420],4));
    mean(X([60:80],5)), mean(X([120:140],5)), mean(X([400:420],5));
    mean(X([60:80],6)), mean(X([120:140],6)), mean(X([400:420],6));]