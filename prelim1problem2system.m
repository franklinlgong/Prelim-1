function dxdt = prelim1problem2system(t,x,I)

%Define parameters
m1 = x(1);
m2 = x(2);
m3 = x(3);
p1 = x(4);
p2 = x(5);
p3 = x(6);

Avog = 6.022*10^23*10^-9; %copies/nmol Avogadro's for converting copy number into nmol

%Biophysical Constants
KX = 0.24;	%nmol/gDW Transcription Saturation Constant
KL = 454.64;	%nmol/gDW Translation Saturation Constant
Kix = 1/(42/60); %1/min Transcription initiation rate
Kil = 1/(15/60); %1/min Translation initiation rate
Dt = 40;    %min Doubling time
Cm = 2.8*10^-13;    %g/cell Single cell mass
Dw = 0.3;   %gDw/g Percent dry weight of cell
Rxt = (1150/Avog)/(Cm*Dw);  %nmol/gDW Concentration of RNAP
Rib = (45000/Avog)/(Cm*Dw);    %nmol/gDW Concentration of ribosomes
Kep = 60*60;   %nt/min Transcription elongation rate
Klp = 16.5*60; %aa/min Translation elongation rate
Lx = 1000;  %nt Characteristic transcript length
Lt = 333;   %aa Characteristic protein length
mRNA_hl = 2.1;	%mRNA half-life (min)
prot_hl = 24*60;	%Protein half-life (min)
Gc = (200/Avog)/(Cm*Dw); 	%nmol/gDW Copy number of plasmid
Dx = log(2)/mRNA_hl; 	%Degradation rate of mRNA
Dl = log(2)/prot_hl;	%Degradation rate of proteins
B = log(2)/Dt;	%Dilution factor

%Species Parameters
Lx1 = 1200; %Gene 1 Length (nts)
Lx2 = 2400; %Gene 2 Length (nts)
Lx3 = 600;  %Gene 3 Length (nts)
Ll1 = 400; 	%Protein 1 Length (AA)
Ll2 = 800;  %Protein 2 Length (AA)
Ll3 = 200;  %Protein 3 Length (AA)

Ke1 = Kep/Lx1;  %Transcription Elongation constant rate for mRNA1
Ke2 = Kep/Lx2;  %Transcription Elongation constant rate for mRNA2
Ke3 = Kep/Lx3;  %Transcription Elongation constant rate for mRNA3
Kl1 = Klp/Ll1; 	%Translation Elongation rate protein 1
Kl2 = Klp/Ll2; 	%Translation Elongation rate protein 2
Kl3 = Klp/Ll3;	%Translation Elongation rate protein 3

tx1 = Ke1/Kix;	%Tau for gene 1
tx2 = Ke2/Kix; 	%Tau for gene 2
tx3 = Ke3/Kix;  %Tau for gene 3
tl1 = Kl1/Kil; 	%Tau for protein 1
tl2 = Kl2/Kil; 	%Tau for protein 2
tl3 = Kl3/Kil; 	%Tau for protein 3

rx1 = Ke1*Rxt*(Gc/(KX*tx1+Gc*tx1+Gc)); %Transcription rate for gene 1
rx2 = Ke2*Rxt*(Gc/(KX*tx2+Gc*tx2+Gc)); %Transcription rate for gene 2
rx3 = Ke3*Rxt*(Gc/(KX*tx3+Gc*tx3+Gc)); %Transcription rate for gene 3

n = 1.5;
a = round(t);
In = I(a+1,1);

%Binding functions
fi1 = (In^n)/(In^n+0.3^n); 
f12 = (p1^n)/(p1^n+1); 
f13 = (p1^n)/(p1^n+1);
f23 = (p2^10)/(p2^10+1^10); 

%Weights
wi1 = 100;  %Effect of inducer on P1 expression
w11 = 0.0000001; %Background expression of P1
w12 = 10;   %Effect of P1 on P2 expression
w13 = 5;    %Effect of P1 on P3 expression
w22 = 0.0000001; %Background expression of P2
w23 = 25;   %Effect of P2 on P3 expression
w33 = 0.0000001; %Background expression of P3

%Transcription rules
um1 = (w11+(wi1*fi1))/(1+w11+(wi1*fi1)); 
um2 = (w22+(w12*f12))/(1+w22+(w12*f12));
um3 = (w33+(w13*f13))/(1+w33+(w13*f13)+(w23*f23)); 

%Kinetic translation rates
rl1 = Kl1*Rib*(m1/(KL*tl1+m1*tl1+m1)); 
rl2 = Kl2*Rib*(m2/(KL*tl2+m2*tl2+m2)); 
rl3 = Kl3*Rib*(m3/(KL*tl3+m3*tl3+m3)); 

%Setup Differential Equations Matrices
vec=[m1;m2;m3;p1;p2;p3];
r=[rx1*um1;rx2*um2;rx3*um3;rl1;rl2;rl3];
A=[-(B+Dx),0,0,0,0,0;
   0,-(B+Dx),0,0,0,0;
   0,0,-(B+Dx),0,0,0;
   0,0,0,-(B+Dl),0,0;
   0,0,0,0,-(B+Dl),0;
   0,0,0,0,0,-(B+Dl)];
S=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0;
   0,0,0,1,0,0;
   0,0,0,0,1,0;
   0,0,0,0,0,1];

%Create differential equation:
dxdt = (A*vec)+(S*r);


end