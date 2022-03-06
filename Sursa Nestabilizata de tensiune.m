clear all; clc
U1=380;
I21n=6;
I22n=10;
I23n=10;
U21n=8;
U22n=14;
U23n=14;

P2=U21n*I21n+2*U22n*I22n
S2=P2;


%P2[W]=500; ntr=0.94; J=1
ntr=0.94; J=1;

Pg=(P2/2)*(1+(1/ntr))
Sfe=1.2*sqrt(Pg*((1+ntr)/(J*ntr)))*100 %mm^2

l=sqrt(Sfe) %mm
lstelat=60; %mm
%am ales tola E+60
b=(Sfe/lstelat)

bstelat=1.2*b %cm
Ntole=floor(b/0.35)

Sfestelat=lstelat*bstelat/100
lmed=(4*lstelat+2*bstelat)/10

Scu21=I21n/J
Scu22=I22n/J
Scu23=Scu22

d21=sqrt((4*Scu21)/pi)
d22=sqrt((4*Scu22)/pi)
d23=d22

I1=Pg/U1
ScuI=I1/J
d1=sqrt((4*ScuI)/pi)

f=50; %Hz
B=1.1;
e=4.44*f*B*(0.9*Sfestelat)*10^(-4) %[V/spira]

W1=floor(U1/e)
W21=floor(U21n/e)
W22=floor(U22n/e)
W23=W22
rezistivitatea=0.0172 %ohm*mm^2/m
R1=W1*((rezistivitatea*lmed)/(100*ScuI))
R21=W21*((rezistivitatea*lmed)/(100*Scu21))
R22=W22*((rezistivitatea*lmed)/(100*Scu22))
R23=R22

R21stelat=R21+R1*(W21/W1)^2
R22stelat=R22+R1*(W22/W1)^2
R23stelat=R22stelat

RS1n=U21n/I21n
RS2n=U22n/I22n
RS3n=RS2n

e21tilda=(W21/W1)*U1
e22tilda=(W22/W1)*U1
e23tilda=e22tilda

U21=e21tilda-R21stelat*I21n
U22=e22tilda-R22stelat*I22n
U23=U22

W21stelat=floor(W21 + ((R21stelat*I21n)/e))
W22stelat=floor(W22 + ((R22stelat*I22n)/e))
W23stelat=W22stelat
%%

%CALCULUL EXPEDITIV AL TRANFORMATORULUI DE MICA PUTERE FOLOSIND MONOGRAME

clear all; clc
U1=380;
I21n=6;
I22n=10;
I23n=10;
U21n=8;
U22n=14;
U23n=14;
ntr=0.94; J=1;
P2=U21n*I21n+2*U22n*I22n
P=(1/ntr)*P2
Pg=(P2/2)*(1+(1/ntr))
Sfe=12*sqrt(Pg*((1+ntr)/(J*ntr)))
Sm=Sfe


%% A2. Calculul circuitelor de redresare
clear all; clc
U1=380;
I21n=6;
I22n=10;
I23n=10;
U21n=8;
U22n=14;
U23n=14;
sigma=1.5; %coeficient de siguranta

Umax21 = sqrt(2)*U21n
Umax22 = sqrt(2)*U22n
Umax23 = Umax22

Umax21_stelat = 1.5*sqrt(2)*U21n
Umax22_stelat = 1.5*2*sqrt(2)*U22n
Umax23_stelat = Umax22_stelat

I21_barat = 0.7*I21n
I22_barat = 0.7*I22n
I23_barat = 0.7*I23n

I21_barat_stelat = 1.05*I21n
I22_barat_stelat = 1.05*I22n
I23_barat_stelat = 1.05*I23n

%%
clear variables; clc

U1=380;
I21n=6;
I22n=10;
I23n=10;
U21n=8;
U22n=14;
U23n=14;
Sm=28
N0i=1.2
N0ii=1.25
w1=N0i*U1
w21=N0ii*U21n
w22=N0ii*U22n
w23=w22
Scu1=0.89;
Scu21=6;
Scu22=10

diametre = [0.07 0.10 0.12 0.15 0.18 0.2 0.22 0.25 0.28 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.5 0.90 0.95 1 1.2 1.4 1.5 1.6 1.8 2]
SCu1_prim = zeros(1,length(diametre));
SCu21_prim = zeros(1,length(diametre));
SCu22_prim = zeros(1,length(diametre));

for i = 1:length(diametre)
    SCu1_prim(i)=((diametre(i)^2)*pi)/4
    SCu21_prim(i)=((diametre(i)^2)*pi)/4
    SCu22_prim(i)=((diametre(i)^2)*pi)/4
end

%%


%%
clc
clear all 
U1 = 380
micro0=0.5
fu=15;
w=314;
u0=8;
e0=9.03;
ku=u0/e0
Rr=0.13;
R0=1.33;
kr=Rr/R0
U21=8; I21=6;
Rsn1=U21/I21
q=7.444;
Rrt1=1.13;
Rs1n=1.333
Cfi_Stelat=(1600*q*(Rrt1+Rs1n))/(Rrt1*Rs1n)

U1_barat = ku*U21*sqrt(2)
U1_barat_max = sqrt(2)*U21
V0_barat=0.5*U1_barat_max
I01_barat=V0_barat/Rsn1
delta_I_tildaL1 = 0.25*I01_barat
V_tildao1 =0.04*V0_barat
fc=1/(pi*w*2)

delta_Vo_tilda=4*V0_barat
Ua1=1.25*U1_barat_max
Ia1=1.5*I01_barat
U_tranzistor1=2*sqrt(2)*U21
I_tranzistor1=1.5*I01_barat
L1=((U1_barat_max-V0_barat)*micro0)/(fc*delta_I_tildaL1)
C=(delta_I_tildaL1*micro0)/(fc*delta_Vo_tilda)
Cales=15
R3=1/(2*fu*C)
A=8
micro=0.5

ui=micro*2*A-8
k=1/(2*A)
kr=0.2

V0_stelat=micro*kr*U1_barat
Co=19*10^(3)
 %k1=129.5756
ceva_pe_acolo = w*Co*R0
Rrt1=0.13
Cf=19*10^(-3)
k1=1/(Rrt1*Cf)
k2=1/(Rsn1*Cf)
%% C calculul regulatorului si simularea functionarii surselor controlate de curent
A=8;
C=46.24;
U21=8;
L=0.52;
Rsn1=1.33;
Hf=tf([(1/(2*A))*sqrt(2)*U21],[L*C L/(Rsn1) 1])
%%
A=8
U1_barat=8;
U21=8
U1max=sqrt(2)*U21
Vo_barat=0.5*U1_barat
I01=Vo_barat/Rs1n
deltaL=0.25*I01
Vo_tilda=4*Vo_barat
tildaL=0.25*I01
L1=((U1_barat_max-Vo_barat)*micro0)/(fc*tildaL)
C1=(tildaL*0.5)/(fu*Vo_tilda)
UA=1.25*U1max
IA=1.5*I01
Utranistor=2*sqrt(2)*8
Itranzistor=1.5*I01
%%
A=8;
U21=8
Hc=tf(407.6,[1 470.1 650.6])
Hr=0.2;
Hdes=Hf*Hr*Hc;
Hd=Hf*Hc;
H0=Hd/(1+Hdes)
step(H0);
bode(Hdes)
% hr=tf([4.93^2 1],[4.93 0])
% Hdes=series(hr,Hf)
% Ho=feedback(Hdes,1)
% figure,
% step(Ho); grid on
% figure,
% bode(Hdes); grid on
