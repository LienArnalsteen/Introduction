%% numerical solution slide 8

rmax=0.48;
dens=20;
area=30528;
K=dens*area;
N(1)=50000;
teind=12;
deltat=1;
tijd=deltat:deltat:teind;


for t=2:numel(tijd);
    N(t)=rmax*((K-N(t-1))/K)*deltat*N(t-1)+N(t-1);

    
    
end
plot(tijd,N)
hold on
%% Analytical solution slide 9
pop(1)=N(1)
pop(1,2:numel(tijd))= ((-K*N(1)*exp(rmax*tijd(1,1:numel(tijd)-1)))/(N(1)-K))/(1-(N(1)*exp(rmax*tijd(1,1:numel(tijd)-1)))/(N(1)-K))
plot(tijd,pop)


%% hunting: back to logaritmic 
clearvars
clc
close all

rmax=0.48;
dens=20;
area=30528;
K=dens*area;
N(1)=50000;
teind=100;
deltat=1;
tijd=deltat:deltat:teind;
alfa=rmax+rmax*0.2;  %stable population if this equal r, +/-rmx*0.2 is if we hunt 20% more/less
H=alfa*N;


for t=2:numel(tijd);
    N(t)=((rmax*N(t-1)*deltat)-H)+N(t-1);
    H=alfa*N(t);
end
plot(tijd,N)

%% wolves
clearvars
clc
close all
r=0.48;
c=0.01;
d=0.24;
e=0.005;

N=100;
W=25;

alfa=0.1;  %stable population if this equal r, +/-rmx*0.2 is if we hunt 20% more/less
H=alfa*N

teind=100;
deltat=1;
tijd=deltat:deltat:teind;


for t=2:numel(tijd);
    Never=(r*N(t-1)*deltat+N(t-1)-c*W(t-1)*N(t-1)*deltat);%met hunting erbij: -alfa*N(t-1)
    Wolf=e*W(t-1)*Never*deltat+W(t-1)-d*W(t-1)*deltat;
    
    N = [N,Never];
    W = [W,Wolf];
end
plot(tijd,N)
hold on
plot(tijd,W)


%% les 3
%rmse bepalen
N0=100;
W0=25;
% r = 0.48;
% c = 0.01;
% d = 0.24;
% e = 0.005;

H_True=LotkaVolterra(N0,W0,r,c,d,e);

H_Mod=LotkaVolterra(N,W,r,c,d,e);

error=rmse(H_Mod,H_True);

error_procent=(error/H_True)*100;

figure
plot(Nboars)
hold on
plot(Nwolves)

figure
plot(H_Mod)





%% ode with roe deer
clearvars
clc
close all

N= [100,80,25];%boars,roe deer, wolf
c= 0.01;
e=0.005;
d=0.24;
r=0.48;
dt=0.1;

alfa=0.1;

T=[0 100]; %tijdsperiode
%%
options=odeset('Refine',20,'NormControl','on') %refine: meer stappen gebruiken
F= @(t,y)[y(1)*r-c*y(3)*y(1);y(2)*r-c*y(3)*y(2);y(3)*(-d)+e*y(3)*y(1)+e*y(3)*y(2)];
[t,r,w]=ode45(F,T,N,options)




B=r(:,1)
R=r(:,2)
W=r(:,3)


plot3(W,B,R)

colormap winter

%% hunting:
F= @(t,y)[y(1)*r-c*y(3)*y(1)-alfa;y(2)*r-c*y(3)*y(2);y(3)*(-d)+e*y(3)*y(1)+e*y(3)*y(2)]
[t,r,w]=ode45(F,T,N)






