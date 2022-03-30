

function EvolutionDiscrete(V0,BC0,Kab,PCspf,PCgc,Rg,WeakMem)

time=[1:1601];
dt=1/((length(time)-1)/16);

fitness=zeros(11*11*11,1);

Viruslocinit= 0*ones(1,3);
for ii=1:11
    for jj=1:11
        for kk=1:11
            fitness((ii-1)*121+(jj-1)*11+kk)=max(0,1-sqrt((Viruslocinit(1)-(0.2*kk-1.2))^2+(Viruslocinit(2)-(0.2*jj-1.2))^2+(Viruslocinit(3)-(0.2*ii-1.2))^2)/sqrt(3));        end
    end
end
fitness=fitness'; GCfitness=fitness;

r=2;



% RsArray=[2:-.04:0]*r;
% RsArray=[.7:.05:1.3]*r;
% RsArray=[.7 1 1.3]*r;
RsArray=[.9: .1 :1.2]*r;
% RsArray=r;
RgArray=[1.4:.1:3.5];
% RgArray=[.7 1 1.3]*r;
% RgArray=r;
% RsArray=[2,4];
% RgArray=[2,4];
% RsArray=[.7,1.3]*r;
% RgArray=[0.3:0.3:1.5];
% RgArray=[1.48:-.04:0]*r;
% RgArray=[2,1]*r;
% RgArray=[2:-.04:0]*r;
% Costarray=zeros(length(RsArray),length(RsArray));
% VirusCell=cell(length(RsArray),1);
Virus2Cell=cell(length(RsArray),length(RgArray),6);
DistCell=cell(length(RgArray));
GC2DistCell=cell(length(RsArray),length(RgArray),6);
SPF2DistCell=cell(length(RsArray),length(RgArray),6);
k=0

for Rg=RgArray
    k=k+1;
fitness=zeros(11*11*11,1);
Virusloc= Viruslocinit;
for ii=1:11
    for jj=1:11
        for kk=1:11
            fitness((ii-1)*121+(jj-1)*11+kk)=max(0,1-sqrt((Virusloc(1)-(0.2*kk-1.2))^2+(Virusloc(2)-(0.2*jj-1.2))^2+(Virusloc(3)-(0.2*ii-1.2))^2)/sqrt(3));        end
    end
end
fitness=fitness'; GCfitness=fitness;

    
    
Virus=zeros(length(time),1); BCs=zeros(length(time),length(fitness)); GCBCs=zeros(length(time),length(GCfitness));
BCs(1,:)=1/length(fitness); GCBCs(1,:)=1/length(fitness);
BCs(1,:)=BC0*BCs(1,:)*10^-6; GCBCs(1,:)=BC0*GCBCs(1,:); %need to multiply by a small number to not break it
Abs=zeros(size(BCs));
Virus(1)=V0;
Rs=0;
for t=time
if Kab*Virus(t)*dt^2>= .1

dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,.001*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
for i=1:999
dy=virusSim(BCs(t+1,:),GCBCs(t+1,:),Abs(t+1,:),Virus(t+1),fitness,Rs,.001*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t+1,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t+1,:)+dy(2,:); Abs(t+1,:)=Abs(t+1,:)+dy(3,:); Virus(t+1)=Virus(t+1)+dy(4,1);

end
    elseif Kab*Virus(t)*dt^2 >= 1*10^-3

dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,.01*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
for i=1:99
dy=virusSim(BCs(t+1,:),GCBCs(t+1,:),Abs(t+1,:),Virus(t+1),fitness,Rs,.01*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t+1,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t+1,:)+dy(2,:); Abs(t+1,:)=Abs(t+1,:)+dy(3,:); Virus(t+1)=Virus(t+1)+dy(4,1);

end
elseif Kab*Virus(t)*dt^2 >= 1*10^-5

    dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,.1*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
for i=1:9
dy=virusSim(BCs(t+1,:),GCBCs(t+1,:),Abs(t+1,:),Virus(t+1),fitness,Rs,.1*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t+1,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t+1,:)+dy(2,:); Abs(t+1,:)=Abs(t+1,:)+dy(3,:); Virus(t+1)=Virus(t+1)+dy(4,1);

end
else
dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
end

if Virus(t)<1
    Virus(t)=0;
    Virus(t+1)=0;
%     if Virus(t-300)==0
%     break
%     end
elseif Virus(t) >=10^9
    break
elseif isnan(Virus(t))
    Virus(t)=1.1*10^9;
    break
end


end
% BCDist(1,:)=BCs(t,:)/sum(BCs(t,:));
%  figure
%  plot(log(Virus)/log(10))
% keyboard;

if WeakMem==1
     GCBCs(301,GCfitness>=0.5)=0; GCBCDist(1,:)=GCBCs(301,:)/sum(GCBCs(301,:));
elseif WeakMem==2
    
GCBCs(301,GCfitness>=0.5)=GCBCs(301,GCfitness>=0.5).*(1-GCfitness(GCfitness>=0.5));
GCBCDist(1,:)=GCBCs(301,:)/sum(GCBCs(301,:));
elseif WeakMem==3
GCBCDist(1,:)=GCBCs(301,:)/sum(GCBCs(301,:));
else
    GCBCDist(1,:)=GCBCs(t,:)/sum(GCBCs(t,:));
end


% VirusCell=Virus;
DistCell{k}=GCBCDist;


%hit Virus 2
%    graphmaker(BCDist,GCBCDist,fitness,Rg)
m=0;
for Rs=RsArray
    m=m+1;
j=0;
% for Vinit=[0:0.2:1];
for Vinit=[0, 0.8];
    j=j+1;
Virus=zeros(length(time),1); BCs=zeros(length(time),length(fitness)); GCBCs=zeros(length(time),length(GCfitness));
BCs(1,:)=BC0*GCBCDist; GCBCs(1,:)=BC0*GCBCDist;
Abs=zeros(size(BCs));
Virus(1)=V0*100;

GCBCs(1,GCBCs(1,:)<1)=0; BCs(1,BCs(1,:)<1)=0; %force discretization

fitness=zeros(11*11*11,1);
Virusloc= Vinit*[1.0 1.0 1.0];
for ii=1:11
    for jj=1:11
        for kk=1:11
            fitness((ii-1)*121+(jj-1)*11+kk)=max(0,1-sqrt((Virusloc(1)-(0.2*kk-1.2))^2+(Virusloc(2)-(0.2*jj-1.2))^2+(Virusloc(3)-(0.2*ii-1.2))^2)/sqrt(3));
        end
    end
end

fitness=fitness'; GCfitness=fitness;
GCMeanFit=sum(GCBCDist.*fitness);
% SPFMeanFit=sum(BCDist.*fitness);
for t=time
%  maxfit(t)=max(GCBCs(t,GCBCs(t,:)<1));

%     GCBCs(t,GCBCs(t,:)<1)=0; %force discretization at every timestep?  Could maybe do better?

    if Kab*Virus(t)*dt^2>= .1

dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,.001*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
for i=1:999
dy=virusSim(BCs(t+1,:),GCBCs(t+1,:),Abs(t+1,:),Virus(t+1),fitness,Rs,.001*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t+1,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t+1,:)+dy(2,:); Abs(t+1,:)=Abs(t+1,:)+dy(3,:); Virus(t+1)=Virus(t+1)+dy(4,1);

end
    elseif Kab*Virus(t)*dt^2 >= 1*10^-3

dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,.01*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
for i=1:99
dy=virusSim(BCs(t+1,:),GCBCs(t+1,:),Abs(t+1,:),Virus(t+1),fitness,Rs,.01*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t+1,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t+1,:)+dy(2,:); Abs(t+1,:)=Abs(t+1,:)+dy(3,:); Virus(t+1)=Virus(t+1)+dy(4,1);

end
elseif Kab*Virus(t)*dt^2 >= 1*10^-5

    dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,.1*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
for i=1:9
dy=virusSim(BCs(t+1,:),GCBCs(t+1,:),Abs(t+1,:),Virus(t+1),fitness,Rs,.1*dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t+1,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t+1,:)+dy(2,:); Abs(t+1,:)=Abs(t+1,:)+dy(3,:); Virus(t+1)=Virus(t+1)+dy(4,1);

end
else
dy=virusSim(BCs(t,:),GCBCs(t,:),Abs(t,:),Virus(t),fitness,Rs,dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf);
BCs(t+1,:)=BCs(t,:)+dy(1,:); GCBCs(t+1,:)=GCBCs(t,:)+dy(2,:); Abs(t+1,:)=Abs(t,:)+dy(3,:); Virus(t+1)=Virus(t)+dy(4,1);
end

if Virus(t)<1
    Virus(t)=0;
    Virus(t+1)=0;
    break
elseif Virus(t) >=10^9
    break
elseif isnan(Virus(t))
    Virus(t)=1.1*10^9;
    break
end


end


Virus2Cell{m,k,j}=Virus;

if WeakMem==1
     GCBCs(301,GCfitness>=0.5)=0; GCBCDist2(1,:)=GCBCs(301,:)/sum(GCBCs(301,:));
elseif WeakMem==2
    GCBCs(301,GCfitness>=0.5)=GCBCs(301,GCfitness>=0.5).*(1-GCfitness(GCfitness>=0.5));
GCBCDist2(1,:)=GCBCs(301,:)/sum(GCBCs(301,:));
elseif WeakMem==3
GCBCDist2(1,:)=GCBCs(301,:)/sum(GCBCs(301,:));
else
    GCBCDist2(1,:)=GCBCs(t,:)/sum(GCBCs(t,:));
end
BCDist2(1,:)=BCs(t,:)/sum(BCs(t,:));

GC2DistCell{m,k,j}=GCBCDist2;
SPF2DistCell{m,k,j}=BCDist2;
end
end
end

Flagstr=['FunctionEv_V0',num2str(V0),'_BC0',num2str(BC0),'_Kab',num2str(Kab),'_Ps',num2str(PCspf),'_Pg',num2str(PCgc),'mu=0.4'];

if WeakMem==1
save([Flagstr,'WeakMemory.mat'],'RsArray','RgArray','Virus2Cell','DistCell','GC2DistCell','SPF2DistCell');
elseif WeakMem==2
 save([Flagstr,'WeakMemoryLinear.mat'],'RsArray','RgArray','Virus2Cell','DistCell','GC2DistCell','SPF2DistCell');
 elseif WeakMem==3
  save([Flagstr,'TimeMem.mat'],'RsArray','RgArray','Virus2Cell','DistCell','GC2DistCell','SPF2DistCell');
   
else
save([Flagstr,'.mat'],'RsArray','RgArray','Virus2Cell','DistCell','GC2DistCell','SPF2DistCell');

end

end

function [y]=virusSim(BCs,GCBCs,Abs,Virus,fitness,Rs,dt,t,GCfitness,Rg,Kab,r,PCgc,PCspf)
mu=.4  ; 

k1BC=BCs.*fitness*Rs/(sum(BCs.*fitness)/sum(BCs));
k2BC=(BCs+k1BC*dt*.5).*fitness*Rs/(sum((BCs+k1BC*dt*.5).*fitness)/sum((BCs+k1BC*dt*.5)));
k3BC=(BCs+k2BC*dt*.5).*fitness*Rs/(sum((BCs+k2BC*dt*.5).*fitness)/sum((BCs+k2BC*dt*.5)));
k4BC=(BCs+k3BC*dt).*fitness*Rs/(sum((BCs+k3BC*dt).*fitness)/sum((BCs+k3BC*dt)));
y=dt/6*(k1BC+2*k2BC+2*k3BC+k4BC);

GCsumfitness=sum(GCBCs.*GCfitness);
    GCmeanfitness=GCsumfitness/sum(GCBCs);
% GCBCs(t+1,:)=GCBCs+(0.5*GCBCs.*GCfitness/GCmeanfitness+0.5/(length(GCfitness)-1)*(GCsumfitness-GCBCs.*GCfitness)/GCsumfitness)*dt*Rg;
k1GCBC=((1-mu)*GCBCs.*GCfitness/GCmeanfitness+mu/(length(GCfitness)-1)*(GCsumfitness-GCBCs.*GCfitness)/GCmeanfitness)*Rg;
 k2GCsumfitness=sum((GCBCs+dt*0.5*k1GCBC).*GCfitness); k2GCmeanfitness=k2GCsumfitness/sum(GCBCs+dt*0.5*k1GCBC);
 k2GCBC=((1-mu)*(GCBCs+dt*0.5*k1GCBC).*GCfitness/k2GCmeanfitness+mu/(length(GCfitness)-1)*(k2GCsumfitness-(GCBCs+dt*0.5*k1GCBC).*GCfitness)/k2GCmeanfitness)*Rg;
k3GCsumfitness=sum((GCBCs+dt*0.5*k2GCBC).*GCfitness); k3GCmeanfitness=k3GCsumfitness/sum(GCBCs+dt*0.5*k2GCBC);
 k3GCBC=((1-mu)*(GCBCs+dt*0.5*k2GCBC).*GCfitness/k3GCmeanfitness+mu/(length(GCfitness)-1)*(k3GCsumfitness-(GCBCs+dt*0.5*k2GCBC).*GCfitness)/k3GCmeanfitness)*Rg;
k4GCsumfitness=sum((GCBCs+dt*k3GCBC).*GCfitness); k4GCmeanfitness=k4GCsumfitness/sum(GCBCs+dt*k3GCBC);
 k4GCBC=((1-mu)*(GCBCs+dt*k3GCBC).*GCfitness/k4GCmeanfitness+mu/(length(GCfitness)-1)*(k4GCsumfitness-(GCBCs+dt*k3GCBC).*GCfitness)/k4GCmeanfitness)*Rg;
y=[y;dt/6*(k1GCBC+2*k2GCBC+2*k3GCBC+k4GCBC)];

    k1Abs=PCspf*BCs+PCgc*GCBCs-Kab*Abs.*fitness*Virus;
    k1Vir=Virus*(r-Kab*sum(Abs.*fitness));

    k2Abs=PCspf*(BCs+0.5*dt*k1BC)+PCgc*(GCBCs+0.5*dt*k1GCBC)-Kab*(Abs+dt*0.5*k1Abs).*fitness*(Virus+dt*0.5*k1Vir);
    k2Vir=(Virus+0.5*dt*k1Vir)*(r-Kab*sum((Abs+dt*0.5*k1Abs).*fitness));
    k3Abs=PCspf*(BCs+0.5*dt*k2BC)+PCgc*(GCBCs+0.5*dt*k2GCBC)-Kab*(Abs+dt*0.5*k2Abs).*fitness*(Virus+dt*0.5*k2Vir);
    k3Vir=(Virus+0.5*dt*k2Vir)*(r-Kab*sum((Abs+dt*0.5*k2Abs).*fitness));
    k4Abs=PCspf*(BCs+dt*k3BC)+PCgc*(GCBCs+dt*k3GCBC)-Kab*(Abs+dt*k3Abs).*fitness*(Virus+dt*k3Vir);
    k4Vir=(Virus+dt*k3Vir)*(r-Kab*sum((Abs+dt*k3Abs).*fitness));
   y=[y;dt/6*(k1Abs+2*k2Abs+2*k3Abs+k4Abs)];
   y=[y; zeros(1,length(y))];

   y(4,1)=dt/6*(k1Vir+2*k2Vir+2*k3Vir+k4Vir);
end

