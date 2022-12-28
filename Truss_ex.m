%Static Analysis of a truss planar
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
A=[30e-4];          %Areas
P=[100e3 -200e3];            %Concentrated Load
qo=0;               %uniform load
Cood=[0 0;100 0;100 100]*1e-2;  %coordinates of nodes
%% Meshing
ne=3;
nne=2;
nn=ne*nne-ne;       %for truss different formulae
dof=2;
dofe=dof*nne;
tdof=nn*dof;
CONN=zeros(ne,dofe);
NCONN=zeros(ne,dof);
theta=zeros(ne,1);
k=1;
for i=1:ne
    for j=1:dofe
        if k>tdof
            k=1;
        end
        CONN(i,j)=k;
        k=k+1;
    end
    k=k-2;
          
end
k=1;
for i=1:ne
    for j=1:nne
        if k>nn
            k=1;
        end
        NCONN(i,j)=k;
        k=k+1;
    end
    k=k-1;   
end

for i=1:ne
    le(i)=sqrt((Cood(NCONN(i,2),2)-Cood(NCONN(i,1),2))^2 + (Cood(NCONN(i,2),1)-Cood(NCONN(i,1),1))^2);
    theta(i)=atan2d(Cood(NCONN(i,2),2)-Cood(NCONN(i,1),2), Cood(NCONN(i,2),1)-Cood(NCONN(i,1),1));
end   

%% Assembly of Element Matrices
KG=zeros(tdof,tdof);
FG=zeros(tdof,1);
KeG=zeros(dofe,dofe);
FeG=zeros(dofe,1);
Te=zeros(nne,dofe);

for i=1:ne
    KeL=E*A/le(i)*[1 -1; -1 1];       
    FeL=qo*le(i)/2*[1;1];
    Te=[cosd(theta(i)) sind(theta(i)) 0 0; 0 0 cosd(theta(i)) sind(theta(i))];
    KeG=transpose(Te)*KeL*Te;
    FeG=transpose(Te)*FeL;
    
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+KeG(j,k);
        end
        FG(CONN(i,j),1)=FG(CONN(i,j),1)+FeG(j,1);
    end
end
FGC=zeros(tdof,1);
FGC(5,1)=FGC(5,1)+P(1);
FGC(6,1)=FGC(6,1)+P(2);

FG=FG+FGC;

KGr=KG;
FGr=FG;

%% Applications of BCs( u1,v1,u2,v2=0)

for kk=[4,3,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    FG(kk,:)=[];
end
KG
FG
UG=inv(KG)*FG

UGr=zeros(tdof,1);
UGr(5)=UG(1);
UGr(6)=UG(2);

R=KGr*UGr-FGr;
R
