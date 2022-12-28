%Static Analysis of a Euler Beam
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
I=[4e6]*1e-12;
P=0;            %Concentrated Load
qo=[0 12]*1e3;               %uniform load
Cood=[0 0;1 0;2 0];  %coordinates of nodes
%% Meshing
ne=2;
nne=2;
nn=ne*nne-ne+1;       
dof=2;
dofe=dof*nne;
tdof=nn*dof;
CONN=zeros(ne,dofe);
NCONN=zeros(ne,nne);
theta=zeros(ne,1);
CONN=[1 2 3 4;3 4 5 6];
NCONN=[1 2;2 3];
for i=1:ne
    le(i)=sqrt((Cood(NCONN(i,2),2)-Cood(NCONN(i,1),2))^2 + (Cood(NCONN(i,2),1)-Cood(NCONN(i,1),1))^2);
    theta(i)=atan2d(Cood(NCONN(i,2),2)-Cood(NCONN(i,1),2), Cood(NCONN(i,2),1)-Cood(NCONN(i,1),1));
end

le
theta


%% Assembly of Element Matrices
KG=zeros(tdof,tdof);
FG=zeros(tdof,1);
KeG=zeros(dofe,dofe);
FeG=zeros(dofe,1);
Te=zeros(nne,dofe);

for i=1:ne
    KeG=E*I/le(i)^3*[12 6*le(i) -12 6*le(i);6*le(i) 4*le(i)^2 -6*le(i) 2*le(i)^2;-12 -6*le(i) 12 -6*le(i);6*le(i) 2*le(i)^2 -6*le(i) 4*le(i)^2];       
    FeG=qo(i)*le(i)/2*[1;le(i)/6;1;le(i)^2/6];
    %Te=[cosd(theta(i)) sind(theta(i)) 0 0; 0 0 cosd(theta(i)) sind(theta(i))];
    %KeG=transpose(Te)*KeL*Te;
    %FeG=transpose(Te)*FeL;
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+KeG(j,k);
        end
        FG(CONN(i,j),1)=FG(CONN(i,j),1)+FeG(j,1);
    end
end
KG
FG
FGC=zeros(tdof,1);

FG=FG+FGC;

KGr=KG
FGr=FG;

%% Applications of BCs( u1,v1,v3=0)

for kk=[5,3,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    FG(kk,:)=[];
end
UG=inv(KG)*FG;
UGr=zeros(tdof,1);
UGr(4)=UG(1);
UGr(6)=UG(2);
UGr





