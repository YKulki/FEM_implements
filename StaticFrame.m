%Static Analysis of a Static Frame
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
I=[20e6]*1e-12;
A=[4e4]*1e-6;
P=0;            %Concentrated Load
qoA=[0 0 0];
qoB=[0 10 0]*1e3;               %uniform load
Cood=[0 0;0 40;25 40;25 0]*1e-2;  %coordinates of nodes
%% Meshing
ne=3;
nne=2;
nn=ne*nne-ne+1;       
dof=3;
dofe=dof*nne;
tdof=nn*dof;
CONN=zeros(ne,dofe);
NCONN=zeros(ne,nne);
theta=zeros(ne,1);
CONN=[1 2 3 4 5 6;4 5 6 7 8 9;7 8 9 10 11 12];
NCONN=[1 2;2 3;3 4];
for i=1:ne
    le(i)=sqrt((Cood(NCONN(i,2),2)-Cood(NCONN(i,1),2))^2 + (Cood(NCONN(i,2),1)-Cood(NCONN(i,1),1))^2);
    theta(i)=atan2d(Cood(NCONN(i,2),2)-Cood(NCONN(i,1),2), Cood(NCONN(i,2),1)-Cood(NCONN(i,1),1));
end



%% Assembly of Element Matrices
KG=zeros(tdof,tdof);
FG=zeros(tdof,1);
KeG=zeros(dofe,dofe);
FeG=zeros(dofe,1);
Te=zeros(dofe,dofe);

for i=1:ne
    m=cosd(theta(i));
    n=sind(theta(i));
    KeL=E*I/le(i)^3*[A*le(i)^2/I 0 0 -A*le(i)^2/I 0 0;0 12 6*le(i) 0 -12 6*le(i);0 6*le(i) 4*le(i)^2 0 -6*le(i) 2*le(i)^2;-A*le(i)^2/I 0 0 A*le(i)^2/I 0 0;0 -12 -6*le(i) 0 12 -6*le(i);0 6*le(i) 2*le(i)^2 0 -6*le(i) 4*le(i)^2];       
    FeL=le(i)/2*[qoA(i);qoB(i);qoB(i)*le(i)/6;qoA(i);qoB(i);-qoB(i)*le(i)/6];
    Te=[m n 0 0 0 0;-n m 0 0 0 0;0 0 1 0 0 0;0 0 0 m n 0;0 0 0 -n m 0;0 0 0 0 0 1];
    KeG=transpose(Te)*KeL*Te;
    FeG=transpose(Te)*FeL;
    
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

for kk=[12,11,10,3,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    FG(kk,:)=[];
end
UG=inv(KG)*FG
UGr=zeros(tdof,1);
UGr(4)=UG(1);
UGr(5)=UG(2);
UGr(6)=UG(3);
UGr(7)=UG(4);
UGr(8)=UG(5);
UGr(9)=UG(6);

UGr



