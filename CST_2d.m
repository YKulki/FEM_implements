%Static Analysis of a Static Frame
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
Cood=[0 0;1 0;1 1;0 1]; %coordinates of nodes
Ae=0.5;
He=0.02;
V=0.25;
P1=-1e6;
bx=[0 0];
by=[0 0];
Tx=[0 0;0 0];
Ty=[0 0;0 0];
LT=[1;1];

%% Meshing
ne=2;
nne=3;
nn=ne*nne-ne;       
dof=2;
dofe=dof*nne;
tdof=nn*dof;
CONN=zeros(ne,dofe);
NCONN=zeros(ne,nne);
CONN=[1 2 3 4 5 6;1 2 5 6 7 8];
NCONN=[1 2 3;1 3 4];


%% Assembly of Element Matrices
KG=zeros(tdof,tdof);
FG=zeros(tdof,1);
Ke=zeros(dofe,dofe);
Fe=zeros(dofe,1);
D=zeros(3,3);
Be=zeros(3,6);

for i=1:ne
    x12=Cood(NCONN(i,1),1)-Cood(NCONN(i,2),1);
    x13=Cood(NCONN(i,1),1)-Cood(NCONN(i,3),1);
    x23=Cood(NCONN(i,2),1)-Cood(NCONN(i,3),1);

    y12=Cood(NCONN(i,1),2)-Cood(NCONN(i,2),2);
    y13=Cood(NCONN(i,1),2)-Cood(NCONN(i,3),2);
    y23=Cood(NCONN(i,2),2)-Cood(NCONN(i,3),2);

    J = x12*y23 - x23*y13;

    D=E/(1-V^2)*[1 V 0;V 1 0; 0 0 (1-V)/2];
    Be=(1/J)*[y23 0 -y13 0 y12 0;0 -x23 0 x13 0 -x12;-x23 y23 x13 -y13 -x12 y12];
    Ke=transpose(Be)*D*Be*Ae*He;
    FeB=(Ae*He/3)*[bx(i);by(i);bx(i);by(i);bx(i);by(i)];
    FeTx=He*LT(i)*[(Tx(i,1)/3)+(Tx(i,2)/6);0;(Tx(i,1)/6)+(Tx(i,2)/3);0;0;0];
    FeTy=He*LT(i)*[0;(Ty(i,1)/3)+(Ty(i,2)/6);0;(Ty(i,1)/6)+(Ty(i,2)/3);0;0];
    Fe=FeB+FeTy+FeTy;
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
        end
        FG(CONN(i,j),1)=FG(CONN(i,j),1)+Fe(j,1);
    end
end

KG;
FGC=zeros(tdof,1);
FGC(6)=P1;


FG=FG+FGC;


%% Applications of BCs

for kk=[8,7,4,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    FG(kk,:)=[];
end
KG;
FG;
UG=inv(KG)*FG

