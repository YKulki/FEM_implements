%Static Analysis of a bar
clc 
clear 

%% Geometry and Material Properties
L=[50e-2];                %length of bar
E=[200e9];     %Moduli
A=[32.5e-4 22.5e-4]; %Areas
P=5e3;            %Concentrated Load
qo=0;                %uniform load

%% Meshing
ne=2;
nne=2;
nn=ne*nne-ne+1;
dof=1;
dofe=dof*nne;
tdof=nn*dof;
le=10e-4;
CONN=zeros(ne,dofe);
k=1;
for i=1:ne
    for j=1:dofe
        CONN(i,j)=k;
        k=k+1;
    end
    k=k-1;
end

%% Assembly of Element Matrices

KG=zeros(tdof,tdof);
FGU=zeros(tdof,1);

for i=1:ne
    Ke=E*A(i)/le*[1 -1; -1 1];       
    Fe=qo*le/2*[1;1];
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
        end
        FGU(CONN(i,j),1)=FGU(CONN(i,j),1)+Fe(j,1);
    end
end
FGC=zeros(tdof,1);
nodeno=2;
FGC(nodeno)=P;

FG=FGU+FGC;

KGr=KG
FGr=FG


%% Applications of BCs( u1,u3=0)

for kk=[1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    FG(kk,:)=[];
end


UG=inv(KG)*FG
UGr=[0;UG(1,1);UG(2,1)]

KGr
UGr
FGr

R=zeros(tdof,1);
R=KGr*UGr - FGr

        