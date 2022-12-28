%Static Analysis of a Static Frame
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
Cood=[0 0;1 0;1 1;0 1]; %coordinates of nodes
Ae=0.5;
He=0.02;
V=0.25;
rho=[7.324e-4];

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
Ke=zeros(dofe,dofe);
MG=zeros(tdof,tdof);
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
    Me=rho*He*Ae/12*[2 0 1 0 1 0;0 2 0 1 0 1;1 0 2 0 1 0;0 1 0 2 0 1;1 0 1 0 2 0;0 1 0 1 0 2];

    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
            MG(CONN(i,j),CONN(i,k))=MG(CONN(i,j),CONN(i,k))+Me(j,k);

        end
    end
end

KG
MG



%% Applications of BCs

for kk=[8,7,4,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    MG(kk,:)=[];
    MG(:,kk)=[];
end

[v,d]=eig(KG,MG)
[omega,index]=sort(sqrt(diag(d)));
[frequency]=omega/(2*3.14)