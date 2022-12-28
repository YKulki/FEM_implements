%Static Analysis of a truss planar
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
A=[30e-4];          %Areas
rho=[7.324e-4];
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

le
theta

%% Assembly of Element Matrices
KG=zeros(tdof,tdof);
KeG=zeros(dofe,dofe);
MG=zeros(tdof,tdof);
MeG=zeros(dofe,dofe);
Te=zeros(nne,dofe);

for i=1:ne
    m=cosd(theta(i));
    n=sind(theta(i));
    KeL=E*A/le(i)*[1 -1; -1 1];  
    MeL=rho*A*le(i)/6*[2 1;1 2];
    Te=[m n 0 0; 0 0 m n];
    KeG=transpose(Te)*KeL*Te;
    MeG=transpose(Te)*MeL*Te;    
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+KeG(j,k);
            MG(CONN(i,j),CONN(i,k))=MG(CONN(i,j),CONN(i,k))+MeG(j,k);

        end
    end
end


%% Applications of BCs( u1,v1,u2,v2=0)

for kk=[4,3,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    MG(kk,:)=[];
    MG(:,kk)=[];
end

[v,d]=eig(KG,MG)
[omega,index]=sort(sqrt(diag(d)));
[frequency]=omega/(2*3.14)