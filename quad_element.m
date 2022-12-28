%Static Analysis of a Static Frame
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
Cood=[0 0;1 0;2 0; 2 1; 2 2;1 2; 0 2;0 1; 1 1]; %coordinates of nodes
Ae=1;
he=0.02;
v=0.3;
P1=-30e3;
bx=0;
by=0;
Tx=[];
Ty=[];
ng=4;


%% Meshing
ne=4;
nne=4;
nn=9;       
dof=2;
dofe=dof*nne;
tdof=nn*dof;
CONN=zeros(ne,dofe);
NCONN=zeros(ne,nne);
NCONN=[1 2 9 8;2 3 4 9;4 5 6 9;6 7 8 9];

CONN=zeros(ne,dofe);
for i=1:ne
    for j=1:4
       CONN(i,2*j-1)=2*NCONN(i,j)-1;
       CONN(i,2*j)= 2*NCONN(i,j);
    end
end

%% Gauss quadrature

p=gauss_quad(ng);
Gp=p(1,:);
W=p(2,:);

D=(E/(1-v^2))*[1 v 0;v 1 0;0 0 (1-v)/2];

KG=zeros(tdof,tdof);
FG=zeros(tdof,1);

for i=1:ne
       de=zeros(8,1);

    for j=1:4
        de(2*j-1,1)=Cood(NCONN(i,j),1);
        de(2*j,1)=Cood(NCONN(i,j),2);
    end
    Ke=zeros(dofe,dofe);
    Feb=zeros(dofe,1);

    for j=1:ng
        for k=1:ng
            xi=Gp(j);    
            eta=Gp(k);

            phi1=0.25*(1-xi)*(1-eta);
            phi2=0.25*(1+xi)*(1-eta);
            phi3=0.25*(1+xi)*(1+eta);
            phi4=0.25*(1-xi)*(1+eta);

            psi1=[phi1 0 phi2 0 phi3 0 phi4 0];
            psi2=[0 phi1 0 phi2 0 phi3 0 phi4];
            
            psi1_xi=0.25*[-(1-eta) 0 (1-eta) 0 (1+eta) 0 -(1+eta) 0];
            psi2_xi=0.25*[0 -(1-eta) 0 (1-eta) 0 (1+eta) 0 -(1+eta)];
            psi1_eta=0.25*[-(1-xi) 0 -(1+xi) 0 (1+xi) 0 (1-xi) 0];
            psi2_eta=0.25*[0 -(1-xi) 0 -(1+xi) 0 (1+xi) 0 (1-xi)];

            J11=psi1_xi*de;
            J12=psi2_xi*de;
            J21=psi1_eta*de;
            J22=psi2_eta*de;
            
            J=J11*J22-J12*J21;

            Be=1/J*[J22 -J12 0 0;0 0 -J21 J11;-J21 J11 J22 -J12]*[psi1_xi;psi2_xi;psi1_eta;psi2_eta];

            Ke=Ke+he*transpose(Be)*D*Be*J*W(j)*W(k);
            Feb=Feb+he*(transpose(psi1)*bx+transpose(psi2)*by)*J*W(j)*W(k);
        end
    end
     for j=1:dofe
        for k=1:dofe
          KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+Ke(j,k);
        end
        FG(CONN(i,j),1)=FG(CONN(i,j),1)+Feb(j,1);
    end
end

KGr=KG;
FGC=zeros(tdof,1);
FGC(10)=P1;


FG=FG+FGC;
FGr=FG;


%% Applications of BCs

for kk=[16,15,14,13,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    FG(kk,:)=[];
end
KG;
FG;
UG=inv(KG)*FG;

UGr=zeros(tdof,1);
UGr(3:12)=UG(1:10);
UGr(17:18)=UG(11:12);
UGr

R=KGr*UGr-FGr

