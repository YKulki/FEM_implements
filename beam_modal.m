%Static Analysis of a Euler Beam
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
I=[2000]*1e-12;
A=[240 240]*1e-6;
rho=7840;
Cood=[0 0;300 0;600 0]*1e-3;  %coordinates of nodes
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
KeG=zeros(dofe,dofe);
MG=zeros(tdof,tdof);
MeG=zeros(dofe,dofe);

for i=1:ne
    KeG=E*I/le(i)^3*[12 6*le(i) -12 6*le(i);6*le(i) 4*le(i)^2 -6*le(i) 2*le(i)^2;-12 -6*le(i) 12 -6*le(i);6*le(i) 2*le(i)^2 -6*le(i) 4*le(i)^2];   
    MeG=rho*A(i)*le(i)/420*[156 22*le(i) 54 -13*le(i);22*le(i) 4*le(i)^2 13*le(i) -3*le(i)^2;54 13*le(i) 156 -22*le(i);-13*le(i) -3*le(i)^2 -22*le(i) 4*le(i)^2];
    Te=[cosd(theta(i)) sind(theta(i)) 0 0; 0 0 cosd(theta(i)) sind(theta(i))];
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+KeG(j,k);
            MG(CONN(i,j),CONN(i,k))=MG(CONN(i,j),CONN(i,k))+MeG(j,k);

        end
    end
end
KG
MG


%% Applications of BCs

for kk=[2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    MG(kk,:)=[];
    MG(:,kk)=[];

end

[v,d]=eig(KG,MG)
[omega,index]=sort(sqrt(diag(d)));
[frequency]=omega/(2*3.14)

firstmodeshape=zeros(3,1);
firstmodeshape(2,1)=v(1,1);
firstmodeshape(3,1)=v(3,1);


secondmodeshape(2,1)=v(1,2);
secondmodeshape(3,1)=v(3,2);

x=[0,10,15];
plot(x,secondmodeshape)



