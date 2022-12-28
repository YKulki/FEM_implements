%Static Analysis of a Static Frame
clc 
clear 

%% Geometry and Material Properties
E=[200e9];          %Moduli
I=[20e6]*1e-12;
A=[4e4]*1e-6;
rho=[7.324e-4];
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
KeG=zeros(dofe,dofe);
MG=zeros(tdof,tdof);
MeG=zeros(dofe,dofe);
Te=zeros(dofe,dofe);

for i=1:ne
    m=cosd(theta(i));
    n=sind(theta(i));
    KeL=E*I/le(i)^3*[A*le(i)^2/I 0 0 -A*le(i)^2/I 0 0;0 12 6*le(i) 0 -12 6*le(i);0 6*le(i) 4*le(i)^2 0 -6*le(i) 2*le(i)^2;-A*le(i)^2/I 0 0 A*le(i)^2/I 0 0;0 -12 -6*le(i) 0 12 -6*le(i);0 6*le(i) 2*le(i)^2 0 -6*le(i) 4*le(i)^2];       
    a=rho*A*le(i)/6;
    b=rho*A*le(i)/420;
    MeL=[2*a 0 0 a 0 0 ; 0 156*b 22*b*le(i)^2 0 54*b -13*b*le(i);0 22*b*le(i)^2 4*b*le(i)^2 0 13*b*le(i) -3*b*le(i)^2;a 0 0 2*a 0 0;0 54*b 13*b*le(i) 0 156*b -22*b*le(i);0 -13*b*le(i) -3*b*le(i)^2 0 -22*b*le(i) 4*b*le(i)^2];
    Te=[m n 0 0 0 0;-n m 0 0 0 0;0 0 1 0 0 0;0 0 0 m n 0;0 0 0 -n m 0;0 0 0 0 0 1];
    KeG=transpose(Te)*KeL*Te;
    MeG=transpose(Te)*MeL*Te;
    
    for j=1:dofe
        for k=1:dofe
            KG(CONN(i,j),CONN(i,k))=KG(CONN(i,j),CONN(i,k))+KeG(j,k);
            MG(CONN(i,j),CONN(i,k))=MG(CONN(i,j),CONN(i,k))+MeG(j,k);

        end
    end
end
KG
MG


%% Applications of BCs( u1,v1,v3=0)

for kk=[12,11,10,3,2,1]
    KG(kk,:)=[];
    KG(:,kk)=[];
    MG(kk,:)=[];
    MsG(:,kk)=[];
end

[v,d]=eig(KG,MG)
[omega,index]=sort(sqrt(diag(d)));
[frequency]=omega/(2*3.14)

