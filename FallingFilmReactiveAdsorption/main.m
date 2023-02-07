clc
clear

m=40;  %Number of elements in y direction
n=200; %Number of elements in z direction

b = 1e-3; %Film Thickness, [m]
L= 0.5;   %Film Length, [m]

dy = b/m ; %size of elements in y direction
dz = L/n;  %size of elements in z direction

y = -(0:dy:b);
z = 0:dz:L;

mu = 24e-3; %fluid viscosity [Pa.s]
rho = 1017; %fluid density [kg/m^3]
theta = 30; %degree of wall inclination 
g = 9.8; %g constant [m/s^2]

v_max = rho*g*b^2*cosd(theta)/2/mu; %maximum velocity of the film [m/s]
v = v_max*(1-(y/b).^2); %velocity distribution in the film [m/s]

D = 1e-9; %Diffustion Coeff. [m^2/s]

k=0.01; %Kinetics First Coeff. [lit/mol.s]

K = 1000; % Kinetics Second Coeff. [lit/mol]
%---------------------------------%

La1 = D*dz/dy^2; %First Constant in the study
La2 = k*dz;      %Second Constant in the study
La3 = k*dz/K;    %Third Constant in the study

Theta = D*z/v_max/b^2; %Based on Akanksha's study
alpha = k*b^2/D;       %Based on Akanksha's study

A = zeros((m-1)*(n-1));
B = zeros((m-1)*(n-1),1);

Ap = zeros((m-1)*(n-1));
Bp = zeros((m-1)*(n-1),1);


CC0 = zeros ((m-1)*(n-1),1);

CC = zeros ((m-1)*(n-1),1);
s = 0;
while max(abs(CC0- CC))>1e-6 || s==0 %check whether the Hypothetical and the Real C_C are the same by an accuracy of 1e-6
    if s>0
        CC0 = CC; %make the hypothetical C_C equal to calculated C_C in the last loop
    end
    %solve C_A with hypothetical C_C
    for j=1:(m-1)*(n-1)
        if ceil(j/(n-1))== n-1
            A(j,j) = -(1+La1/v(ceil(j/(n-1)))+La2/v(ceil(j/(n-1)))*(1-CC0(j)));
        else
            A(j,j) = -(1+2*La1/v(ceil(j/(n-1)))+La2/v(ceil(j/(n-1)))*(1-CC0(j)));
        end
        
        if j+n-1 < (m-1)*(n-1)+1
            A(j,j+n-1) = La1/v(ceil((j)/(n-1)));
            A(j+n-1,j) = La1/v(ceil((j+n-1)/(n-1)));
        end
        
        if rem(j,(n-1))~=0
            A(j+1,j) = 1;
        end

        if j<n
            B(j) = -0.1 * La1/v(1)-La3/v(ceil(j/(n-1)))*(CC0(j));
        else
            B(j) = -La3/v(ceil(j/(n-1)))*(CC0(j));
        end
    end

    CA = linsolve (A,B);
    
    %solve C_C 
    for q=1:(m-1)*(n-1)
        
        if q<n || ceil(q/(n-1))== n-1
            Ap(q,q) = -(1+La1/v(ceil(q/(n-1)))) - La2*CA(q)/v(ceil(q/(n-1))) -La3/v(ceil(q/(n-1)));
        else
            Ap(q,q) = -(1+2*La1/v(ceil(q/(n-1)))) - La2*CA(q)/v(ceil(q/(n-1))) -La3/v(ceil(q/(n-1)));
        end
            
        if q+n-1 < (m-1)*(n-1)+1
            Ap(q,q+n-1) = La1/v(ceil(q/(n-1)));
            Ap(q+n-1,q) = La1/v(ceil((q+n-1)/(n-1)));
        end

        if rem(q,(n-1))~=0
            Ap(q+1,q) = 1;
        end

        
        if rem(q,n-1)==1
            Bp(q) = -La2/v(ceil(q/(n-1)))*CA(q);
        else
            Bp(q) = -La2/v(ceil(q/(n-1)))*CA(q);
        end
    end
    CC = linsolve(Ap,Bp);
    s=s+1 %corrector counter
end
CAF = zeros (m+1,n+1);
CAF(1,:)=0.1;

CCF = zeros (m+1,n+1);

CCF(:,1)=0;

for i=1:m-1
    
    CAF (i+1,2:end-1) = CA(1+(i-1)*(n-1):i*(n-1));
    CCF (i+1,2:end-1) = CC(1+(i-1)*(n-1):i*(n-1));
    
end

CAF (:,end) = CAF(:,end-1);
CAF (end,:) = CAF(end-1,:);

CCF (:,end) = CCF(:,end-1);
CCF (end,:) = CCF(end-1,:);
CCF(1,:)=CCF(2,:);

CBF = 1-CCF;

%(mean(CCF(:,end-1))+mean(CAF(:,end-1)))*mean(v(:))


sums=0;
% sumt=0;
for ii=1:m-1
    sums = sums + dy/3 * ((CAF(ii,end)+CCF(ii,end))*v(ii)+4*((CAF(ii+1,end)+CCF(ii+1,end))*v(ii+1))+(CAF(ii+2,end)+CCF(ii+2,end))*v(ii+2));
%     sumt = sumt + dy*((CAF(ii,end)+CAF(ii+1,end))/2+(CCF(ii,end)+CCF(ii+1,end))/2)*(v(ii)+v(ii+1))/2;
end
Fluxs = sums/L
% Fluxt = sumt/L