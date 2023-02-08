% clc
clear

m=40;
n=200;

b = 1e-3; %meter
L= 0.5; %meter

dy = b/m ; %size of elements in y direction
dz = L/n; %size of elements in z direction

y = -(0:dy:b);
z = 0:dz:L;

mu = 24e-3; %fluid viscosity
rho = 1017; %fluid density
theta = 30; %degree of wall inclination 
g = 9.8; %g constant

v_max = rho*g*b^2*cosd(theta)/2/mu; %maximum velocity of film
v = v_max*(1-(y/b).^2);

CA_Y0 = 0.1;
CA_Z0 = 0;

CB_Y0 = 1;
CB_Z0= 1;

CC_Y0 = 0;
CC_Z0 = 0;

D = 1e-9;

k=0.01;

K = 1000;
%---------------------------------%

La = D*dz/dy^2;

A = zeros((m-1)*(n-1));
B = zeros((m-1)*(n-1),1);

for j=1:(m-1)*(n-1)
    if ceil(j/(n-1))== n-1
        A(j,j) = -(1+La/v(ceil(j/(n-1))));
    else
        A(j,j) = -(1+2*La/v(ceil(j/(n-1))));
    end
    if j+n-1 < (m-1)*(n-1)+1
        A(j,j+n-1) = La/v(ceil((j)/(n-1)));
        A(j+n-1,j) = La/v(ceil((j+n-1)/(n-1)));
    end
    
    if rem(j,(n-1))~=0
        A(j+1,j) = 1;
    end
    
    if j<n
        B(j) = -0.1 * La/v(1);
    end
end
CA = linsolve (A,B);

CAF = zeros (m+1,n+1);
CAF(1,:)=0.1;

for i=1:m-1
    
        
    CAF (i+1,2:end-1) = CA(1+(i-1)*(n-1):i*(n-1));
    
end
CAF (:,end) = CAF(:,end-1);
CAF (end,:) = CAF(end-1,:);

% mean(CAF(:,end-1))
CCF =zeros(m+1,n+1);
sums=0;
sumt=0;
for ii=1:m-1
    sums = sums + dy/3 * ((CAF(ii,end)+CCF(ii,end))*v(ii)+4*((CAF(ii+1,end)+CCF(ii+1,end))*v(ii+1))+(CAF(ii+2,end)+CCF(ii+2,end))*v(ii+2));
%     sumt = sumt + dy*((CAF(ii,end)+CAF(ii+1,end))/2+(CCF(ii,end)+CCF(ii+1,end))/2)*(v(ii)+v(ii+1))/2;
end
Fluxs = sums/L
% Fluxt = sumt/L