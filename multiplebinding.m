
function dydt = multiplebinding(t,y)


%% parameters

dydt = zeros(25,1);
%parameters
dm = 7.24;
mu = 1728;
P0 = 1;
phi = y(25)/y(1);
g = dm + (mu*phi/(P0+phi));
ka = 0.8847;
kd = 0.533;
l = 1200;
L = l*y(2);
kin = 881;
pl = 0.5;
pf = 0.5;
n = 20;
eta = 0.1037;
c1 = 44.8;
c2 = 44.8;
c3 = 83.28;
df = 43.2;
db = 3.3122;
v = 4;
dp = 1738;

%% ODES


%cells
dydt(1) = y(1)*(1-y(1)) - g*y(1);

%free BPVs
dydt(2) = -ka*y(2)*y(23) + kd*y(3)/l;

%BPVs with one bond
dydt(3) = ka*L*y(23) -kd*y(3)-kin*y(3) - ka*(l*pl - 1)*pf*y(3)*y(23) + 2*kd*y(4);

%BPVs with i bonds
j = 2;
for i = 4:21
    dydt(i) = ka*(pl*l -(j-1))*pf*y(23)*y(i-1) - j*kd*y(i)-kin*y(i)-ka*(pl*l-j)*pf*y(23)*y(i)+(j+1)*kd*y(i+1);
    j = j+1;
end
   
%BPV with n bonds
dydt(22) = ka*(pl*l-(n-1))*pf*y(23)*y(n-1)-kin*y(22)-n*kd*y(n);



%Receptors
sum1 = 0;
j = 1;
for i = 3:22
    sum1 = sum1 + (pl*l-j)*pf*y(23)*y(i);
    j = j+1;
end

sum2 = 0;
j = 1;
for i = 3:22
    sum2 = sum2 + j*y(i);
    j = j+1;
end

btot = sum2/y(1);
R = c1+(c2*btot/(c3+btot));
dydt(23) = -eta*ka*y(23)*L-eta*ka*sum1+eta*kd*sum2-df*y(23)+R*y(1);


%internalised BPVS
sum3 = 0;
for i = 3:22
    sum3 = sum3 + y(i);
end
dydt(24) = kin*sum3/l - db*y(24);

%released drug
dydt(25) = v*db*y(24)-dp*y(25);


end