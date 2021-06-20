function dydt = reducedSystem(t,y)


%% parameters
dydt = zeros(6,1);
dm = 7.24;
mu = 1728;
ka = 0.8847;
l = 1200;
L = l*y(2);
kin = 881;
eta = 0.1037;
c1 = 44.8;
c2 = 44.8;
df = 43.2;
db = 3.3122;
v = 4;
dp = 1728;

%% ODES

%cells
dydt(1) = y(1)*(1-y(1)) - (mu + dm)*y(6);

%free BPVs
dydt(2) = -ka*y(2)*y(4);

%BPVs with one bond
dydt(3) = ka*L*y(4)-kin*y(3);

%Receptors
b = y(3)/y(1);
R = c1+c2;
dydt(4) = -eta*ka*y(4)*L-df*y(4)+R*y(1);

%internalised BPVS
dydt(5) = kin*y(3)/l - db*y(5);

%released drug
dydt(6) = v*db*y(5)-dp*y(6);


end