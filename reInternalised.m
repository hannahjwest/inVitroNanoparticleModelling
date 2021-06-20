%numerical vs analytical solutions
clear;

%% Timespan

%calculate timespan for simulation (non dim
r = 6.944e-4;
max_t = 60*r;%time for 1hour
tspan = [0, max_t];
%% Numerical Solution

%Initial conditions
y0_reduced = [1 1 0 1 0 0];
%Numerical solution
[t_reduced,y_reduced] = ode45(@reducedSystem, tspan, y0_reduced);

%% Analytical solution

%Parameters
V0 = 1;% 
ka = 0.8847;
F0 = 1;
eta = 0.1037;
l = 1200;
df = 43.2;
epsilon = 0.1;
c1 = 44.8;
c2 = 44.8;
beta = ka*eta*l*V0 +c1;
kin = 881;
db = 3.3122;

%Caclulate variables
R = c1+c2;
C = F0-(R/beta);
C3 = -ka*l*V0*((R/(kin*beta)) + (C/(kin-beta)));
gamma1 = (ka*V0*R/beta);
gamma2 = kin*ka*V0*C/(kin - beta);
gamma3 = kin*C3/l;
C4 =-( (gamma1/db)+(gamma2/(db - beta))+ (gamma3/(db-kin)));

time = linspace(0,max_t);
Bin = zeros(1, length(t_reduced));
%Calculate analytical solution
for i = 1: length(t_reduced)
    t = t_reduced(i);
    
    a = (gamma1/db);
    
    b = ((gamma2*exp(-beta*t))/(db-beta));
     
    c = ((gamma3*exp(-kin*t))/(db-kin));
    
    e=exp(-db*t);
    d = e*C4;
    
    Bin0_t =(a+b+c+d);

    Bin(i) = Bin0_t;
end

%% Convert back to dimensional
t_reduced = t_reduced/r;
time = time/r;

l=1200;
V0 = 1e10;
y_reduced = y_reduced*V0/l;%BPVs per cm^3
Bin = Bin*V0/l;%BPVs per cm^3


%% Relative error
relative_error = abs(((Bin- y_reduced(:,5)')./y_reduced(:,5)')*100);

%% Plot
hold on;
scatter(t_reduced, relative_error);
xlabel('Time (mins)');
xlim([0, 60]);
ylabel('Relative error (%)');
hold off;
set(gca,'fontsize', 18); 
set(gca,'fontsize', 18); 