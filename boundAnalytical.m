%numerical vs analytical solutions
clear;

%% Timespan
%calculate timespan for simulation (non dim
r = 6.944e-4;
max_t = 7*r;%time for 7 mins
tspan = [0, max_t];

%% Numerical Solution

%Initial conditions
y0_reduced = [1 1 0 1 0 0];
%Numerical solution
[t_reduced,y_reduced] = ode45(@reducedSystem, tspan, y0_reduced);

%% Analytical solution

V0 = 1;
ka = 0.8847;
F0 = 1;
eta = 0.1037;
l = 1200;
df = 43.2;
beta = ka*eta*l*V0 +df;
kin = 881;
c1 = 44.8;
c2 = 44.8;
R = c1 + c2;

C = F0-(R/beta);
C3 = -ka*l*V0*((R/(kin*beta)) + (C/(kin-beta)));

time = linspace(0,max_t);
B = zeros(1, length(t_reduced));
%for i = 1: length(time)
for i = 1: length(t_reduced)
    t = t_reduced(i);
    B0_t = ka*l*V0*((R/(kin*beta)) + (C*exp(-t*beta)/(kin-beta)))+C3*exp(-kin*t);
    B(i) = B0_t;
end


%% Convert back to dimensional

t_reduced = t_reduced/r;

l=1200;
V0 = 1e10;
y_reduced = y_reduced*V0/l;%BPVs per cm^3
B = B*V0/l;%BPVs per cm^3


%% Plot

hold on;
plot(t_reduced, y_reduced(:, 3));
scatter(t_reduced, B);
xlabel('Time (mins)');
ylabel({'Bound Polymersomes'; '(Polymersomes cm^{-3})'});
hold off;
legend('Numerical', 'Analytical', 'Location','northwest');
set(gca,'fontsize', 18); 