clear;

%% Timespan
%calculate timespan for simulation (non dim)
r =6.944e-4;
max_t = 30*60*r;%time for 30 hours mins
tspan = [0, max_t];

%% ODES

%%%%%%%%%%%%%%%%%%%%%%%%%
%multiple binding
%order of equations:
%1. Cells
%2. Free BPVs
%3. BPVs with one bond 
%4-21. BPVs with multiple bonds
%22. BPVs with n bonds
%23. Receptors
%24. Internalised BPVs
%25. Released drug
n_minus2 = 18; %i.e n is 20
y0_multiple = [1 1 0];
y0_bound = zeros(1, n_minus2);
y0_multiple = [y0_multiple y0_bound 0 1 0 0];
[t,y_multiple] = ode45(@multiplebinding, tspan, y0_multiple);

%% Convert back to dimensional units

%convert time  from non dim to hours
%convert to min first
t = t/r;
%convert to hours
t = t/60;

%convert cells from non dim to density
m = y_multiple(:,1);
m = m *1e8;%cells per cm^3

%convert Free BPVs
V0 = 1e10;
V = y_multiple(:,2);
V = V* V0; %BPVs per cm^3

%convert BPV bound with 1 bond
l=1200;
B1 = y_multiple(:,3);
B1 = B1*V0/l;%BPVs per cm^3

%convert F
F = y_multiple(:,23);
K = 1e8;
f0 = 1.6e-20;
NA = 6.022e23;
F = F*f0*K*NA; %receptors per cm3

%convert Bin
Bin = y_multiple(:,24);
Bin = Bin*V0/l;%BPVs per cm^3

%convert P
P0 = 4.15e-14;
P = y_multiple(:,25);
P = P*P0*K;


%% Plot

subplot(3, 2, 1);
hold on;
plot(t, m);
hold off;
xlabel('Time (hours)');
ylabel({'Cell density'; '(cells cm^{-3})'});
set(gca,'fontsize', 18); 

subplot(3, 2, 2);
hold on;
plot(t, F);
hold off;
xlabel('Time (hours)');
ylabel({'Free Receptors'; '(Receptors cm^{-3})'});
set(gca,'fontsize', 18); 

subplot(3, 2, 3);
hold on;
plot(t, V);
hold off;
xlabel('Time (hours)');
ylabel({'Free'; 'polymersomes'; '(polymersomes cm^{-3})'});
set(gca,'fontsize', 18); 


subplot(3, 2, 4);
hold on;
plot(t, B1);
hold off;
xlabel('Time (hours)');
ylabel({'Polymersomes with'; 'one bond'; '(Polymersomes cm^{-3})'});
set(gca,'fontsize', 18); 

subplot(3, 2, 5);
hold on;
plot(t, Bin);
hold off;
xlabel('Time (hours)');
ylabel({'Internalised'; 'Polymersomes'; '(Polymersomes cm^{-3})'});
set(gca,'fontsize', 18); 

subplot(3, 2, 6);
hold on;
plot(t, P);
hold off;
xlabel('Time (hours)');
ylabel({'Released drug';'(mol cm^{-3})'});
set(gca,'fontsize', 18); 
