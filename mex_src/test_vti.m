

clear

g.xmin = -20;
g.zmin = 0;
g.dx = 1;
g.dz = 1;
g.nx = 40;
g.nz = 20;
g.nsx = 15;
g.nsz = 15;

xmax = g.xmin + g.nx*g.dx;
zmax = g.zmin + g.nz*g.dz;

Vp0 = 3000*ones(g.nx*g.nz, 1);
Vs0 = 2000*ones(g.nx*g.nz, 1);
epsilon = 0.5*ones(g.nx*g.nz, 1);
delta = -0.25*ones(g.nx*g.nz, 1);
gamma = 0.5*ones(g.nx*g.nz, 1);

r = 20;
dt= 5*pi/180;
theta = -pi/2:dt:pi/2;
x = r*sin(theta);
z = r*cos(theta);
Rx = [x(:) z(:)];


Tx = zeros(size(Rx));


[ttp, raisp] = ttcr2dvti_psv(Vp0, Vs0, epsilon, delta, g, Tx, Rx, true);
[ttsv, raissv] = ttcr2dvti_psv(Vp0, Vs0, epsilon, delta, g, Tx, Rx, false);
[ttsh, raissh] = ttcr2dvti_sh(Vs0, gamma, g, Tx, Rx);


t = theta*180/pi;
figure(1)
plot(t, 1000*ttp)
hold on 
plot(t, 1000*ttsv,'g')
plot(t, 1000*ttsh,'r')
hold off
legend('P','SV','SH')
xlabel('Angle (deg)')
ylabel('Temps (ms)')

