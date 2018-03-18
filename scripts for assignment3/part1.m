% ELEC4700 Assignment 3
% By Huanyu Liu 100986552
% Part1
clear;
clc;

% Initialize the parameters
n=50; % number of particles
T=300; % temperture of the backgound
L=200e-9; % length of the frame (figure 1)
H=100e-9; % height of the frame
tao=0.2e-12; % the given mean time between collisions
m0=9.109e-31; % mass of a particle
mn=0.26*m0; % effective mass
kb=1.38e-23; % constant coeffient
vth=sqrt(2*kb*T/mn); % average speed of each particle
con=1e11; % The electron concentration
% Initialize the positions of each particle
Pox = L*rand(1,n);
Poy = H*rand(1,n);

% New parameters for assignment3
V=0.1;
E=V/L;
e=1.60217662e-19;
F=E*e;
a=F/m0;

% Initialize the speed of each particle and measure the initial temperature
for num=1:n
Vx(num) = randn()*vth/sqrt(2);
Vy(num) = randn()*vth/sqrt(2);
end
% draw the first locations of the particles and the blocks
figure(1)
plot(Pox,Poy,'.');
xlim([0 L]);
ylim([0 H]);
hold on
% more parameters that will be used in the loop
TStop = 1e-12; % max running time
t=0; % start time
dt=1e-14; % step time
intervals=round(TStop/dt); % number of steps 
Vz=zeros(1,intervals); % initial the size of all changing speed (will be used in hist)
ddt = 0; % time since last timestop
collisions=0; % number of timestops
time=0; % initialize the duration between collisions
path=zeros(1,n); % initialize the size of path length
j=0; % before the model runs, the current is 0
current(1)=0;
 while t < TStop
     z=round(1+t/dt); % index, the z-th interval between collisions
     Pscat = 1-exp(-ddt/tao); % scattering posibility
     if Pscat > rand % if scatter
         time=time+ddt; % total time when scattering occur
         ddt=0; % reset the parameter for the possibility as required
         collisions=collisions+1; % one more collision occurs
         Vx = randn(1,n).*vth/sqrt(2);
         Vy = randn(1,n).*vth/sqrt(2); % velocity changes (in maxwell-boltzmann distribution)
         average_path_length(collisions)=sum(path)/n; % average path length for this interval
         path=zeros(1,n); % reset the path length
     else % nothing happens, same speed the next duration of time step
         path=path+sqrt(Vx.^2+Vy.^2).*dt; % add the next timestep's path length to the total path length
         ddt=ddt+dt; % add the timestep size to the parameter
     end
         Vact=sqrt(sum(Vx.^2+Vy.^2)/n); % the average speed of all the particles
         Vz(z)=Vact; % will be used to get the distribution in hist
         
     Vx = Vx + a.*t;
     tPx = Pox + Vx.*dt; % predict the position 
     tPy = Poy + Vy.*dt;
% when the particles go to the right and left border
     px1 = Pox >= L;
     Pox(px1) = Pox(px1) - L;
     px2 = Pox <= 0;
     Pox(px2) = Pox(px2) + L;
     
     py1 = tPy <= 0;
     Vy(py1) = Vy(py1) .* (-1);
     py2 = tPy >= H;
     Vy(py2) = Vy(py2) .* (-1);
         % now all velocity have been modified to the correct direction,
         % update the position
     PreviousPox = Pox;
     PreviousPoy = Poy;
     Pox = Pox + Vx.*dt;
     Poy = Poy + Vy.*dt;
     
     j=j+1;
     current(j)=H*e*con*sum(Vx)/n; % current formula
     figure(1)
     for i=1:n
     plot([PreviousPox(i),Pox(i)],[PreviousPoy(i),Poy(i)]);
     end
     xlim([0 L]);
     ylim([0 H]);
     hold on
     

     pause(0.01)
     t=t+dt;
     pre=current;
 end

     figure(2)
     recordt=0:dt:TStop;
     plot(recordt,current);
     xlabel('time');
     ylabel('current (A)');
     title('current vs. time');
 
 n1=Pox<0.2*L;
 n2=Pox<0.4*L;
 n3=Pox<0.6*L;
 n4=Pox<0.8*L;
 n5=Poy<0.25*H;
 n6=Poy<0.5*H;
 n7=Poy<0.75*H;
 Den=zeros(5,4);
 temper=zeros(5,4);
 Den(1,1)=sum(n1&n5);
 temper(1,1)=sum(Vx(n1&n5).^2+Vy(n1&n5).^2).*mn./(2*kb*Den(1,1));
 Den(1,2)=sum(n1&n6&(~n5));
 temper(1,2)=sum(Vx(n1&n6&(~n5)).^2+Vy(n1&n6&(~n5)).^2).*mn./(2*kb*Den(1,2));
 Den(1,3)=sum(n1&n7&(~n6));
 temper(1,3)=sum(Vx(n1&n7&(~n6)).^2+Vy(n1&n7&(~n6)).^2).*mn./(2*kb*Den(1,3));
 Den(1,4)=sum(n1&(~n7));
 temper(1,4)=sum(Vx(n1&(~n7)).^2+Vy(n1&(~n7)).^2).*mn./(2*kb*Den(1,4));
 Den(2,1)=sum((~n1)&n2&n5);
 temper(2,1)=sum(Vx((~n1)&n2&n5).^2+Vy((~n1)&n2&n5).^2).*mn./(2*kb*Den(2,1));
 Den(2,2)=sum((~n1)&n2&n6&(~n5));
 temper(2,2)=sum(Vx((~n1)&n2&n6&(~n5)).^2+Vy((~n1)&n2&n6&(~n5)).^2).*mn./(2*kb*Den(2,2));
 Den(2,3)=sum((~n1)&n2&n7&(~n6));
 temper(2,3)=sum(Vx((~n1)&n2&n7&(~n6)).^2+Vy((~n1)&n2&n7&(~n6)).^2).*mn./(2*kb*Den(2,3));
 Den(2,4)=sum((~n1)&n2&(~n7));
 temper(2,4)=sum(Vx((~n1)&n2&(~n7)).^2+Vy((~n1)&n2&(~n7)).^2).*mn./(2*kb*Den(2,4));
 Den(3,1)=sum((~n2)&n3&n5);
 temper(3,1)=sum(Vx((~n2)&n3&n5).^2+Vy((~n2)&n3&n5).^2).*mn./(2*kb*Den(3,1));
 Den(3,2)=sum((~n2)&n3&n6&(~n5));
 temper(3,2)=sum(Vx((~n2)&n3&n6&(~n5)).^2+Vy((~n2)&n3&n6&(~n5)).^2).*mn./(2*kb*Den(3,2));
 Den(3,3)=sum((~n2)&n3&n7&(~n6));
 temper(3,3)=sum(Vx((~n2)&n3&n7&(~n6)).^2+Vy((~n2)&n3&n7&(~n6)).^2).*mn./(2*kb*Den(3,3));
 Den(3,4)=sum((~n2)&n3&(~n7));
 temper(3,4)=sum(Vx((~n2)&n3&(~n7)).^2+Vy((~n2)&n3&(~n7)).^2).*mn./(2*kb*Den(3,4));
 Den(4,1)=sum((~n3)&n4&n5);
 temper(4,1)=sum(Vx((~n3)&n4&n5).^2+Vy((~n3)&n4&n5).^2).*mn./(2*kb*Den(4,1));
 Den(4,2)=sum((~n3)&n4&n6&(~n5));
 temper(4,2)=sum(Vx((~n3)&n4&n6&(~n5)).^2+Vy((~n3)&n4&n6&(~n5)).^2).*mn./(2*kb*Den(4,2));
 Den(4,3)=sum((~n3)&n4&n7&(~n6));
 temper(4,3)=sum(Vx((~n3)&n4&n7&(~n6)).^2+Vy((~n3)&n4&n7&(~n6)).^2).*mn./(2*kb*Den(4,3));
 Den(4,4)=sum((~n3)&n4&(~n7));
 temper(4,4)=sum(Vx((~n3)&n4&(~n7)).^2+Vy((~n3)&n4&(~n7)).^2).*mn./(2*kb*Den(4,4));
 Den(5,1)=sum((~n4)&n5);
 temper(5,1)=sum(Vx((~n4)&n5).^2+Vy((~n4)&n5).^2).*mn./(2*kb*Den(5,1));
 Den(5,2)=sum((~n4)&n6&(~n5));
 temper(5,2)=sum(Vx((~n4)&n6&(~n5)).^2+Vy((~n4)&n6&(~n5)).^2).*mn./(2*kb*Den(5,2));
 Den(5,3)=sum((~n4)&n7&(~n6));
 temper(5,3)=sum(Vx((~n4)&n7&(~n6)).^2+Vy((~n4)&n7&(~n6)).^2).*mn./(2*kb*Den(5,3));
 Den(5,4)=sum((~n4)&(~n7));
 temper(5,4)=sum(Vx((~n4)&(~n7)).^2+Vy((~n4)&(~n7)).^2).*mn./(2*kb*Den(5,4));
 [X,Y]=meshgrid(H/4:H/4:H,L/5:L/5:L);
 figure(3)
 subplot(1,2,1),surf(X,Y,Den);
 title('Density');
 subplot(1,2,2),surf(X,Y,temper);
 title('Temperature');
     fprintf(' The electrical field is: %g V/m\n',E);
     fprintf(' The force on each particle is: %g N\n',F);
     fprintf(' The acceleration of each particle is: %g m/s^2\n', a);