% Aryeh Leibish Laks
%
% Depicts the path of a projectile being affected by the gravity of the
% Earth and the Moon using Euler's Method. It also contains an if-structure
% that changes the size of the time-step depending on the speed of the
% projectile. This improves the accuracy when the projectile is moving
% quickly and eliminates redundancy when the object is moving slowly

clc,clear all,close all

tic

% Independant Variables
vtotal=11100; % Initial velocity vector of the projectile, m/s
angle=39; % angle from the positive x-axis, degrees
s=500; % time-step variator (should be in the same 10s-place as dmode - the mode of the differential distance)
N=1500000; % total number of iterations

% Given celestial distances
r_E=6.371*10^6; % radius of earth, meters
r_M=1.737*10^6; % radius of moon, meters
d=384.4*10^6; % distance from center of earth to center of moon in meters
fprintf('Assuming the radius of the Earth is %g million meters,\nthe radius of the Moon is %g million meters,\nthe distance from the Earth to the Moon is %g million meters,\n',r_E/10^6,r_M/10^6,d/10^6)

% 1st graph

subplot(2,1,1)
% Model of the Earth
n=100;
xE=linspace(-r_E,r_E,n);
yE1=sqrt(r_E^2-xE.^2);
yE2=-sqrt(r_E^2-xE.^2);
plot(xE,yE1,'r',xE,yE2,'r')

hold on
% Model of the Moon
xM1=linspace(-r_M+d,r_M+d,n);
xM2=linspace(-r_M,r_M,n);
yM1=sqrt(r_M^2-xM2.^2);
yM2=-sqrt(r_M^2-xM2.^2);
plot(xM1,yM1,'g',xM1,yM2,'g')

% Initial position vector of the Earth and the Moon
rE(1,:)=[0,r_E]; % meters
rM(1,:)=rE(1,:)-[d,0]; % meters
rEnorm(1)=norm(rE); % meters. Extra for 3rd graph

% Turning the initial velocity vector of the projectile into components
vnorm(1)=vtotal; % meters/seconds
vx=vtotal*cosd(angle); % meters/seconds
vy=vtotal*sind(angle); % meters/seconds
v(1,:)=[vx,vy]; % meters/seconds

% Constant masses and gravitational constant
G=6.673*10^-11; % universal gravitational constant, m^3*kg^-1*s^-2
m_Earth=5.972*10^24; % mass of earth, kg
m_Moon=7.348*10^22; % mass of moon, kg (7.34767309)
fprintf('the mass of the Earth is %g kg,\nthe mass of the Moon is %g kg,\nG is %g m^3*kg^-1*s^-2.\n',m_Earth,m_Moon,G)

fprintf('\nThe initial velocity is %g m/s at %g degrees from the positive x-axis,\n%g million meters from the Earth"s center...\n\n',vtotal,angle,norm(rE(1,:))/10^6)

% Empty matrix to speed up computation
v(2:N+1,2)=zeros;
rE(2:N+1,2)=zeros;
rM(2:N+1,2)=zeros;
vnorm(2:N+1)=zeros;
rEnorm(2:N+1)=zeros;

% Variable time-steps
h1=.0001*s; h2=.001*s; h3=.01*s; h4=.1*s; h5=1*s;
h6=10*s; h7=100*s; h8=1000*s; h9=10000*s; h10=100000*s;
% For statistics
H=[h1,h2,h3,h4,h5,h6,h7,h8,h9,h10];

% Keeps track of the variable time-steps in order to calculate the total amount of seconds
j1=0; j2=0; j3=0; j4=0; j5=0;
j6=0; j7=0; j8=0; j9=0; j10=0;

% Euler's Method Loop
for i=1:N
    % Helps one not give up hope in reaching the end of the loop
    if i==N/2
        disp('You have reached the half-way point in the loop. Thank you for your patience.')
    end
    
    % Determines which variable time-step to use
    if vnorm(i)>=10000
        h=h1; % (if the speed is below 100,000) the distance-step is ~ 1*s
        j1=j1+1; % used later to calculate the total time and the most common time-step
    elseif vnorm(i)>=1000 & vnorm(i)<=10000
        h=h2; j2=j2+1;
    elseif vnorm(i)>=100 & vnorm(i)<=1000
        h=h3; j3=j3+1;
    elseif vnorm(i)>=10 & vnorm(i)<=100
        h=h4; j4=j4+1;
    elseif vnorm(i)>=1 & vnorm(i)<=10
        h=h5; j5=j5+1;
    elseif vnorm(i)>=.1 & vnorm(i)<=1
        h=h6; j6=j6+1;
    elseif vnorm(i)>=.01 & vnorm(i)<=.1
        h=h7; j7=j7+1;
    elseif vnorm(i)>=.001 & vnorm(i)<=.01
        h=h8; j8=j8+1;
    elseif vnorm(i)>=.0001 & vnorm(i)<=.001
        h=h9; j9=j9+1;
    elseif vnorm(i)>=.00001 & vnorm(i)<=.0001
        h=h10; j10=j10+1;
    end
    
    % Euler's Method Executed
    v(i+1,:)=v(i,:)+h.*(-rE(i,:).*G.*m_Earth./(norm(rE(i,:)).^3)-rM(i,:).*G.*m_Moon./(norm(rM(i,:)).^3)); % m/s
    vnorm(i+1)=norm(v(i+1)); % m/s, for 2nd graph (and 3rd graph)
    rE(i+1,:)=rE(i,:)+h.*v(i,:); % distance vector from the earth, m
    rM(i+1,:)=rE(i+1,:)-[d,0]; % distance vector from the moon, m
    rEnorm(i+1)=norm(rE(i+1,:)); % m. Extra for 3rd graph% supposed to stop the projectile from shooting off after hitting the center of the earth or moon
    % supposed to stop the projectile from shooting off after hitting the center of the earth or moon
    if norm(rE(i+1,:))<0 || norm(rM(i+1,:))<0
        break % didn't work
    end
end

% 1st graph continued
x=rE(:,1); % x-component of the position vector, m
y=rE(:,2); % y-component of the position vector, m
plot(x,y) % no point of using comet if the time-steps are not constant
legend('Earth','Earth','Moon','Moon','Projectile Path')
xlabel('x (meters)')
ylabel('y (meters)')
title('Path of a Projectile Being Affected by the Gravity of the Earth and the Moon')

% 2nd graph

subplot(2,1,2)
plot3(x,y,vnorm)
xlabel('x (m)')
ylabel('y (m)')
zlabel('velocity (m/s)')
title('Position Vector vs. Velocity')

% 3rd graph

figure
semilogy(rEnorm,vnorm)
axis([0 max(rEnorm)*65/64 100 max(vnorm)*5/4])
xlabel('absolute distance from the earth (m)')
ylabel('velocity (m/s)')
title('Absolute Distance from the Earth vs. Velocity')

% Statistics

% Total amount of time
sec=j1*h1+j2*h2+j3*h3+j4*h4+j5*h5+j6*h6+j7*h7+j8*h8+j9*h9+j10*h10;
days=sec/(60^2*24);
% Most and least common variable time-step
J=[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10];
Hmax=H(find(J==max(J))); % seconds
Hmin=H(find(J==min(J))); % seconds
Hmin(2:end)=[];
% Average time
t_avg=sec/N; % seconds

% Total and average distance
displacement=abs(diff(rEnorm)); % meters
dtotal=sum(displacement); % meters
d_avg=dtotal/N; % meters
dmode=mode(displacement); % meters

% Maximum velocity and its location
vmax=max(vnorm); % meters/seconds
I=find(vnorm==vmax);
rvmax_E=norm(rE(I,:)); % meters
rvmax_M=norm(rM(I,:)); % meters
% Minimum velocity and its location
vmin=min(vnorm); % meters/seconds
I=find(vnorm==vmin);
rvmin_E=norm(rE(I,:)); % meters
rvmin_M=norm(rM(I,:)); % meters
% Average velocity and mode
v_avg=dtotal/sec; % meters/seconds
vmode=mode(vnorm); % meters/seconds

fprintf('\n\nThe time-step (h) changes, the average time-step is %g seconds.\n',t_avg)
fprintf('The total amount of time is %g seconds, or %g days.\n',sec,days)
fprintf('The most common time-step was %g seconds and the least common time-step was %g seconds.\n',Hmax,Hmin)

fprintf('\nThe total distance traveled is %g million meters.\n',dtotal/10^6)
fprintf('The average distance between time-steps is %g meters and the mode is %g meters.\n',d_avg,dmode)

fprintf('\nThe maximum velocity is %g m/s at %g million meters away from the Earth and %g million meters away from the Moon.\n',vmax,rvmax_E/10^6,rvmax_M/10^6)
fprintf('The minimum velocity is %g m/s at %g million meters away from the Earth and %g million meters away from the Moon.\n',vmin,rvmin_E/10^6,rvmin_M/10^6)
fprintf('The average velocity is %g meters/second and the mode is %g meters/second.\n',v_avg,vmode)

total_time=toc;
fprintf('\nIt took %g seconds for entire program to run.\n',total_time)