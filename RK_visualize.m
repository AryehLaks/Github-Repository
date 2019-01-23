function [t1,y1] = RK_visualize(func,y0,t0,h,sq,sl)
%RK_visualize Creates an intuitive visualization of the Classical Fourth-Order Runge-Kutta Method
%
%RK_visualize(func,y0,t0,h)
%[t1,y1] = RK_visualize(func,y0,t0,h)
%[t1,y1] = RK_visualize(func,y0,t0,h,sq,sl)
%
%   This program is designed to help one learn the Classical Fourth-Order
%   Runge-Kutta Method. Note that the program only displays one timestep.
%   (Created by Aryeh Laks)
%
%Input Arguments
%   func = function, where dy/dt = f(t,y)
%   y0 = y_i, the initial condition
%   t0 = t_i, the initial time value
%   h = the desired timestep
%
%Output Arguments
%   y1 = y_(i+1), the next value of y created by the Runge-Kutta method
%   t1 = t_(i+1), the next value of t (in other words, t1 = t0 + h)
%
%Options. Note that RK_visualize(func,y0,t0,h) will run as usual and
%display a slope-field and a legend even though no values are entered into
%the spaces containing sq and sl
%
%   If the space containing sq = 'no', no slope-field will be created. This
%   may be useful when quiver, which is the primary function used to create
%   the slope field, will not perform well. When using this option, a
%   warning "Ignoring extra legend entries" will be displayed on the
%   command window.
%   Example: RK_visualize(func,y0,t0,h,'no') will not display a
%   slope-field.
%
%   If the space containing sl = 'no', no legend will be displayed. This
%   may be useful if the legend is obscuring useful information.
%   Example: RK_visualize(func,y0,t0,h,0,'no') will not display a legend.

close all
% "Actual" solution of the differential equation
[t,y]=ode45(func,[t0,h*1.5],y0); % h*1.5 is used to make the graph look nicer
range=[min(y),max(y)]; % used later to correctly size the quiver plot height
plot(t,y,'y')
% Overall graphics
whitebg('k')
xlabel('t')
ylabel('y')
title('\color{black}Visualization of Runge-Kutta')
hold on

% Visualization of the 4 slopes used to create Runge-Kutta

t=t0;
y=y0;

% 1st slope (m_1)
m1=func(t,y); % 1st slope
plot(t,y,'og') % Vis. for "(t,y)=(t_i,y_i)"
% Vis. for "m_1=f(t_i, y_i)"
tplot=linspace(t,h/2,20);
yplot=m1.*tplot+y;
range(end+1:end+2)=[min(yplot);max(yplot)]; % for quiver
plot(tplot,yplot,'-.g')

% 2nd slope (m_2)
m2=func(t+h/2,y+h/2*m1); % 2nd slope
plot(t+h/2,y+h/2*m1,'or') % Vis. (visualization) for "(t,y)=(t_i+h\2,y_i+h\2*m_1)"
% Vis. for "m_2 @ (t_i, y_i) => y=y_i+t*m_2"
tplot=linspace(t,h/2,20);
yplot=m2.*tplot+y;
plot(tplot,yplot,'-.r')
range(end+1:end+2)=[min(yplot);max(yplot)]; % for quiver
% Vis. for "m_2=f(t_i+h\2,y_i+h\2*m_1)"
plot(tplot+h/2,yplot+h/2*m1,'--r')
range(end+1:end+2)=[min(yplot)+h/2*m1;max(yplot)+h/2*m1]; % for quiver

% 3rd slope (m_3)
m3=func(t+h/2,y+h/2*m2); % 3rd slope
plot(t+h/2,y+h/2*m2,'oc') % Vis. for "(t,y)=(t_i+h\2,y_i+h\2*m_2)"
% Vis. for "m_3 @ (t_i, y_i) => y=y_i+t*m_3"
tplot=linspace(t,h,20);
yplot=m3.*tplot+y;
plot(tplot,yplot,'-.c')
% Vis. for "m_3=f(t_i+h\2,y_i+h\2*m_2)
tplot(12:end)=[];
yplot(12:end)=[];
range(end+1:end+2)=[min(yplot);max(yplot)]; % for quiver
plot(tplot+h/2,yplot+h/2*m2,'--c')
range(end+1:end+2)=[min(yplot)+h/2*m2;max(yplot)+h/2*m2]; % for quiver

% 4th slope (m_4)
m4=func(t+h,y+h*m3); % 4th slope
plot(t+h,y+h*m3,'om') % Vis. for "(t,y)=(t+h,y+h*m3)
% Vis. for "m_4 @ (t_i, y_i) => y=y_i+t*m_4"
tplot=linspace(t,h/2,20);
yplot=m4.*tplot+y;
plot(tplot,yplot,'-.m')
range(end+1:end+2)=[min(yplot);max(yplot)]; % for quiver
% Vis. for "m_4=f(t_i+h,y_i+h*m_3)
plot(tplot+h,yplot+h*m3,'--m')
range(end+1:end+2)=[min(yplot)+h*m3;max(yplot)+h*m3]; % for quiver

% Classical 4th-Order Runge-Kutta slope

ode45m=m1./6+m2./3+m3./3+m4./6; % Runge-Kutta slope
% Vis. for "(t,y)=(t_(i+1),y_(i+1))"
plot(t+h,y+h.*ode45m,'*w','MarkerSize',15)
% Vis. for "y_(i+1)=y_i+h*[m_1/6+m_2/3+m_3/3+m_4/6]"
tplot=linspace(t,h,20);
yplot=ode45m.*tplot+y;
plot(tplot,yplot,'w')

% Creation of slope-field by using quiver
if ~exist('sq','var')
    sq = 0;
end
if sq ~= 'no'
    H=h/8; % controls the spacing of the slope-field
    x=0:H:h*3/2;
    y=linspace(min(range),max(range),length(x));
    [X,Y]=meshgrid(x,y);
    % dy/dt=f(t,y)/1 => dy=f(t,y), dt=1
    v=func(X,Y); % dy=f(t,y)
    u=ones(size(v)); % dt=1
    quiver(X,Y,u,v,'b:')
end

% Legend
if ~exist('sl','var')
    sl = 0;
end
if sl ~= 'no'
    lgd = legend('"actual" value (uses ode45)',...
        '(t, y) = (t_i, y_i)','m_1 = f(t_i, y_i)',...
        '(t, y) = (t_i + h \div 2, y_i + h \div 2 \bullet m_1)','m_2 @ (t_i, y_i) \Rightarrow y = y_i + t \bullet m_2','m_2 = f(t_i + h \div 2, y_i + h \div 2 \bullet m_1)',...
        '(t, y) = (t_i + h \div 2, y_i + h \div 2 \bullet m_2)','m_3 @ (t_i, y_i) \Rightarrow y = y_i + t \bullet m_3','m_3 = f(t_i + h \div 2, y_i + h \div 2 \bullet m_2)',...
        '(t, y) = (t_i + h, y_i + h \bullet m_3)','m_4 @ (t_i, y_i) \Rightarrow y = y_i + t \bullet m_4','m_4 = f(t_i + h, y_i + h \bullet m_3)',...
        '(t, y) = (t_{i + 1}, y_{i + 1})','y_{i + 1} =  y_i + h \bullet [m_1/6 + m_2/3 + m_3/3 + m_4/6]',...
        'slope-field');
    lgd.FontSize = 7;
end

% Outputs
t1=t0+h; % the next value of t
y1=y0+h.*ode45m; % the next value of y created by the Runge-Kutta method
end