%engmae146 final
load('finaldata.mat');
orbit.stm.mass=planet.mars.mass;
orbit.ste.mass=planet.earth.mass;
orbit.stm.periapsis=orbitpostransfer(orbit.stm,orbitpos(findr(orbit.stm,0),0));
orbit.ste.periapsis=orbitpostransfer(orbit.ste,orbitpos(findr(orbit.ste,0),0));

juliannow=juliandate(datetime(2032,1,1));

Marsanonow=caltruano(orbit.stm,juliannow);
Marsposnow=findr(orbit.stm,Marsanonow);
tratime=timebetween(orbit.stm,2.5,3.6);
%%
angle=findangle([1,1,0],[0,0,1])

%%
%get 3d vector r of an orbit
function r3d=orbvec3d(orbit,t)
trueanomaly=caltruano(orbit,t);
absr=findr(orbit,trueanomaly);
r2d=orbitpos(absr,trueanomaly);
r3d=orbitpostransfer(orbit,r2d);

end

%%
function r3d=orbvec3dano(orbit,trueanomaly)
absr=findr(orbit,trueanomaly);
r2d=orbitpos(absr,trueanomaly);
r3d=orbitpostransfer(orbit,r2d);

end
%%
function vpqw=findvpqw(orbit1,ano)
vpqw=[-sqrt(orbit1.u/orbit1.p)*sin(ano);sqrt(orbit1.u/orbit1.p)*(orbit1.e+cos(ano));0];
end


%%
%finding trueanomaly with orbit and time
function truanomaly=caltruano(orbit,timenow)
%orbit,timenow
e=orbit.e;
p=orbit.p;
a=orbit.a;
t=timenow-orbit.peritime;
M=2.*pi.*t./orbit.T;
f = @(x) x - e*sin(x)-M;
E=fzero(f,1);
truanomaly=2.*atan(sqrt((1+e)/(1-e)).*tan(E/2));
end

%%
%finding |r| with true anomaly
function r=findr(orbit,trueanomaly)
r=orbit.p/(1+orbit.e*cos(trueanomaly));
end 

%%
% %finding r position on orbit
function rpos=orbitpos(r,trueanomaly)
rpos=r*[cos(trueanomaly);sin(trueanomaly);0];
end



%%
%project r to plane
function vector2=projectionvec(vector1,plane)
vector2=vector1-(dot(vector1,plane)/dot(plane,plane))*plane;
end

%%
%finding angle
function theta=findangle(vector1,vector2)
theta=acos(dot(vector1,vector2)./norm(vector1)./norm(vector2));
end

%%
%time between 2 point
function t=timebetween(orbit,angle1,angle2)
angle=[angle1,angle2];
e=orbit.e;
T=orbit.T;
E=2.*atan(sqrt((1-e)/(1+e)).*tan(angle./2));
M=E-e.*sin(E);
time=M.*T./2./pi;
t=time(2)-time(1);
if t>0
    return
else
    t=T+t;
    return
end

end


%%
%hyper periapsis
%e=1+rp*velinfinity^2/u, rp(parking orbit decided)

%%
%transfer |r| into 3D vector

function rvector=orbitpostransfer(orbit,r)
Omega=orbit.raan;
i=orbit.inc;
w=orbit.aop;


R1 = [cos(Omega), -sin(Omega), 0;sin(Omega),  cos(Omega), 0;0,0,1];

R2 = [1,  0,   0;
      0, cos(i), -sin(i);
      0, sin(i),  cos(i)];

R3 = [cos(w), -sin(w), 0;
      sin(w),  cos(w), 0;
      0,      0,      1];

rvector= R1 * R2 * R3*r;
end

%%
%gravity field
%R_soi
function r=gravifield(r1,planetmass,majormass)
r=r1*(planetmass/majormass)^(2/5);
end

%%
%find luanch window
function [launchwin,orblist]=findnextlaunchwin(orbit1,orbit2,t,planet)
i=0;
timenow=t;
mass1=orbit1.mass;
mass2=planet.mass;
orbit.trans.u=planet.u;
oblist=[];
matrix1=[];
while i<1000
    anor1=caltruano(orbit1,timenow);
    r1=findr(orbit1,anor1);
    rsoi=gravifield(r1,mass1,mass2);
    ra=findr(orbit1,(anor1+pi));
    vecra=orbitpostransfer(orbit1,orbitpos(ra,anor1+pi));
    pjv=projectionvec(vecra,orbit2.plane);
    anor2=findangle(orbit2.periapsis,pjv);
    trans.evector=orbvec3d(orbit1,timenow);
    trans.nodevect=cross(trans.evector,[0;0;1]);
    trans.raan=findangle([1;0;0],trans.nodevect);
    trans.aop=findangle(trans.nodevect,trans.evector);
    r2=findr(orbit2,anor2);
    orbit.trans.a=(r1+r2+rsoi)/2;
    trans.e=(r2-r1-rsoi)/(2*orbit.trans.a);
    orbit.trans.T=findT(orbit.trans);
    orbit.trans.inc=findangle(vecra,pjv);
    transfertime=orbit.trans.T/2/3600/24;
    marsnow=caltruano(orbit2,timenow);

   
    planetarrivaltime=timebetween(orbit2,marsnow,anor2);
    
    phasetime=transfertime-planetarrivaltime;
    
    timenow=timenow+1;

    


    i=i+1;
    oblist(i,:)=[orbit.trans.a,trans.e,orbit.trans.inc,trans.raan,pi/2];
    matrix1(i,:)=[i,transfertime,planetarrivaltime,phasetime,findangle(vecra,pjv),anor2];

    
%    disp(orbit.trans)
%    disp([r1,r2,rsoi])
    
  

end

launchwin=matrix1;
orblist=oblist;



end    

%%
% find period
function T=findT(orbit)
u=orbit.u;
a=orbit.a;
T=(2*pi/sqrt(u))*a^(3/2);
end


%%
%r=orbitpostransfer(orbit.stm,orbitpos(Marsposnow,Marsanonow))
%%
[etmlwin,transorblist1]=findnextlaunchwin(orbit.ste,orbit.stm,juliannow,planet.sun);
[mtelwin,transorblist2]=findnextlaunchwin(orbit.stm,orbit.ste,juliannow,planet.sun);
etmrows=abs(etmlwin(:,4))<5;
disp(etmlwin(etmrows,:))
disp(transorblist1(etmrows,2:5))
mterows=abs(mtelwin(:,4))<1;
disp(mtelwin(mterows,:))
disp(transorblist2(mterows,2:5))
orblist1=transorblist1(etmrows,:);
winlist1=etmlwin(etmrows,:);

etmtransorb=orb(orblist1(1,1),orblist1(1,2),planet.sun.u,orblist1(1,3),orblist1(1,4),orblist1(1,5),(juliannow+winlist1(1,1)));
trans.v0=norm(findvpqw(etmtransorb,0));
earth.v=norm(findvpqw(orbit.ste,etmtransorb.raan+pi/2));
vinfi(1,1)=abs(earth.v-trans.v0);
trans.vt=norm(findvpqw(etmtransorb,pi));
mars.v=norm(findvpqw(orbit.stm,etmlwin(1,6)));
trans.vinc=norm(findvpqw(etmtransorb,pi/2));
vinfi(2,1)=trans.vt-mars.v;
energy=vinfi.^2./2;
pu=[planet.earth.u;planet.mars.u];
parkingr=[orbit.ep.a;planet.stm.radius+500];
parkingv=sqrt(pu./parkingr)
vp=sqrt((energy+pu./parkingr).*2)
dv=vp-parkingv
dv(3)=2*trans.vinc*sin(etmtransorb.inc/2)
%%
%phase=timebetween(orbit.ste,0,pi)

%%
t=1;
for i=0:800
    rstm(:,t)=orbvec3d(orbit.stm,t);
    rste(:,t)=orbvec3d(orbit.ste,t);

    t=t+1;
end

%%

%plot3(r(1,:),r(2,:),r(3,:),'b','LineWidth',2)
hold on
plot3(rste(1,:),rste(2,:),rste(3,:))
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3d trace')
grid on
axis equal
%%
%grid on;
%for i=0:366
%    h=plot3([],[],[],LineWidth=2);
%    set(h,'Xdata',r(1,:),'Ydata',r(2,:),'Zdata',r(3,:));
%    pause(0.05);
%end


%%
% 生成 3D 轨迹数据
t = linspace(0, 10, 200);
x = rstm(1,:);
y = rstm(2,:);
z = rstm(3,:);
x2=rste(1,:);
y2=rste(2,:);
z2=rste(3,:);

% 创建图形窗口
figure;
h = plot3(NaN, NaN, NaN, '-o', 'LineWidth', 2); % 先创建空轨迹
hold on
h2 = plot3(NaN, NaN, NaN, '-r', 'LineWidth', 2);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('动态 3D 轨迹');
grid on;

% 设置 X/Y/Z 轴比例
axis equal;
pbaspect([1 1 1]);
view(3);

% 动画循环
for i = 1:length(rstm)
    % 更新轨迹数据
    set(h, 'XData', x(1:i), 'YData', y(1:i), 'ZData', z(1:i));
    
    set(h2, 'XData', x2(1:i), 'YData', y2(1:i), 'ZData', z2(1:i));
    pause(0.02); % 控制更新速度
end
%%
reperi=orbvec3dano(orbit.ste,2)
ermperi=orbvec3dano(orbit.stm,0)
repjv=projectionvec(reperi,orbit.ste.plane)
%r1=[4;5;1]
%r1pjv=projectionvec(r1,orbit.ste.plane)
angle=rad2deg(findangle(reperi,repjv))


%%
%load('finaldata.mat');
%juliannow=juliandate(datetime(2025,3,10));




