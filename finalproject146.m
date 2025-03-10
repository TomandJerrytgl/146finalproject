%engmae146 final
load('finaldata.mat');
orbit.stm.mass=planet.mars.mass;
orbit.ste.mass=planet.earth.mass;
orbit.stm.periapsis=orbitpostransfer(orbit.stm,orbitpos(findr(orbit.stm,0),0))
juliannow=juliandate(datetime(2025,3,10));

Marsanonow=caltruano(orbit.stm,juliannow);
Marsposnow=findr(orbit.stm,Marsanonow);
tratime=timebetween(orbit.stm,2.5,3.6);

angle=findangle([1,1,0],[0,0,1]);

%%

%%
%finding trueanomaly with orbit and time
function truanomaly=caltruano(orbit,timenow)
%orbit,timenow
e=orbit.e;
p=orbit.p;
a=orbit.a;
t=timenow-orbit.peritime;
M=2*pi*t/orbit.T;
f = @(x) x - e*sin(x)-M;
E=fzero(f,1);
truanomaly=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
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
theta=acos(norm(vector1.*vector2)./norm(vector1)./norm(vector2));
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
Omega=deg2rad(orbit.raan);
i=deg2rad(orbit.inc);
w=deg2rad(orbit.aop);


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
function launchwin=findnextlaunchwin(orbit1,orbit2,t,u,m)
i=0;
timenow=t;
mass1=orbit1.mass;
mass2=m;
orbit.trans.u=u;
phasetime=1000;
matrix1=[];
while i<1000&&abs(phasetime)>1
    anor1=caltruano(orbit1,timenow);
    r1=findr(orbit1,anor1);
    rsoi=gravifield(r1,mass1,mass2);
    ra=findr(orbit1,(anor1+pi));
    vecra=orbitpostransfer(orbit1,orbitpos(ra,anor1+pi));
    anor2=findangle(orbit2.periapsis,projectionvec(vecra,orbit2.plane));
    
    r2=findr(orbit2,anor2);
    orbit.trans.a=(r1+r2+rsoi)/2;
    orbit.trans.T=findT(orbit.trans);
    transfertime=orbit.trans.T/2/3600/24;
    marsnow=caltruano(orbit2,timenow);
   
    planetarrivaltime=timebetween(orbit2,marsnow,anor2);
    
    phasetime=transfertime-planetarrivaltime;
    timenow=timenow+1;

    


    i=i+1;
    matrix1(i,:)=[i,transfertime,planetarrivaltime,phasetime];

    
%    disp(orbit.trans)
%    disp([r1,r2,rsoi])
    
  

end
disp(matrix1)
launchwin=timenow;



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
nextluanchwin=findnextlaunchwin(orbit.ste,orbit.stm,juliannow,planet.sun.u,planet.sun.mass)
%%
%phase=timebetween(orbit.ste,0,pi)