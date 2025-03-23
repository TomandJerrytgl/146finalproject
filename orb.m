classdef orb
    properties
        a
        e
        u
        inc
        raan
        aop
        T
        p
        peritime
    end
    
    methods
        % 构造函数
        function obj = orb(a,e,u,inc,raan,aop,t)
            if nargin > 0 % 确保有输入参数
                obj.a = a;
                obj.e = e;
                obj.u = u;
                obj.inc=inc;
                obj.raan=raan;
                obj.aop=aop;
                obj.peritime=t;
                obj.T=(2*pi/sqrt(u))*a^(3/2)/3600/24;
                obj.p=a*(1-e^2);
            end
        end

        function trueanomaly=ttrueano(obj,timenow)
            t=timenow-obj.peritime;
           
            M=2.*pi.*t./obj.T;
            f = @(x) x - obj.e*sin(x)-M;
            E=fzero(f,1);
            trueanomaly=2.*atan(sqrt((1+obj.e)/(1-obj.e)).*tan(E/2));
        end


        function rvect=trvec(obj,t)
            ano=obj.ttrueano(t);
            rnormal=obj.p/(1+obj.e*cos(ano));
            rvect=[rnormal*cos(ano);rnormal*sin(ano);0];
        end

        function vvect=tvvec(obj,t)
            ano=obj.ttrueano(t);
            vvect=[-sqrt(obj.u/obj.p)*sin(ano);sqrt(obj.u/obj.p)*(obj.e+cos(ano));0];
        end

        function transmatrix=pqw2xyz(obj)
            Omega=obj.raan;
            i=obj.inc;
            w=obj.aop;


            R1 = [cos(Omega), -sin(Omega), 0;sin(Omega),  cos(Omega), 0;0,0,1];

            R2 = [1,  0,   0;
                 0, cos(i), -sin(i);
                   0, sin(i),  cos(i)];

            R3 = [cos(w), -sin(w), 0;
                sin(w),  cos(w), 0;
                0,      0,      1];
            transmatrix=R1*R2*R3;




        end




    end
end

