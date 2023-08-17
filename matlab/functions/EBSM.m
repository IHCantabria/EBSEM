classdef EBSM
    properties
        time
        Hs
        Hb
        depthb
        SS
        AT
        Tp
        d50
        Omega
        MD
        Y09
        SF
    end
    methods
        function obj=init(obj, time, Hs, SS, AT, Tp, d50)
            obj.time=time;
            obj.Hs=Hs;
            obj.SS=SS;
            obj.AT=AT;
            obj.Tp=Tp;
            obj.d50=d50;
            
        end
        function obj=LinearBreak(obj, theta, depth, angleBathy)
            [obj.Hb, ~, obj.depthb] = BreakingPropagation(obj.Hs,obj.Tp,theta,depth,angleBathy);
            obj.depthb(obj.Hb==0) = 0.5/0.78;
            obj.Hb(obj.Hb==0) = 0.5;
            obj.Omega = obj.Hb ./ (wMOORE(obj.d50) .* obj.Tp);
        end
        function obj=MillerDean(obj, dt, Yi, Hberm, flagP)
            obj.MD.dt = dt;
            obj.MD.Yi = Yi;
            obj.MD.Hberm = Hberm;
            obj.MD.flagP = flagP;
            obj.MD.sl = obj.AT + obj.SS;
        end
        function obj=Yates09(obj, dt, S0)
            obj.Y09.dt = dt;
            obj.Y09.S0 = S0;
            obj.Y09.E = obj.Hb .^2;
        end
        function obj=ShoreFor(obj, dt, Sini)
            obj.SF.dt = dt;
            obj.SF.Sini = Sini;
            obj.SF.P = 1 ./ 16 .* 1025 .* 9.81 .* obj.Hb.^2 .* (9.81 .* obj.Hb./.78).^.5;
        end
    end
end