classdef ButterflyFrame < handle
    %BUTTERFLYFRAME Butterfly frame class
    %   Contains dynamics and simple simulation software
    
    properties
        J = 1.8e-3/2;    % Moment of intertia of the frame
        F = 0.005      % Friction coefficient
        K = 1          % Motor coupling coefficient
        kp = 1;        % PD proportional gain
        kd = 1;        % PD derivative gain
        ks1 = 1;       % Sliding mode gain 1
        ks2 = 1;       % Sliding mode gain 2
        eps =1;        % Epsilon for sat(s/eps) for SMC type 2
        SMCt=1;        % Sliding mode type. 1: sign(s); 2: sat(s)
        Klqr = eye(1)  % LQR gain matrix
        ref = 1;       % Reference signal
        a = pi/3;      % Amplitude
        w = 1*pi;      % Frequency
        o = 0*pi/2;    % Offset
        cont= 1;       % Controller type
        fricComp = 0;  % Add friction compensation (true/false);
        ffwd = 0;      % Add feed forward of angular acceleration reference
    end    
    methods
        % It is possible to add a constructor should you so wish.
%         function obj=ButterflyFrame()
%         end
        function [r,Dr,DDr]  = getReference_r(obj,t,x)
            switch obj.ref
                case 1 
                    r  = pi/2;
                    Dr = 0;
                    DDr= 0;
                case 2
                    r  = obj.a*sin(obj.w*t+obj.o);
                    Dr = obj.a*obj.w*cos(obj.w*t+obj.o);
                    DDr=-obj.a*obj.w^2*sin(obj.w*t+obj.o);
            end
        end
        function u = getController_u(obj,t,x)
            theta  = x(1)-sign(x((1)))*floor(abs(x(1))/pi)*pi; % Same as the "excess" function
            Dtheta = x(2);
            % Get reference signal
            [r,Dr,DDr]= getReference_r(obj,t,x);
            % Compute errors
            e = r-theta;
            De= Dr-Dtheta;
            % Choose controller  
            switch obj.cont
                case 1 % PD controller                    
                    u = obj.kp*e+obj.kd*De; % Critically damped if kd=2*sqrt(kp)/sqrt(K/J)
                case 2 % Sliding mode control                    
                    s = e*obj.ks1+De;                    
                    switch obj.SMCt
                        case 1 % Sliding mode with sign
                            u = (obj.J/obj.K)*(obj.ks1*De)+(obj.F+obj.J*obj.a*obj.w^2+obj.ks2)*sign(s);                        
                        case 2 % Sliding mode with sat                            
                            s_sat=s/obj.eps;
                            if abs(s/obj.eps) > 1
                                s_sat=sign(s/obj.eps);
                            end
                            u = (obj.J/obj.K)*(obj.ks1*De)+(obj.F+obj.J*obj.a*obj.w^2+obj.ks2)*s_sat;
                    end
                case 3 %Linear quadratic regulator (LQR)                    
                    u = obj.Klqr*[e;De]; 
            end
            if obj.fricComp % Friction compensation
                u = u + (1/obj.K)*obj.F*sign(Dtheta);
            end
            if obj.ffwd     % Feed-forward of anuglar acceleration reference
                u = u + (obj.J/obj.K)*DDr;
            end
            % Saturation: abs(u)<=0.1
            if u > 0.1
                u = 0.1;
            elseif u < -0.1
                u = -0.1;
            end
        end
        function [t,x,u] = simEOM(obj,t0,tEnd,x0)
            options = odeset( 'RelTol', 1e-8, 'AbsTol', 1e-8);
            tint=linspace(t0,tEnd,round(tEnd)*100);
            [t,x] = ode23(@(t,y)EOM(obj,t,y),tint,x0,options);
            %[t,x] = ode45(@(t,y)EOM(obj,t,y),[t0,tEnd],x0,options);
            [u,~,~]=get_u_r_Dr(obj,t,x);
            
        end
        function dx = EOM(obj,t,x)
            %theta  = x(1);
            Dtheta = x(2);
            u      = getController_u(obj,t,x);
            DDtheta= get_DDtheta(obj,t,x,u);
            dx = [Dtheta;DDtheta];
        end
        function DDtheta = get_DDtheta(obj,t,x,u)
            %theta = x(1);
            Dtheta = x(2);           
            DDtheta= (1/obj.J)*(obj.K*u-obj.F*sign(Dtheta));
        end
        function [u,r,Dr] = get_u_r_Dr(obj,t,x)
            n = length(x);
            u = zeros(n,1);
            r = zeros(n,1);    
            Dr= zeros(n,1);
            for i=1:n                
                u(i)    = getController_u(obj,t(i),x(i,:));
                [r(i),Dr(i),~]= getReference_r(obj,t(i),x(i,:));
            end
        end
        function animate(obj,t,x,plotSpeed)
            % Shape of the frame is given by  \phi \in [0,2\pi],
            % \delta(\phi)=a-b\cos(2\phi)
            ap  = 0.1095;   
            bp  = 0.0405;    
            fig20 = figure(20);clf(20);
            grid on; 
            axis equal;
            hold on; 
            axis( [-0.2 0.2 -0.2 0.2]); 
            % The integrated data 
            N = length(x); 
            figure(fig20);             
            % The animation 
            for i=1:4:N
                q1  = x(i,1); 
                q1r =  obj.a*sin(obj.w*t(i)+obj.o);
                % Frame
                phi = 0:0.01:2*pi;
                % Actual position
                rf = ap-bp*cos(2*phi);
                xf = rf.*sin(phi-q1);
                yf = rf.*cos(phi-q1);
                % Desired position
                xfr = rf.*sin(phi-q1r);
                yfr = rf.*cos(phi-q1r);
                % Update figure
                if i == 1
                    frame = plot(xf,yf,'Color','black');
                    frame_ref = plot(xfr,yfr,'--r');
                    legend('Actual','Desired');
                    pause(1);
                else
                    set(frame,'XData',xf,'YData',yf);
                    set(frame_ref,'XData',xfr,'YData',yfr);
                end
                % Show the time
                %pause( t(end)/(length(t)/4)/plotSpeed ); 
                pause(0.03/plotSpeed);
                title(strcat('t=',num2str(t(i))));
            end
            hold off;
        end
    end
end