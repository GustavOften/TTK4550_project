classdef butterfly_robot_phi_varphi
    %Differential equation:
    %M(q)dot_dot_q + C(q, dot_q)dot_q + G(q) = [u;0]
    properties
        a = 0.1095; %Constant for describing the butterfly frame
        b = 0.0405; %Constant for describing the butterfly frame
        %c = 0.49; %Constant for virtual holonomic constraint.
        g = 9.81; %Gravity
        J_s = 5.48e-7; %Mass moment of inertia of the ball
        J_f = 1.581e-2; %Mass moment of inertia of the frame, from article.
        %J_f = 1.8e-3/2; %From lab in modrob.
        m_b = 3.0e-3; %Mass of the ball
        %m_f = 0.15; %Mass of the frame
        r_f = 16.5e-3; %Half of the distance between both frame plates
        R_b = sqrt(16.55e-3^2-12.5e-3^2); %Radius of the ball
        delta
        tau_delta
        rho
        length_rho
        bad_rho
        kappa_frame
        %alpha
        curve_for_phi        
        ds
        dds
        tau
        kappa
        p
        d_p
        R
        diff_R
        tau_diff 
        theta
        d_theta
        dd_theta
        K
        function_for_X
        function_for_dphi
        Gamma = 1
        Q = @(t)[1 0 0;0 1 0;0 0 1];
        X
        phi_test
        A
        B
        Phi
        dPhi
        ddPhi
        dddPhi
        normal_vector
    end
    methods
        function obj = butterfly_robot_phi_varphi(calculate_riccati)
            %% phi = f(varphi), varphi = g(varphi) 
            delta_curve_fit = @(phi)(obj.a - obj.b*cos(2*phi))*[sin(phi);cos(phi);0];
            tau_curve_fit = @(phi)((-2*obj.b*sin(2*phi))*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0])/ ...
                norm((-2*obj.b*sin(2*phi))*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]);
            
            range_for_functions = linspace(0,2*pi,500);
            function_g = @(x) atan2([1 0 0]*delta_curve_fit(x) - obj.R_b*[0 1 0]*tau_curve_fit(x), ...
                                [0 1 0]*delta_curve_fit(x)+obj.R_b*[1 0 0]*tau_curve_fit(x));
            res_fun_f = zeros(length(range_for_functions),1);
            for i = 1:length(range_for_functions)
                res_fun_g(i) = function_g(range_for_functions(i));
            end
            res_fun_g = unwrap(res_fun_g);
            figure
            plot(res_fun_g,range_for_functions)
            hold on;
            
            k = spline(res_fun_g, range_for_functions)
            dk = fnder(k, 1)
            ddk = fnder(k, 2)
            dddk = fnder(k, 3)
            result_spline = @(x)ppval(k,mod(x,2*pi));
            result_dspline = @(x)ppval(dk,mod(x,2*pi));
            result_ddspline = @(x)ppval(ddk,mod(x,2*pi));
            result_dddspline = @(x)ppval(dddk,mod(x,2*pi));

            plot(res_fun_g, result_spline(res_fun_g));
            plot(res_fun_g, result_dspline(res_fun_g));
            plot(res_fun_g, result_ddspline(res_fun_g));
            plot(res_fun_g, result_dddspline(res_fun_g));
            legend("true", "spline", "d_spline", "dd_spline", "ddd_spline");
            grid on;
            obj.Phi = result_spline;
            obj.dPhi = result_dspline;
            obj.ddPhi = result_ddspline;
            obj.dddPhi = result_dddspline;
            
            
            syms phi(varphi) theta Phi DPhi DDPhi DDDPhi real 
            
            dpdv = diff(phi(varphi),varphi);
            ddpdvv = diff(phi(varphi),varphi, varphi);
            dddpdvvv = diff(phi,varphi,varphi,varphi);
            %% Delta: parametrization of the shape of the frame.
            delta = (obj.a - obj.b*cos(2*phi))*[sin(phi);cos(phi);0];
            obj.delta = matlabFunction(subs(delta,phi,Phi), 'Vars',[Phi;varphi]);
            l_delta = obj.a - obj.b*cos(2*phi);
            l_d_delta = 2*obj.b*sin(2*phi)*diff(phi,varphi);
            
            %% Tau: 
            tau = (2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]) ...
                /sqrt(sum((2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2));
            obj.tau = matlabFunction(subs(tau, phi,Phi),'Vars',[Phi;varphi]);
            
            %% Kappa frame:
            kappa_frame = norm((4*obj.b*cos(2*phi)*[sin(phi);cos(phi);0]+2*obj.b*sin(2*phi)*[cos(phi);-sin(phi);0] + ...
                2*obj.b*sin(2*phi)*[cos(phi);-sin(phi);0]+(obj.a - obj.b*cos(2*phi))*[-sin(phi);-cos(phi);0])/ ...
                sum((2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2));
            obj.kappa_frame = matlabFunction(subs(kappa_frame,phi,Phi),'Vars',[Phi;varphi]);
            
            %% varphi = g(phi) 
            g = subs(atan2([1 0 0]*delta-[0 1 0]*tau,[0 1 0]*delta+[1 0 0]*tau),phi,Phi);
            dg = diff(g,Phi);
            ddg = diff(dg,Phi);
            
            %% p and dp 
            p = 1/(obj.R_b*(-1+obj.R_b*kappa_frame));
            d_p = diff(p, varphi);
            obj.p = matlabFunction(subs(p, phi,Phi),'Vars',[Phi;varphi]);
            obj.d_p = matlabFunction(subs(d_p, {phi,dpdv},{Phi,DPhi}),'Vars',[Phi;DPhi;varphi]);
            
            %% Rho: Vector to center of ball, given in bodyframe.
            bad_rho = (obj.a-obj.b*cos(2*phi)+obj.R_b)*[sin(phi);cos(phi);0];
            obj.bad_rho = matlabFunction(bad_rho);
            rho = delta + obj.R_b*[[0 -1 0]*tau;[1 0 0]*tau;0];
            obj.normal_vector = matlabFunction(subs([[0 -1 0]*tau;[1 0 0]*tau;0],phi,Phi),'Vars',[Phi;varphi]);
            obj.rho = matlabFunction(subs(rho, phi,Phi),'Vars',[Phi;varphi]);
            obj.length_rho = matlabFunction(subs(norm(rho),phi,Phi),'Vars',[Phi;varphi]);
            
            %% s: Arclength of the balls path.
            ds = norm(diff(rho,varphi));
            dds = diff(ds, varphi);
            obj.ds = matlabFunction(subs(ds, {phi,dpdv,ddpdvv},{Phi,DPhi,DDPhi}),'Vars',[Phi;DPhi;DDPhi;varphi]);
            obj.dds = matlabFunction(subs(dds, {phi,dpdv,ddpdvv},{Phi,DPhi,DDPhi}),'Vars',[Phi;DPhi;DDPhi;varphi]);

            %% Kappa: Curvature of curve.
            kappa = diff(rho,varphi,varphi)/sum(diff(rho,varphi).^2);
            obj.kappa = matlabFunction(subs(kappa, {phi,dpdv,ddpdvv, dddpdvvv},{Phi,DPhi,DDPhi,DDDPhi}),'Vars',[Phi;DPhi;DDPhi;DDDPhi;varphi]);

            %% Rotational matrices:
            obj.R = matlabFunction([cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]);
            obj.diff_R = matlabFunction([-sin(theta) -cos(theta) 0;cos(theta) -sin(theta) 0;0 0 0]);
            
            %% Control stuff
            Rtau = obj.R(varphi)*obj.tau(varphi)
            %% Theta used in Case-study non-prehensile %%%%%%%
            a = -0.03;
            numerator = a*sin(2*varphi-pi)*[1 0 0]*Rtau-[0 1 0]*Rtau;
            denominator = a*sin(2*varphi-pi)*[0 1 0]*Rtau+[1 0 0]*Rtau;
            %Theta = varphi - 0.49*sin(2*varphi)
            Theta = varphi+atan(numerator/denominator)
            dTheta = diff(Theta,varphi)
            ddTheta = diff(dTheta,varphi)
            %% Theta used in Internship report
            %Theta = phi-1.3*sin(2*phi);
            %Theta = phi-0.2*sin(Rtau(1));
            obj.theta = matlabFunction(subs(Theta, phi,Phi),'Vars',[Phi;varphi]);
            obj.d_theta = matlabFunction(subs(dTheta, {phi,dpdv},{Phi,DPhi}),'Vars',[Phi;DPhi;varphi]);
            obj.dd_theta = matlabFunction(subs(ddTheta, {phi,dpdv,ddpdvv},{Phi,DPhi,DDPhi}),'Vars',[Phi;DPhi;DDPhi;varphi]);

            %% Plots phase plane of alpha betta gamma function
             a = @(x) obj.alpha_beta_gamma(x);
             f_plane = @(t,x) [x(2);-1/(a(x(1))'*[1;0;0])*((a(x(1))'*[0;1;0])*x(2)^2+a(x(1))'*[0;0;1])];
%             phase_plot_2_interactive(f_plane,[0 6*pi;-1 10],10,'',[100,100],0.1)
            %% Finding solution to periodic Riccati equation
            if calculate_riccati
                options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
                [ts,ys] = ode45(f_plane,[0,10],[0;1], options);
                i = 1;
                phi_dot = [0 0];
                length(ys)
                while ys(i,1) <= 2*pi
                    phi_dot(i,:) = ys(i,:);
                    i = i+1;
                    if i >= length(ys)
                        break
                    end
                end
                phi_dot(i,:) = ys(i,:)
                %% Interpolated 
                interpolation_for_dphi = @(x) interp1(phi_dot(:,1),phi_dot(:,2), x);
                obj.function_for_dphi = @(phi) interpolation_for_dphi(mod(phi,pi));
                % Using the curvfitted function since ode45 has problems with
                % interpolations.
                [A, B] = obj.get_linearization(phi, obj.function_for_dphi,true)          
                obj.A = @(phi) A(phi)
                obj.B = @(phi) B(phi)
                A(0)
                B(0)
                [X,phi] = sdp_riccati(obj.A,obj.B,obj.Q,obj.Gamma,0,pi,700,7,3);
                obj.X = X; obj.phi_test = phi;
                %% Interpolation for Riccati solution
                interpolation_for_X = @(x)[interp1(phi,reshape(X(1,1,:),1,length(X)),x) interp1(phi,reshape(X(1,2,:),1,length(X)),x) interp1(phi,reshape(X(1,3,:),1,length(X)),x);
                                          interp1(phi,reshape(X(2,1,:),1,length(X)),x) interp1(phi,reshape(X(2,2,:),1,length(X)),x) interp1(phi,reshape(X(2,3,:),1,length(X)),x);
                                          interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x);]
%                    iX = cell(length(X), length(X));
%                 for i = 1:length(X)
%                     for j = 1:length(X)
%                         iX{i,j} = @(x) interp1(phi,reshape(X(i,j,:),1,length(X)),x);
%                     end
%                 end                      
%                 obj.function_for_X = @(x) [iX{1,1}(x) iX{1,2}(x) iX{1,3}(x) iX{1,4}(x) iX{1,5}(x) iX{1,6}(x);
%                                            iX{2,1}(x) iX{2,2}(x) iX{2,3}(x) iX{2,4}(x) iX{2,5}(x) iX{2,6}(x);
%                                            iX{3,1}(x) iX{3,2}(x) iX{3,3}(x) iX{3,4}(x) iX{3,5}(x) iX{3,6}(x);
%                                            iX{4,1}(x) iX{4,2}(x) iX{4,3}(x) iX{4,4}(x) iX{4,5}(x) iX{4,6}(x);
%                                            iX{5,1}(x) iX{5,2}(x) iX{5,3}(x) iX{5,4}(x) iX{5,5}(x) iX{5,6}(x);
%                                            iX{6,1}(x) iX{6,2}(x) iX{6,3}(x) iX{6,4}(x) iX{6,5}(x) iX{6,6}(x)];
                obj.function_for_X = @(x) interpolation_for_X(mod(x,pi));
            end
        end
        
        function delta = get_delta(obj, varphi)
            delta = obj.delta(obj.Phi(varphi),varphi);
        end
        function tau = get_tau(obj,varphi)
            tau = obj.tau(obj.Phi(varphi),obj.dPhi(varphi),varphi);
        end
        function normal_vector = get_normal_vector(obj,varphi)
            normal_vector = obj.normal_vector(obj.Phi(varphi),obj.dPhi(varphi),varphi);
        end
        function rho = get_rho(obj, varphi)
           rho = obj.rho(obj.Phi(varphi),obj.dPhi(varphi), varphi);
        end
        function l_rho = get_l_rho(obj, varphi)
            l_rho = obj.length_rho(obj.Phi(varphi), obj.dPhi(varphi),varphi);
        end
        function kappa = get_kappa(obj, varphi)
            kappa = obj.kappa(obj.Phi(varphi), obj.dPhi(varphi), obj.ddPhi(varphi), ...
                obj.dddPhi(varphi), varphi);
        end
        function ds = get_ds(obj, varphi)
            ds = obj.ds(obj.Phi(varphi), obj.dPhi(varphi), obj.ddPhi(varphi),varphi);
        end
        function dds = get_dds(obj, varphi)
            dds = obj.dds(obj.Phi(varphi), obj.dPhi(varphi), obj.ddPhi(varphi),varphi);
        end
        function p = get_p(obj, varphi)
            p = obj.p(obj.Phi(varphi), varphi);
        end
        function dp = get_dp(obj, varphi)
            dp = obj.d_p(obj.Phi(varphi), obj.dPhi(varphi), varphi);
        end
        function theta = get_theta(obj, varphi)
            theta = obj.theta(obj.Phi(varphi),varphi);
        end
        function d_theta = get_dtheta(obj, varphi)
            d_theta = obj.d_theta(obj.Phi(varphi),obj.dPhi(varphi),varphi);
        end
        function dd_theta = get_ddtheta(obj, varphi)
            dd_theta = obj.dd_theta(obj.Phi(varphi), obj.dPhi(varphi), obj.ddPhi(varphi),varphi);
        end
        
        function M = get_M(obj, q)
            varphi = q(2);
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi), obj.get_tau(varphi));
            m11 = obj.get_l_rho(varphi)^2+obj.J_f/obj.m_b+obj.J_s/obj.m_b;
            m12 = obj.get_ds(varphi)*(rhoxtau+obj.J_s/obj.m_b*obj.get_p(varphi));
            m22 = (1+obj.J_s/obj.m_b*obj.get_p(varphi)^2)*obj.get_ds(varphi)^2;
            M =   obj.m_b*[m11 m12;
                   m12 m22];
        end
        
        function C = get_C(obj, q, dq)
            varphi = q(2);
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi), obj.get_tau(varphi));
            taudotrho = obj.get_tau(varphi)'*obj.get_rho(varphi);
            rhoxkappa = [0 0 1]*cross(obj.get_rho(varphi),obj.get_kappa(varphi));
            c11 = obj.m_b*obj.get_ds(varphi)*taudotrho*dq(2);
            
            c12 = obj.m_b*(taudotrho*obj.get_ds(varphi)*dq(1) ...
                +(obj.get_dds(varphi)*(rhoxtau+obj.J_s/obj.m_b*obj.get_p(varphi)) ...
                +obj.get_ds(varphi)^2*rhoxkappa ... 
                +obj.get_ds(varphi)*obj.J_s/obj.m_b*obj.get_dp(varphi))*dq(2));
            
            c21 = -obj.m_b*obj.get_ds(varphi)*taudotrho*dq(1);
            
            c22 = obj.m_b*(1+obj.J_s/obj.m_b*obj.get_p(varphi)^2)*obj.get_ds(varphi)*obj.get_dds(varphi)*dq(2) + ...
                obj.get_ds(varphi)*obj.J_s*obj.get_p(varphi)*obj.get_dp(varphi);
            
            C = [c11 c12;
                 c21 c22];
        end
        
        function G = get_G(obj, q)
            varphi = q(2);
            G = [obj.m_b*[0;obj.g;0]'*obj.diff_R(q(1))*obj.get_rho(varphi);
                 obj.m_b*[0;obj.g;0]'*obj.R(q(1))*obj.get_tau(varphi)*obj.get_ds(varphi)];
        end
        
        function ddq = calculate_ddq(obj, q, dq, u)
            C = obj.get_C(q,dq);
            G = obj.get_G(q);
            M = obj.get_M(q);
            ddq = (M^-1*(-C*dq-G+[u;0]));
        end
  
        function plot_constraints(obj)
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            figure 
            k = linspace(0,pi,100);
            constraint_eq_1 = obj.get_ds(0)*(1+obj.J_s/obj.m_b*obj.p(0)^2)/([0 0 1]*cross(obj.rho(0),obj.get_tau(0))+obj.J_s/obj.m_b*obj.get_p(0));
            constraint_eq_2 = obj.get_ds(pi/2)*(1+obj.J_s/obj.m_b*obj.get_p(pi/2)^2)/([0 0 1]*cross(obj.get_rho(pi/2),obj.get_tau(pi/2))+obj.J_s/obj.m_b*obj.get_p(pi/2));
            not_eq_constraint = @(x) -obj.get_ds(x)*(1+obj.J_s/obj.m_b*obj.get_p(x)^2)/([0 0 1]*cross(obj.get_rho(x),obj.get_tau(x))+obj.J_s/obj.m_b*obj.get_p(x));
            %constraint_eq_1 = obj.ds(0)*(1+obj.J_s/obj.m_b/obj.R_b^2)/([0 0 1]*cross(obj.rho(0),obj.tau(0))-obj.J_s/obj.m_b/obj.R_b);
            %constraint_eq_2 = obj.ds(pi/2)*(1+obj.J_s/obj.m_b/obj.R_b^2)/([0 0 1]*cross(obj.rho(pi/2),obj.tau(pi/2))+obj.J_s/obj.m_b/obj.R_b);
            %not_eq_constraint = @(x) -obj.ds(x)*(1+obj.J_s/obj.m_b/obj.R_b^2)/([0 0 1]*cross(obj.rho(x),obj.tau(x))-obj.J_s/obj.m_b/obj.R_b);
            hold on;
            results = zeros(length(k),1);
            parfor i = 1:length(k)
                results(i,1) = not_eq_constraint(k(i));
            end
            plot(k,results);
            scatter([0;pi],[constraint_eq_1;constraint_eq_1],'x','red');
            scatter(pi/2, constraint_eq_2, 'x','blue');
            grid on;
            legend("Constraint all $\varphi$", "Constraint eq. point  $\varphi = n\pi$", "Constraint eq.point $\varphi = n\frac{\pi}{2}$");
            title("Constraints for $\Theta'(\varphi)$")
            set(gca,'XTick',0:pi/4:pi) 
            set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$'})
        end

        function normal_force = get_normal_force(obj, q)
            %Using a simplified model. 
            % Assuming the normal vector is 90 deg anti-clockwise to tau.
            n = [0 -1 0;1 0 0;0 0 1]*obj.get_tau(q(2));
            normal_force = [0;9.81;0]'*obj.R(q(1))*n;
        end
        function position_ball = get_position_ball(obj, q)
            ball_position = obj.get_rho(q(2));
            position_ball = obj.R(q(1))*ball_position;
        end
        function velocity_ball = get_velocity_ball(obj, q, dq)
            txRrho = cross([0;0;dq(1)],obj.R(q(1))*obj.get_rho(q(2)));
            velocity_ball = txRrho + obj.R(q(1))*obj.get_tau(q(2))*obj.get_ds(q(2))*dq(2);
        end
        function M = make_movie_butterfly_robot(obj, results, fps)
            h = figure;
            time = results.tout(end);
            number_of_frames = round(time*fps);
            x_frame = zeros(1,length(results.tout));
            y_frame = zeros(1,length(results.tout));
            x_ball = obj.R_b*cos(linspace(0, 2*pi,20));
            y_ball = obj.R_b*sin(linspace(0, 2*pi,20));
            xlim([-0.2 0.2])
            ylim([-0.2 0.2])
            axis equal;
            grid on;
            M(number_of_frames) = struct('cdata',[],'colormap',[]);
            h.Visible = 'off';
            j = 0;
            for i = linspace(0,2*pi)
                j = j + 1;
                position = obj.delta(i);
                x_frame(j) = position(1);
                y_frame(j)= position(2);
            end
            frame = hgtransform;
            ball = hgtransform;
            patch('XData',x_frame,'YData',y_frame,'FaceColor','yellow','Parent',frame) 
            patch('XData',x_ball,'YData',y_ball,'FaceColor','red','Parent',ball) 
            start_falling = length(results.after_fall(:,1))-length(find(results.after_fall(:,1)));
            acumulated_time = 0;
            current_frame = 0;
            for i = 2:start_falling
                acumulated_time = acumulated_time + results.tout(i) - results.tout(i-1);
                if acumulated_time >= 1/fps
                    ball_position = obj.rho(results.q(i,2));
                    ball_position_inertial = obj.R(results.q(i,1))*ball_position;
                    frame.Matrix = makehgtform('zrotate',results.q(i,1));
                    ball.Matrix = makehgtform('translate', ball_position_inertial);
                    drawnow
                    current_frame = current_frame + 1
                    M(current_frame) = getframe;
                    acumulated_time = acumulated_time - 1/fps;
                end
            end
            for i = start_falling+1:length(results.after_fall(:,1))
               acumulated_time = acumulated_time + results.tout(i)-results.tout(i-1);
               if acumulated_time >= 1/fps
                    frame.Matrix = makehgtform('zrotate',results.after_fall(i,1));
                    ball.Matrix = makehgtform('translate', [results.after_fall(i,2);results.after_fall(i,3);0]);
                    drawnow
                    current_frame = current_frame + 1
                    M(current_frame) = getframe;
                    acumulated_time = acumulated_time - 1/fps;
                end
            end
            h.Visible = 'on';
            movie(M,1,fps);
        end
        
        
        
        %% Controller Stuff
        
        function a = alpha_beta_gamma(obj, varphi)
            %alpha = [obj.get_dtheta(varphi) 1]*obj.get_M([obj.get_theta(varphi);varphi])*[0;1];
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            rhodottau = obj.get_rho(varphi)'*obj.get_tau(varphi);
            alpha = obj.m_b*obj.get_ds(varphi)*((rhoxtau + obj.J_s*obj.get_p(varphi)/obj.m_b)*obj.get_dtheta(varphi)+...
                obj.get_ds(varphi)*(1+obj.J_s/obj.m_b*obj.get_p(varphi)^2));
            
            %beta = [0 1]*obj.get_M([obj.get_theta(varphi);varphi])*[obj.get_ddtheta(varphi);0] + ...
                %[0 1]*obj.get_C([obj.get_theta(varphi);varphi],[obj.get_dtheta(varphi);1])*[1;1];
            beta = obj.m_b*obj.get_ds(varphi)*((rhoxtau+obj.J_s*obj.get_p(varphi)/obj.m_b)*obj.get_ddtheta(varphi)...
                -rhodottau*obj.get_dtheta(varphi)^2 + obj.get_dds(varphi)*(1+obj.J_s/obj.m_b*obj.get_p(varphi)^2) + ...
                obj.J_s/obj.m_b*obj.get_p(varphi)*obj.get_dp(varphi));
            gamma = obj.m_b*[0 obj.g 0]*obj.R(obj.get_theta(varphi))*obj.get_tau(varphi)*obj.get_ds(varphi);
            
            a = [alpha;beta;gamma];
        end
        function g = get_g_w_y_on_trajectory(obj, varphi, d_varphi)
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            g_w = -obj.m_b*obj.get_ds(varphi)*(rhoxtau+obj.J_s/obj.m_b*obj.get_p(varphi));
            g_y_dot = obj.m_b*obj.get_ds(varphi)*obj.tau(varphi)'*obj.rho(varphi)*(2*obj.get_dtheta(varphi)*d_varphi);
            g_y = -obj.get_ds(varphi)*obj.m_b*obj.g*[cos(obj.get_theta(varphi)) -sin(obj.get_theta(varphi)) 0]*obj.get_tau(varphi);
            g = [g_w;g_y_dot;g_y];
        end
        function [A, B] = get_linearization(obj, phi, phi_dot, function_handles)
            if function_handles == true
                abg = @(x)obj.alpha_beta_gamma(x);
                gwy = @(x)obj.get_g_w_y_on_trajectory(x, phi_dot(x));
                A = @(x)[0 1 0;
                         0 0 0;
                        [0 0 1]*gwy(x)/([1 0 0]*abg(x)) [0 1 0]*gwy(x)/([1 0 0]*abg(x)) ([0 0 1]*abg(x)-[0 1 0]*abg(x)*phi_dot(x)^2)/([1 0 0]*abg(x))/phi_dot(x)];
                B = @(x)[0;1;[1 0 0]*gwy(x)/([1 0 0]*abg(x))];
            else
                abg = obj.alpha_beta_gamma(phi);
                gwy = obj.get_g_w_y_on_trajectory(phi, phi_dot);
                A = [0 1 0;
                     0 0 0;
                     gwy(3)/abg(1) gwy(2)/abg(1) (abg(3)-abg(2)*phi_dot^2)/abg(1)/phi_dot];
                B = [0;1;gwy(1)/abg(1)];
            end
        end
        function w = get_w(obj, riccati_sol, q, epsilon)
            w = -(1/obj.Gamma)*obj.B(q(2))'*riccati_sol*epsilon;
        end
        function out = riccati_times_epsilon(obj, q, epsilon)
            out = obj.function_for_X(mod(q(2),2*pi))*epsilon;
        end
        function u = get_u(obj, q, dq, epsilon)
            L = [1 obj.get_dtheta(q(2));
                 0 1];
            N = [obj.get_ddtheta(q(2))*dq(2)^2;
                 0];
            M = obj.get_M([obj.get_theta(q(2));q(2)]);
            C = obj.get_C([obj.get_theta(q(2));q(2)],[obj.get_dtheta(q(2))*dq(2);dq(2)]);
            G = obj.get_G([obj.get_theta(q(2));q(2)]);
            denom = (L^-1)*(M^-1);
            numerator = L^-1*(N+(M^-1)*C*L*[epsilon(2);dq(2)]+(M^-1)*G)*1;
            numerator = numerator(1)+obj.get_w(obj.function_for_X(q(2)),q, epsilon);
            u = numerator/denom(1,1);
        end
        function u = get_my_u(obj, q, dq, epsilon)
            M = obj.get_M([obj.get_theta(q(2));q(2)]);
            C = obj.get_C([obj.get_theta(q(2));q(2)],[obj.get_dtheta(q(2))*dq(2);dq(2)]);
            G = obj.get_G([obj.get_theta(q(2));q(2)]);
            k = [1 0]*M*[obj.get_dtheta(q(2));1]/([obj.get_dtheta(q(2)) 1]*M*[0;1]);
            u = k*([k^-1 -1]*M*[obj.get_w(obj.function_for_X(q(2)),q, epsilon)+obj.get_ddtheta(q(2))*dq(2)^2;0]+[k^-1 -1]*C*[(obj.get_dtheta(q(2)))*dq(2)+epsilon(2);dq(2)]+[k^-1 -1]*G);
        end
        function epsilon = get_epsilon(obj, q, dq)
            epsilon = [q(1)-obj.get_theta(q(2));dq(1)-obj.get_dtheta(q(2))*dq(2);dq(2)-obj.function_for_dphi(q(2))];
        end
    end
end