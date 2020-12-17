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
        dsf
        ddsf
        tau
        kappa
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
        Q = @(t)[20 0 0;0 20 0;0 0 20];
        X
        phi_test
        A
        B
        Phi
        dPhi
        ddPhi
        dddPhi
        normal_vector
        g_fun
        dg_fun
        ddg_fun
        drho_test 
        ys
    end
    methods
        function obj = butterfly_robot_phi_varphi(calculate_riccati)
            %% phi = f(varphi), varphi = g(varphi) 
            delta_curve_fit = @(phi)(obj.a - obj.b*cos(2*phi))*[sin(phi);cos(phi);0];
            tau_curve_fit = @(phi)(2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]) ...
                /sqrt(sum((2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2));
            
            range_for_functions = linspace(0,2*pi,500);
            function_g = @(x) atan2([1 0 0]*delta_curve_fit(x) - obj.R_b*[0 1 0]*tau_curve_fit(x), ...
                                [0 1 0]*delta_curve_fit(x)+obj.R_b*[1 0 0]*tau_curve_fit(x));
            res_fun_g = zeros(length(range_for_functions),1);
            res_fun_dg = zeros(length(range_for_functions),1);
            for i = 1:length(range_for_functions)
                res_fun_g(i) = function_g(range_for_functions(i));
            end
            res_fun_g = unwrap(res_fun_g);
            figure
            plot(res_fun_g,range_for_functions)
            hold on;
            
            k = spline(res_fun_g,range_for_functions);
            result_spline = @(x)ppval(k,mod(x,2*pi));
            subplot(1,2,1)
            plot(res_fun_g,result_spline(res_fun_g));
            
            legend("true", "spline");
            grid on;
            obj.Phi = result_spline;
          
            syms phi theta varphi
            
            %% Delta: parametrization of the shape of the frame.
            delta = (obj.a - obj.b*cos(2*phi))*[sin(phi);cos(phi);0];
            obj.delta = matlabFunction(delta);
            
            %% Tau: 
            tau = (2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]) ...
                /sqrt(sum((2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2));
            %tau = diff(delta,phi)/norm(diff(delta,phi));
            obj.tau = matlabFunction(tau);
            
            %% Kappa frame:
            kappa_frame = (4*obj.b*cos(2*phi)*[sin(phi);cos(phi);0]+2*obj.b*sin(2*phi)*[cos(phi);-sin(phi);0] + ...
                2*obj.b*sin(2*phi)*[cos(phi);-sin(phi);0]+(obj.a - obj.b*cos(2*phi))*[-sin(phi);-cos(phi);0])/ ...
                sum((2*obj.b*sin(2*phi)*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2);
            %kappa_frame = diff(tau,phi)/norm(diff(delta,phi));
            obj.kappa_frame = matlabFunction(kappa_frame);
            
            %% varphi = g(phi) 
            g = atan2([1 0 0]*delta-obj.R_b*[0 1 0]*tau,[0 1 0]*delta+obj.R_b*[1 0 0]*tau);
            dg = diff(g,phi);
            ddg = diff(dg,phi);
            obj.g_fun = matlabFunction(g);
            obj.dg_fun = matlabFunction(dg);
            obj.ddg_fun = matlabFunction(ddg);

            %% Rho: Vector to center of ball, given in bodyframe.
            rho = delta + obj.R_b*[[0 -1 0]*tau;[1 0 0]*tau;0];
            obj.normal_vector = matlabFunction([[0 -1 0]*tau;[1 0 0]*tau;0]);
            obj.rho = matlabFunction(rho);
            obj.length_rho = matlabFunction(norm(rho));
            
            %% s: Arclength of the balls path.
            ds = norm(diff(rho,phi))/dg;
            dds = diff(rho,phi)'*diff(rho,phi,phi)/(norm(diff(rho,phi))*dg^2)-norm(diff(rho,phi))*ddg/dg^3;
            obj.ds = matlabFunction(ds);
            obj.dds = matlabFunction(dds);
            
            %% sf: Arclength of the frame with respect to \varphi
            dsf = norm(diff(delta,phi))/dg;
            ddsf = diff(delta,phi)'*diff(delta,phi,phi)/norm(diff(delta,phi))/dg^2-norm(diff(delta,phi))*ddg/dg^3;
            obj.dsf = matlabFunction(dsf);
            obj.ddsf = matlabFunction(ddsf);
            
            %% Kappa: Curvature of curve.
            kappa = diff(rho,phi,phi)/dg^2/ds^2;%kappa_frame/(1-obj.R_b*norm(kappa_frame));
            obj.kappa = matlabFunction(kappa);

            %% Rotational matrices:
            obj.R = matlabFunction([cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]);
            obj.diff_R = matlabFunction([-sin(theta) -cos(theta) 0;cos(theta) -sin(theta) 0;0 0 0]);
            
            %% Control stuff
            %New variables to differentiate Theta(\varphi, \phi(\varphi)) with respect to \varphi
            syms varphi
            
            %% Theta used in Case-study non-prehensile %%%%%%%
            
%             R = [cos(varphi) -sin(varphi) 0;sin(varphi) cos(varphi) 0;0 0 1];
%             tau_theta = ((-2*obj.b*sin(2*varphi))*[sin(varphi);cos(varphi);0]+(obj.a - obj.b*cos(2*varphi))*[cos(varphi);-sin(varphi);0])/ ...
%                 norm((-2*obj.b*sin(2*varphi))*[sin(varphi);cos(varphi);0]+(obj.a - obj.b*cos(2*varphi))*[cos(varphi);-sin(varphi);0]);
%             Rtau = R*tau_theta;
%             a = 0.03;
%             numerator = a*sin(2*varphi)*[1 0 0]*Rtau -[0 1 0]*Rtau;
%             denominator = a*sin(2*varphi)*[0 1 0]*Rtau +[1 0 0]*Rtau;
%             Theta = varphi+atan(numerator/denominator);
%             dTheta = diff(Theta,varphi);
%             ddTheta = diff(dTheta,varphi);
            %% Denne er ganske god. IKke fjern den!
            %% works with phi_dot_zero = 2.6
            %a= 3.646053e+00; b= -1.373081e-01; c= -3.280457e-02; d= -5.262885e-03;
            %% ##########################################
            %% Also good.
            %a= -3.219353e+00; b= 1.838894e-01; c= 4.369538e-02; d= 3.738777e-03;
            %%
            %% PRETTY GOOD, center in pi/2.
            a= -8.696967e-01; b= 8.285338e-01; c= 2.104423e-01; d= 7.803824e-02; e= 4.994753e-02; f= 2.044127e-02;
            %%%
            %a= -1.058798e+01; b= 5.721937e-02; c= 1.156266e-02; d= 4.636083e-03; e= 1.862895e-03; f= 7.562149e-04;
            %a= -7.479424e-01; b= 9.567741e-01; c= 3.023719e-01; d= 2.335561e-02; e= 2.090285e-02; f= 1.754860e-02;
            %% Pretty good. varphi_dot_zero = 0.89;
            a= -8.017198e-01; b= 8.326337e-01; c= 2.178205e-01; d= 9.341089e-02; e= 2.924495e-02; f= -4.505087e-03;
            
            
            a= -6.362650e-01; b= 1.283998e+00; c= 3.471409e-01; d= 5.587734e-02; e= 4.910347e-03;
            a= -6.361089e-01; b= 1.284578e+00; c= 3.473617e-01; d= 5.586551e-02; e= 4.882655e-03;
            %% Denne kan gå veldig sakte. varphi_dot_zero = 0.224;
            a= -6.012854e-01; b= 1.415366e+00; c= 3.964811e-01; d= 5.645448e-02; e= 2.783807e-03;

            %% This was used for fast moving perpetual, dot_varphi = 1.17;
            a= -6.023940e+01; b= 9.115470e-03; c= 1.941337e-03; d= 1.731620e-03; e= 7.253528e-04;
            
            %% This was used for center in pi/2 a= -2.816689e+00; b= 1.850450e-01; c= -9.372952e-03; d= -3.295735e-02; e= -1.621506e-02;
            
            Theta = a*atan(b*sin(2*varphi)+c*sin(4*varphi)+d*sin(6*varphi)+e*sin(8*varphi))+varphi;
            %a=-3.201623e-01;
            %Theta = a*sin(2*varphi)+varphi;
            %a= -6.073674e-01; b= -1.083216e-01; c= -5.987843e-02; d= -5.190061e-02; e= 2.158077e-02;
            %Theta = a*sin(2*varphi)+b*sin(4*varphi)+c*sin(6*varphi)+d*sin(8*varphi)+e*sin(10*varphi)+varphi;
            %Theta = varphi -0.49*sin(2*varphi);
            %Theta = varphi -1.5076*sin(2*varphi)-0.0063*sin(4*varphi);
            dTheta = diff(Theta,varphi);%-0.2*cos(2*varphi)/(sin(2*varphi)^2+1)+1;
            ddTheta = diff(Theta,varphi,varphi);%0.4*cos(2*varphi)/(cos(2*varphi)^2+1)+0.8*sin(2*varphi)^2*cos(2*varphi)/(cos(2*varphi)^2+1)^2;
            %% Theta used in Internship report
            %Theta = phi-1.3*sin(2*phi);
            %Theta = phi-0.2*sin(Rtau(1));
            obj.theta = matlabFunction(Theta);
            obj.d_theta = matlabFunction(dTheta);
            obj.dd_theta = matlabFunction(ddTheta);
% 
%             %% Plots phase plane of alpha betta gamma function
            a = @(x) obj.alpha_beta_gamma(x);
            figure
            k = linspace(0,pi, 400);
            linear_term = @(x) -obj.get_dgamma(x(1))/([1 0 0]*a(x(1)));
            res = zeros(400,1);
            for i = 1:400
               res(i) = linear_term(k(i));%[0 0 1]*a(k(i))/([1 0 0]*a(k(i))); 
            end
            plot(k,res);
            f_plane = @(t,x) [x(2);-1/(a(x(1))'*[1;0;0])*((a(x(1))'*[0;1;0])*x(2)^2+a(x(1))'*[0;0;1])];
            %f_plane = @(t,x) [x(2);-obj.get_dgamma(0)/([1 0 0]*a(0))*x(1)];
            %phase_plot_2_interactive(f_plane,[0 pi;-2 2],10,'',[10,10],0.1)
            %% Finding solution to periodic Riccati equation
            if calculate_riccati
                options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
                [ts,ys] = ode45(f_plane,[0,10],[0;1.17], options);
                i = 1;
                phi_dot = [0 0];
                length(ys)
                while ys(i,1) <= pi
                    phi_dot(i,:) = ys(i,:);
                    i = i+1;
                    if i >= length(ys)
                        break
                    end
                end
                phi_dot(i,:) = ys(i,:);
                %% Interpolated 
                interpolation_for_dphi = @(x) interp1(phi_dot(:,1),phi_dot(:,2), x);
                obj.function_for_dphi = @(phi) interpolation_for_dphi(mod(phi,pi));
                % Using the curvfitted function since ode45 has problems with
                % interpolations.
                figure()
                plot(ys(:,1),ys(:,2))
                [A, B] = obj.get_linearization(phi, obj.function_for_dphi,true);          
                obj.A = @(phi) A(phi);
                obj.B = @(phi) B(phi);
                [X,phi] = sdp_riccati(obj.A,obj.B,obj.Q,obj.Gamma,0,pi,100,40,3);
                obj.X = X;
                obj.phi_test = phi;
                %% Interpolation for Riccati solution
                interpolation_for_X = @(x)[interp1(phi,reshape(X(1,1,:),1,length(X)),x) interp1(phi,reshape(X(1,2,:),1,length(X)),x) interp1(phi,reshape(X(1,3,:),1,length(X)),x);
                                          interp1(phi,reshape(X(2,1,:),1,length(X)),x) interp1(phi,reshape(X(2,2,:),1,length(X)),x) interp1(phi,reshape(X(2,3,:),1,length(X)),x);
                                          interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x);];
                obj.function_for_X = @(x) interpolation_for_X(mod(x,pi));
            end
        end
        
        function delta = get_delta(obj, varphi)
            delta = obj.delta(obj.Phi(varphi));
        end
        function tau = get_tau(obj,varphi)
            tau = obj.tau(obj.Phi(varphi));
        end
        function normal_vector = get_normal_vector(obj,varphi)
            normal_vector = obj.normal_vector(obj.Phi(varphi));
        end
        function rho = get_rho(obj, varphi)
           rho = obj.rho(obj.Phi(varphi));
        end
        function l_rho = get_l_rho(obj, varphi)
            l_rho = obj.length_rho(obj.Phi(varphi));
        end
        function kappa = get_kappa(obj, varphi)
            kappa = obj.kappa(obj.Phi(varphi));
        end
        function ds = get_ds(obj, varphi)
            ds = obj.ds(obj.Phi(varphi));
        end
        function dds = get_dds(obj, varphi)
            dds = obj.dds(obj.Phi(varphi));
        end
        function dsf = get_dsf(obj, varphi)
            dsf = obj.dsf(obj.Phi(varphi));
        end
        function ddsf = get_ddsf(obj, varphi)
            ddsf = obj.ddsf(obj.Phi(varphi));
        end
        function theta = get_theta(obj, varphi)
            theta = obj.theta(varphi);
        end
        function d_theta = get_dtheta(obj, varphi)
            d_theta = obj.d_theta(varphi);
        end
        function dd_theta = get_ddtheta(obj, varphi)
            dd_theta = obj.dd_theta(varphi);
        end
        function d_alpha = get_dalpha(obj,varphi)
            rxt = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            rxk = [0 0 1]*cross(obj.get_rho(varphi),obj.get_kappa(varphi));
            d_s = obj.get_ds(varphi);d_sf = obj.get_dsf(varphi);
            dd_s = obj.get_dds(varphi);dd_sf = obj.get_ddsf(varphi);
            dth = obj.get_dtheta(varphi);
            ddth = obj.get_ddtheta(varphi);
            JRm = obj.J_s/obj.m_b/obj.R_b;
            d_alpha = obj.m_b*((dd_s*rxt+rxk*d_s^2-dd_sf*JRm)*dth + ...
                (d_s*rxt-d_sf*JRm)*ddth + 2*dd_s + 2*dd_sf*JRm/obj.R_b);
        end
        function energy = get_energy(obj, q, dq)
            v_b = obj.diff_R(q(1))*dq(1)*obj.get_rho(q(2))...
                +obj.R(q(1))*obj.get_tau(q(2))*obj.get_ds(q(2))*dq(2)
            omega_b = [0;0;dq(1)-obj.get_dsf(q(2))/obj.R_b*dq(2)]
            omega_frame = dq(1);
            kinetic = 1/2*obj.m_b*(v_b')*v_b+1/2*obj.J_s*(omega_b')*omega_b+obj.J_s*omega_frame^2;
            potential = obj.m_b*[0 obj.g 0]*obj.R(q(1))*obj.get_rho(q(2));
            energy = [potential;kinetic];
        end
        
        function gamma_diff = get_dgamma(obj,varphi)
            gamma_diff = obj.m_b*obj.g*[0 1 0]*(...
                obj.diff_R(obj.get_theta(varphi))*obj.get_dtheta(varphi)*obj.get_tau(varphi)*obj.get_ds(varphi)+...
                obj.R(obj.get_theta(varphi))*obj.get_kappa(varphi)*obj.get_ds(varphi)^2+...
                obj.R(obj.get_theta(varphi))*obj.get_tau(varphi)*obj.get_dds(varphi));
        end
        
        function M = get_M(obj, q)
            varphi = q(2);
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi), obj.get_tau(varphi));
            m11 = obj.get_l_rho(varphi)^2+obj.J_f/obj.m_b+obj.J_s/obj.m_b;
            m12 = obj.get_ds(varphi)*rhoxtau-obj.get_dsf(varphi)*obj.J_s/obj.m_b/obj.R_b;
            m22 = obj.get_ds(varphi)^2+obj.get_dsf(varphi)^2*obj.J_s/obj.m_b/obj.R_b^2;
            M =   obj.m_b*[m11 m12;
                           m12 m22];
        end
        
        function C = get_C(obj, q, dq)
            varphi = q(2);
            ds = obj.get_ds(varphi);
            dsf = obj.get_dsf(varphi);
            dds = obj.get_dds(varphi);
            ddsf = obj.get_ddsf(varphi);
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi), obj.get_tau(varphi));
            taudotrho = obj.get_tau(varphi)'*obj.get_rho(varphi);
            rhoxkappa = [0 0 1]*cross(obj.get_rho(varphi),obj.get_kappa(varphi));
            c11 = obj.m_b*taudotrho*ds*dq(2);
            
            c12 = obj.m_b*(taudotrho*ds*dq(1) ...
                +(dds*rhoxtau-ddsf*obj.J_s/obj.m_b/obj.R_b ...
                +ds^2*rhoxkappa)*dq(2));
            
            c21 = -obj.m_b*ds*taudotrho*dq(1);
            
            c22 = obj.m_b*(ds*dds+dsf*ddsf*obj.J_s/obj.R_b^2/obj.m_b)*dq(2);
            
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
            k = linspace(0,pi,500);
            %constraint_eq_1 = obj.get_ds(0)*(1+obj.J_s/obj.m_b*obj.get_p(0)^2)/([0 0 1]*cross(obj.get_rho(0),obj.get_tau(0))+obj.J_s/obj.m_b*obj.get_p(0));
            %constraint_eq_2 = obj.get_ds(pi/2)*(1+obj.J_s/obj.m_b*obj.get_p(pi/2)^2)/([0 0 1]*cross(obj.get_rho(pi/2),obj.get_tau(pi/2))+obj.J_s/obj.m_b*obj.get_p(pi/2));
            JRm = obj.J_s/obj.R_b/obj.m_b;
            not_eq_constraint = @(x) -(obj.get_ds(x)^2+obj.get_dsf(x)^2*JRm/obj.R_b)/ ...
                (obj.get_ds(x)*([0 0 1]*cross(obj.get_rho(x),obj.get_tau(x))) - obj.get_dsf(x)*JRm);
            hold on;
            constraint_gamma =@(x) -obj.get_tau(x)'*obj.diff_R(x)'*(obj.R(x)*obj.get_kappa(x)*obj.get_ds(x)-obj.R(x)*obj.get_tau(x)*obj.get_dds(x)/obj.get_ds(x));
            diff_gamma = @(x) [0 1 0]*(obj.diff_R(obj.get_theta(x))*obj.get_dtheta(x)*obj.get_tau(x)*obj.get_ds(x)+...
                obj.R(obj.get_theta(x))*obj.get_kappa(x)*obj.get_ds(x)^2+obj.R(obj.get_theta(x))*obj.get_tau(x)*obj.get_dds(x));
            results = zeros(length(k),1);
            diff_theta = zeros(length(k),1);
            find_eq_points_atan_tau = zeros(length(k),1);
            for i = 1:length(k)
                find_eq_points_atan_tau(i,1) = atan2(-[0 1 0]*obj.get_tau(k(i)),[1 0 0]*obj.get_tau(k(i)));
                results(i,1) = not_eq_constraint(k(i));
                diff_theta(i,1) = obj.get_dtheta(k(i));
            end
            plot(k,results);
            plot(k,obj.get_theta(k));
            plot(k, unwrap(find_eq_points_atan_tau));
            eq_points = [0;pi/16;pi/4;pi/2;pi/2+pi/4;pi-pi/16];
            scatter(eq_points,arrayfun(constraint_gamma, eq_points));
            not_eq_constraint(pi/2)
            constraint_gamma(0)
            constraint_gamma(pi/2)
            plot(k,diff_theta);
            grid on;
            legend("Constraint all $\varphi$", "Constraint eq. point  $\varphi = n\pi$", "$\Theta'(\varphi)$");
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
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            rhodottau = obj.get_rho(varphi)'*obj.get_tau(varphi);
            ds_ = obj.get_ds(varphi); dsf_ = obj.get_dsf(varphi);
            dds_ = obj.get_dds(varphi); ddsf_ = obj.get_ddsf(varphi);
            JRm = obj.J_s/obj.R_b/obj.m_b;
            alpha = obj.m_b*((ds_*rhoxtau-dsf_*JRm)*obj.get_dtheta(varphi)+ds_^2+dsf_^2*JRm/obj.R_b);
            beta = obj.m_b*((ds_*rhoxtau-dsf_*JRm)*obj.get_ddtheta(varphi)-rhodottau*ds_*obj.get_dtheta(varphi)^2 ...
                + ds_*dds_ + dsf_*ddsf_*JRm/obj.R_b);
            gamma = obj.m_b*[0 obj.g 0]*obj.R(obj.get_theta(varphi))*obj.get_tau(varphi)*ds_;
            
            a = [alpha;beta;gamma];
        end
        function g = get_g_w_y_on_trajectory(obj, varphi, d_varphi)
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            g_w = -obj.m_b*(obj.get_ds(varphi)*rhoxtau-obj.get_dsf(varphi)*obj.J_s/obj.m_b/obj.R_b);
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
        function plot_phase_plane(obj, x_0,x_end,y_0,y_end,title,time)
            a = @(x) obj.alpha_beta_gamma(x);
            f_plane = @(t,x) [x(2);-1/(a(x(1))'*[1;0;0])*((a(x(1))'*[0;1;0])*x(2)^2+a(x(1))'*[0;0;1])];
            phase_plot_2_interactive(f_plane,[x_0 x_end;y_0 y_end],time,title,[20,20],0.1);
        end
        function plot_abg(obj)
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            phi = linspace(0,2*pi,1000);
            abg = zeros(3,1,length(phi));
            for i = 1:length(phi)
                abg(:,1,i) = obj.alpha_beta_gamma(phi(i));
            end
            figure
            k = 1
            y_label = ["$\alpha$" "$\beta$" "$\gamma$"]
            for i = 1:3
                subplot(3,1,i)
                hold on;
                plot(phi,reshape(abg(i,1,:),1,length(phi)));
                set(gca,'XTick',0:pi/4:2*pi) 
                set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$' ...
                    ,'$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
                k = k+1;
                grid on;
                xlabel('$\varphi$', 'interpreter', 'latex')
                ylabel(y_label(i), 'interpreter', 'latex')
                
            end
            grid on
            sgtitle('$\alpha, \beta, \gamma$ of $\phi$', 'Interpreter','latex');
        end  
    end
end