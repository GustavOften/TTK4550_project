classdef butterfly_robot
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
        r_f = 12.5e-3; %Half of the distance between both frame plates
        R_b = sqrt(16.55e-3^2-12.5e-3^2); %Radius of the ball
        delta
        tau_delta
        rho
        length_rho
        bad_rho
        alpha
        curve_for_phi        
        ds
        dds
        tau
        kappa
        R
        diff_R
        tau_diff 
        theta
        diff_theta
        diff_diff_theta
        K
        function_for_X
        function_for_dphi
        Gamma = 0.1
        Q = @(t) [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        X
        phi_test
        A
        B
    end
    methods
        function obj = butterfly_robot(calculate_riccati)
            %% varphi = f(phi), phi = g(varphi)
            delta_curve_fit = @(phi)(obj.a - obj.b*cos(2*phi))*[sin(phi);cos(phi);0];
            tau_curve_fit = @(phi)((-2*obj.b*sin(2*phi))*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0])/ ...
                norm((-2*obj.b*sin(2*phi))*[sin(phi);cos(phi);0]+(obj.a - obj.b*cos(2*phi))*[cos(phi);-sin(phi);0]);
            
            range_for_functions = linspace(0,2*pi,100);
            function_f = @(x) atan2([1 0 0]*delta_curve_fit(x) - obj.R_b*[0 1 0]*tau_curve_fit(x), ...
                                [0 1 0]*delta_curve_fit(x)+obj.R_b*[1 0 0]*tau_curve_fit(x));
            result_function_f = zeros(length(range_for_functions),1);
            for i = 1:length(range_for_functions)
                result_function_f(i) = function_f(range_for_functions(i));
            end
            result_function_f = unwrap(result_function_f);
            p = polyfit(result_function_f, range_for_functions, 7);
            figure
            plot(result_function_f,range_for_functions)
            hold on;
            g_varphi = matlabFunction(poly2sym(p))
            plot(result_function_f, g_varphi(result_function_f));
            %phi_polynomial = poly2sym(p, phi);
            
            syms varphi phi(varphi) theta real
            %phi = poly2sym(p, varphi);
            phi = varphi;
            %Delta: parametrization of the shape of the frame.
            
            delta = (obj.a - obj.b*cos(2*phi))*[sin(phi);cos(phi);0];
            obj.delta = matlabFunction(delta);
            tau = diff(delta,varphi)/norm(diff(delta,varphi));
            obj.tau = matlabFunction(tau);
            
            %% Rho: Vector to center of ball, given in bodyframe.
            
            bad_rho = (obj.a-obj.b*cos(2*phi)+obj.R_b)*[sin(phi);cos(phi);0];
            obj.bad_rho = matlabFunction(bad_rho);
            rho = delta + obj.R_b*[-tau(2);tau(1);0];
            rho = rewrite(rho, 'sqrt'); % atan2 makes it complex...
            obj.rho = matlabFunction(rho);
            obj.length_rho = matlabFunction(norm(rho));
            
            
            %% s: Arclength of the balls path.
            % s = integral(|diff(rho,phi)|) from 0 to phi
            % Fundamental theorem of calculus gives the derivative of s
            % with respect to phi.
            ds = norm(diff(rho, varphi))
            dds = diff(ds, varphi)
            obj.ds = matlabFunction(ds);
            obj.dds = matlabFunction(dds);

            %% Kappa: Curvature of curve.
            kappa = diff(rho,varphi,varphi)/norm(diff(rho,varphi))^2
            obj.kappa = matlabFunction(kappa);
            
            %% Rotational matrices:
            obj.R = matlabFunction([cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]);
            obj.diff_R = matlabFunction([-sin(theta) -cos(theta) 0;cos(theta) -sin(theta) 0;0 0 0]);
            
            %%%%%%%%%%%%%%%%%%%%%%%% Control stuff %%%%%%%%%%%%%%%%%%%%%%
            Rtau = obj.R(phi)*obj.tau(phi);
            %%%%%%%%%%%%%% Theta used in Case-study non-prehensile %%%%%%%
            a = -0.03;
            numerator = a*sin(2*phi-pi)*Rtau(1)-Rtau(2)
            denominator = a*sin(2*phi-pi)*Rtau(2)+Rtau(1)
            %Theta = simplify(phi+atan(numerator/denominator))
            a= -1.007046e+00; b= 6.398405e-01; c= 2.507379e-01; d= 1.164677e-01;
            Theta = a*atan(b*sin(2*varphi)+c*sin(4*varphi)+d*sin(6*varphi))+varphi;
            dTheta = simplify(diff(Theta,varphi))
            ddTheta = simplify(diff(dTheta,varphi))
            %% Theta used in Internship report
            %Theta = phi-1.3*sin(2*phi);
            %Theta = phi-0.2*sin(Rtau(1));
            obj.theta = matlabFunction(Theta);
            obj.diff_theta = matlabFunction(dTheta);
            obj.diff_diff_theta = matlabFunction(ddTheta);

            %% Plots phase plane of alpha betta gamma function
            a = @(x) obj.alpha_beta_gamma(x);
            f_plane = @(t,x) [x(2);-1/(a(x(1))'*[1;0;0])*((a(x(1))'*[0;1;0])*x(2)^2+a(x(1))'*[0;0;1])];
            %phase_plot_2_interactive(f_plane,[-pi pi;-1 3],10,'',[100,100],0.1)
            %% Finding solution to periodic Riccati equation
            if calculate_riccati
                options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
                [ts,ys] = ode45(f_plane,[0,10],[0;1.5], options);
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
                obj.A = @(phi) A(phi);
                obj.B = @(phi) B(phi);
                [X,phi] = sdp_riccati(obj.A,obj.B,obj.Q,obj.Gamma,0,pi,500,5,4);
                obj.X = X; obj.phi_test = phi;
                %% Interpolation for Riccati solution
                interpolation_for_X = @(x)[interp1(phi,reshape(X(1,1,:),1,length(X)),x) interp1(phi,reshape(X(1,2,:),1,length(X)),x) interp1(phi,reshape(X(1,3,:),1,length(X)),x) interp1(phi,reshape(X(1,4,:),1,length(X)),x);
                                          interp1(phi,reshape(X(2,1,:),1,length(X)),x) interp1(phi,reshape(X(2,2,:),1,length(X)),x) interp1(phi,reshape(X(2,3,:),1,length(X)),x) interp1(phi,reshape(X(2,4,:),1,length(X)),x);
                                          interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x) interp1(phi,reshape(X(3,4,:),1,length(X)),x);
                                          interp1(phi,reshape(X(4,1,:),1,length(X)),x) interp1(phi,reshape(X(4,2,:),1,length(X)),x) interp1(phi,reshape(X(4,3,:),1,length(X)),x) interp1(phi,reshape(X(4,4,:),1,length(X)),x)];
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
        
        function M = get_M(obj, q)
            rhoxtau = cross(obj.rho(q(2)), obj.tau(q(2)));
            m11 = obj.m_b*obj.length_rho(q(2))^2+obj.J_f+obj.J_s;
            m12 = (obj.m_b*rhoxtau(3)-obj.J_s/obj.R_b)*obj.ds(q(2));
            m22 = (obj.m_b+obj.J_s/obj.R_b^2)*obj.ds(q(2))^2;
            M =   [m11 m12;
                   m12 m22];
        end
        
        function C = get_C(obj, q, dq)
            phi = q(2);
            rhoxtau = cross(obj.rho(q(2)), obj.tau(q(2)));
            taudotrho = obj.tau(phi)'*obj.rho(phi);
            rhoxkappa = cross(obj.rho(phi),obj.kappa(phi));
            c11 = obj.m_b*obj.ds(phi)*taudotrho*dq(2);
            c12 = obj.m_b*(obj.ds(phi)*taudotrho*dq(1)+(obj.ds(phi)^2*rhoxkappa(3)+(rhoxtau(3)-obj.J_s/obj.R_b)*obj.dds(phi))*dq(2));
            c21 = -obj.m_b*obj.ds(phi)*taudotrho*dq(1);
            c22 = obj.m_b*(1+obj.J_s/obj.m_b/obj.R_b^2)*obj.ds(phi)*obj.dds(phi)*dq(2);
            C = [c11 c12;
                 c21 c22];
        end
        
        function G = get_G(obj, q)
            G = [obj.m_b*[0;obj.g;0]'*obj.diff_R(q(1))*obj.rho(q(2));
                 obj.m_b*[0;obj.g;0]'*obj.R(q(1))*obj.tau(q(2))*obj.ds(q(2))];
        end
        
        function ddq = calculate_ddq(obj, q, dq, u)
            C = obj.get_C(q,dq);
            G = obj.get_G(q);
            M_inverse = obj.get_M(q)^-1;
            ddq = M_inverse*(-C*dq-G+[u;0]);
        end
        
        function normal_force = get_normal_force(obj, q)
            %Using a simplified model. 
            % Assuming the normal vector is 90 deg anti-clockwise to tau.
            n = [0 -1 0;1 0 0;0 0 1]*obj.tau(q(2));
            normal_force = [0;9.81;0]'*obj.R(q(1))*n;
        end
        function position_ball = get_position_ball(obj, q)
            ball_position = obj.rho(q(2));
            position_ball = obj.R(q(1))*ball_position;
        end
        function velocity_ball = get_velocity_ball(obj, q, dq)
            txRrho = cross([0;0;dq(1)],obj.R(q(1))*obj.rho(q(2)));
            velocity_ball = txRrho + obj.R(q(1))*obj.tau(q(2))*obj.ds(q(2))*dq(2);
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
        
        function a = alpha_beta_gamma(obj, phi)
            rxt = cross(obj.rho(phi),obj.tau(phi));
            rxk = cross(obj.rho(phi),obj.kappa(phi));
            sum1 = rxt(3)-obj.J_s/obj.m_b/obj.R_b;
            sum2 = 1+obj.J_s/obj.m_b/obj.R_b^2;
            Alpha = obj.m_b*obj.ds(phi)*(sum1*obj.diff_theta(phi)+sum2*obj.ds(phi));
            beta = obj.m_b*obj.ds(phi)*(sum1*obj.diff_diff_theta(phi)-obj.tau(phi)'*obj.rho(phi)*obj.diff_theta(phi)^2+sum2*obj.dds(phi));
            rotated_tau = obj.R(obj.theta(phi))*obj.tau(phi);
            gamma = obj.m_b*[0;obj.g;0]'*rotated_tau*obj.ds(phi);
            a = [Alpha;beta;gamma];
        end
        function g = get_g_w_y_on_trajectory(obj, phi, phi_dot)
            rhoxtau = cross(obj.rho(phi),obj.tau(phi));
            g_w = -obj.m_b*obj.ds(phi)*(rhoxtau(3)-obj.J_s/obj.m_b/obj.R_b);
            g_y_dot = obj.m_b*obj.ds(phi)*obj.tau(phi)'*obj.rho(phi)*(2*obj.diff_theta(phi)*phi_dot);
            g_y = -obj.ds(phi)*obj.m_b*obj.g*[cos(obj.theta(phi)) -sin(obj.theta(phi)) 0]*obj.tau(phi);
            g = [g_w;g_y_dot;g_y];
        end
        function [A, B] = get_linearization(obj, phi, phi_dot, function_handles)
            if function_handles == true
                abg = @(x)obj.alpha_beta_gamma(x);
                gwy = @(x)obj.get_g_w_y_on_trajectory(x, phi_dot(x));
                A = @(x)[0 1 0 0;
                         0 0 0 0;
                        [0 0 1]*gwy(x)/([1 0 0]*abg(x)) [0 1 0]*gwy(x)/([1 0 0]*abg(x)) ([0 0 1]*abg(x)-[0 1 0]*abg(x)*phi_dot(x)^2)/([1 0 0]*abg(x))/phi_dot(x) 0;
                         0 0 1 0];
                B = @(x)[0;1;[1 0 0]*gwy(x)/([1 0 0]*abg(x));0];
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
            L = [1 obj.diff_theta(q(2));
                 0 1];
            N = [obj.diff_diff_theta(q(2))*dq(2)^2;
                 0];
            M = obj.get_M([obj.theta(q(2));q(2)]);
            C = obj.get_C([obj.theta(q(2));q(2)],[obj.diff_theta(q(2))*dq(2);dq(2)]);
            G = obj.get_G([obj.theta(q(2));q(2)]);
            denom = (L^-1)*(M^-1);
            numerator = L^-1*(N+(M^-1)*C*L*[epsilon(2);dq(2)]+(M^-1)*G)*1;
            numerator = numerator(1)+obj.get_w(obj.function_for_X(q(2)),q, epsilon);
            u = numerator/denom(1,1);
        end
        function u = get_my_u(obj, q, dq, epsilon)
            M = obj.get_M([obj.theta(q(2));q(2)]);
            C = obj.get_C([obj.theta(q(2))+epsilon(1);q(2)],[obj.diff_theta(q(2))*dq(2)+epsilon(2);dq(2)]);
            G = obj.get_G([obj.theta(q(2))+epsilon(1);q(2)]);
            k = [1 0]*M*[obj.diff_theta(q(2));1]/([obj.diff_theta(q(2)) 1]*M*[0;1]);
            u = k*([k^-1 -1]*M*[obj.get_w(obj.function_for_X(q(2)),q, epsilon)+obj.diff_diff_theta(q(2))*dq(2)^2;0]+[k^-1 -1]*C*[(obj.diff_theta(q(2)))*dq(2)+epsilon(2);dq(2)]+[k^-1 -1]*G);
            %u = numerator/denom(1,1);

        end
        function epsilon = get_epsilon(obj, q, dq)
            epsilon = [q(1)-obj.theta(q(2));dq(1)-obj.diff_theta(q(2))*dq(2);dq(2)-obj.function_for_dphi(q(2))];
        end
    end
end