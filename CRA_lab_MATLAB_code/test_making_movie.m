h = figure;
time = out.tout(end);
fps = 30;
number_of_frames = round(time*fps);
x_frame = zeros(1,length(out.tout));
y_frame = zeros(1,length(out.tout));
x_ball = bf.R_b*cos(linspace(0, 2*pi,20));
y_ball = bf.R_b*sin(linspace(0, 2*pi,20));
xlim([-0.2 0.2])
ylim([-0.2 0.2])
axis equal;
grid on;
M(number_of_frames) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';
j = 0;
for i = linspace(0,2*pi)
    j = j + 1;
    position = bf.delta(i);
    x_frame(j) = position(1);
    y_frame(j)= position(2);
end
frame = hgtransform;
ball = hgtransform;
patch('XData',x_frame,'YData',y_frame,'FaceColor','yellow','Parent',frame) 
patch('XData',x_ball,'YData',y_ball,'FaceColor','red','Parent',ball) 
start_falling = length(out.after_fall(:,1))-length(find(out.after_fall(:,1)));
acumulated_time = 0;
current_frame = 0;
for i = 2:start_falling
    acumulated_time = acumulated_time + out.tout(i) - out.tout(i-1);
    if acumulated_time >= 1/fps
        ball_position = bf.rho(out.q(i,2));
        ball_position_inertial = bf.R(out.q(i,1))*ball_position;
        frame.Matrix = makehgtform('zrotate',out.q(i,1));
        ball.Matrix = makehgtform('translate', ball_position_inertial);
        drawnow
        current_frame = current_frame + 1
        M(current_frame) = getframe;
        acumulated_time = acumulated_time - 1/fps;
    end
end
for i = start_falling+1:length(out.after_fall(:,1))
   acumulated_time = acumulated_time + out.tout(i)-out.tout(i-1);
   if acumulated_time >= 1/fps
        frame.Matrix = makehgtform('zrotate',out.after_fall(i,1));
        ball.Matrix = makehgtform('translate', [out.after_fall(i,2);out.after_fall(i,3);0]);
        drawnow
        current_frame = current_frame + 1
        M(current_frame) = getframe;
        acumulated_time = acumulated_time - 1/fps;
    end
end
h.Visible = 'on';
movie(M,1,30);