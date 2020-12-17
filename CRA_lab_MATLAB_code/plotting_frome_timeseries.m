function plotting_frome_timeseries(ts)
    %% dvarphi vs varphi
    figure
    hold on;
    xlim([-pi/16 2*pi])
    plot(ts.q(:,2),ts.dq(:,2));
    plot(ts.q(:,2),ts.dphi_star(:))
    set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    grid on;
    legend('$\dot{\varphi}$','$\dot{\varphi}_\star$')
    sgtitle('$\dot{\varphi}$ and $\varphi$')
    xlabel('$\varphi$[rad]');
    ylabel('$\dot{\varphi}[rad/s]$')
    exportgraphics(gcf,'varphi_vs_dvarphi.png')
    
    %% Theta vs varphi
    figure
    hold on;
    xlim([-pi/16 2*pi])
    plot(ts.q(:,2),ts.q(:,1));
    plot(ts.q(:,2),ts.Theta(:))
      set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    grid on;
    legend('$\vartheta$','$\Theta$')
    sgtitle('$\vartheta$ and $\Theta$ ')
    ylabel('$\vartheta$[rad]')
    xlabel('$\varphi$[rad]');
    exportgraphics(gcf,'theta_vs_varphi.png')
    
    %% Epsilon 
    figure
    hold on;
    y_label = ["$y$" "$\dot{y}$" "$z$"];
    for i = 1:3
        subplot(3,1,i)
        hold on;
        plot(ts.q(:,2),ts.epsilon(:,i));
        set(gca,'XTick',0:pi/2:2*pi) 
        set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
        xlim([0 2*pi])
        ylabel(y_label(i), 'interpreter', 'latex')
        grid on;
        xlabel('$\varphi$[rad]');
    end
    sgtitle('Transverse coordinates')
    exportgraphics(gcf,'epsilon.png')
    
    %% Epsilon as function of time
    figure
    hold on;
    y_label = ["$y$" "$\dot{y}$" "$z$"];
    for i = 1:3
        subplot(3,1,i)
        hold on;
        plot(ts.tout(:),ts.epsilon(:,i));
        ylabel(y_label(i), 'interpreter', 'latex')
        grid on;
        xlabel('time[s]');
    end
    sgtitle('Transverse coordinates')
    exportgraphics(gcf,'epsilon_time.png')
    
    %% U, input
    figure
    hold on;
    xlim([-pi/16 2*pi])
    plot(ts.q(:,2),ts.u(:));
    plot(ts.q(:,2),ts.u_nominal(:))
      set(gca,'XTick',0:pi/2:2*pi) 
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
    grid on;
    legend('u', '$u_\star$', 'interpreter', 'latex');
    sgtitle('Input')
    ylabel('Input [Nm]')
    xlabel('$\varphi$[rad]');
    exportgraphics(gcf,'input.png')
    
    %% U, input functions of time
    figure
    hold on;
    %xlim([-pi/16 2*pi])
    plot(ts.tout(:),ts.u(:));
    plot(ts.tout(:),ts.u_nominal(:))
    grid on;
    legend('u', '$u_\star$', 'interpreter', 'latex');
    sgtitle('Input')
    ylabel('Input [Nm]')
    xlabel('time[s]');
    exportgraphics(gcf,'input_time.png')
    
    
    %% varphi as a function of time 
    figure
    hold on;
    plot(ts.tout(:),ts.q(:,2));
    grid on;
    sgtitle('$\varphi(t)$')
    xlabel('time[s]');
    ylabel('$\varphi(t)$')
    exportgraphics(gcf,'varphi_function_of_time.png')
end