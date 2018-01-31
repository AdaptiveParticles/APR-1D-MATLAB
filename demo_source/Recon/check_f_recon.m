function [linear_inf,pc_inf,wc_inf] = check_f_recon(apr,plot_flag)

y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^apr.l_max):apr.s_dom(2);

%% Reconstruct and check

f_r = interp1(apr.y_p,apr.f_p,y);
f_gt = apr.f(y);

linear_inf = max(abs(f_gt-f_r)/apr.scale(1));

if(plot_flag == true)
    figure;plot(y,f_r)
    
    hold on
    plot(y,f_gt)
    xlabel('y')
    ylabel('f(y)')
    title('Linear Reconstruciton')
    
    disp(['Desired: ',num2str(apr.E),' Linear Interpolation E_obs: ',num2str(linear_inf)])

end

%% Nearest Neighbour Reconstruction

f_r = interp1(apr.y_p,apr.f_p,y,'nearest');
f_gt = apr.f(y);

pc_inf = max(abs(f_gt-f_r)/apr.scale(1));

if(plot_flag == true)
    figure;plot(y,f_r)
    
    hold on
    plot(y,f_gt)
    xlabel('y')
    ylabel('f(y)')
    title('NN Piece-wise Constant Reconstruciton')
    
    disp(['Desired: ',num2str(apr.E),' Nearest Neighbour E_obs: ',num2str(pc_inf)])
end





%% Worse Case Reconstruction

[max_rc,min_rc,y] = worst_case(apr);

recon_max = max(abs(apr.f(y)-max_rc)/apr.scale(1));
recon_min = max(abs(apr.f(y)-min_rc)/apr.scale(1));

wc_inf = max(recon_max,recon_min);

if(plot_flag == true)
    figure;plot(y,max_rc);hold on
    plot(y,min_rc);
    xlabel('y')
    ylabel('f(y)')
    title('Worst Case Reconstruciton Bounds')
    
    disp(['Desired: ',num2str(apr.E),' Worst Case E_obs: ',num2str(wc_inf)])
end






end