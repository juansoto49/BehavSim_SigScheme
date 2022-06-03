function ir=get_CTLE_ir_basic(fz,fp1,fp2,g_dc,t)
% Get CTLE impulse response for classitc 2-pole-1-zero transfer function
% This is based on the analytical solution.

w_p1=2*pi*fp1;
w_p2=2*pi*fp2;
w_z=2*pi*fz;

d0=w_p1*w_p2/w_z;
a0=w_z*10^(g_dc/20);

a=(a0-w_p1)/(w_p2-w_p1);
b=(w_p2-a0)/(w_p2-w_p1);
    
ir=d0*(a*exp(-w_p1*t)+b*exp(-w_p2*t))*(t(2)-t(1));