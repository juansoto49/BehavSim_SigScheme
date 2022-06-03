function ir=get_CTLE_gen6_ir_basic(fz1,fz3,fp1,fp2,fp3,fp4,fp5,fp6,g_dc,t1)
% Get CTLE impulse response for classitc 2-pole-1-zero transfer function
% This is based on the analytical solution.

wp1=2*pi*fp1; %c
wp2=2*pi*fp2; %d
wp3=2*pi*fp3; %e
wp4=2*pi*fp4; %f
wp5=2*pi*fp5; %g
wp6=2*pi*fp6; %h
wz1=2*pi*fz1; %a
wz2=wp2*10^(g_dc/20); %k
wz3=2*pi*fz3; %b

%d0 = (wp1*wp3*wp4*wp5*wp6)/(wz1*wz3);

%a = (wz1*wz3*wp1-wz1*wz3*wz2-wz1*wp1^2+wz1*wp1*wz2-wz3*wp1^2+wz3*wp1*wz2+wp1^3-wp1^2*wz2)/((wp1-wp2)*(wp1-wp3)*(wp1-wp4)*(wp1-wp5)*(wp1-wp6));
%b = (-wz1*wz3*wp2+wz1*wz3*wz2+wz1*wp2^2-wz1*wp2*wz2+wz3*wp2^2-wz3*wp2*wz2-wp2^3+wp2^2*wz2)/((wp1-wp2)*(wp2-wp3)*(wp2-wp4)*(wp2-wp5)*(wp2-wp6));
%c = (-wz1*wz3*wp3+wz1*wz3*wz2+wz1*wp3^2-wz1*wp3*wz2+wz3*wp3^2-wz3*wp3*wz2-wp3^3+wp3^2*wz2)/((wp1-wp3)*(wp3-wp2)*(wp3-wp4)*(wp3-wp5)*(wp3-wp6));
%d = (-wz1*wz3*wp4+wz1*wz3*wz2+wz1*wp4^2-wz1*wp4*wz2+wz3*wp4^2-wz3*wp4*wz2-wp4^3+wp4^2*wz2)/((wp1-wp4)*(wp4-wp2)*(wp4-wp3)*(wp4-wp5)*(wp4-wp6));
%e = (-wz1*wz3*wp5+wz1*wz3*wz2+wz1*wp5^2-wz1*wp5*wz2+wz3*wp5^2-wz3*wp5*wz2-wp4^3+wp5^2*wz2)/((wp1-wp5)*(wp5-wp2)*(wp5-wp3)*(wp5-wp4)*(wp5-wp6));
%f = (-wz1*wz3*wp6+wz1*wz3*wz2+wz1*wp6^2-wz1*wp6*wz2+wz3*wp6^2-wz3*wp6*wz2-wp4^3+wp6^2*wz2)/((wp1-wp6)*(wp6-wp2)*(wp6-wp3)*(wp6-wp4)*(wp6-wp5));

%ir=d0*(a*exp(-wp1*t)+b*exp(-wp2*t)+c*exp(-wp3*t)+d*exp(-wp4*t)+e*exp(-wp5*t)+f*exp(-wp6*t))*(t(2)-t(1));

syms s t
H = ((wp1*wp3*wp4*wp5*wp6)*(s+wz1)*(s+wz2)*(s+wz3))/(wz1*wz3*(s+wp1)*(s+wp2)*(s+wp3)*(s+wp4)*(s+wp5)*(s+wp6));
h = ilaplace(H);
ir = double(subs(h,t,t1))*(t1(2)-t1(1));

