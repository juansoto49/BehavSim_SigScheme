function [p_out,p_out_mat,cursorUI,cursorRow]=pulse_centering(t_in,p_in,steps_per_UI,opt)

% t_in = input time
% p_in = input puasle
% opt = sampling method, 1=peak, 2=Mueller-Muller, 3=max(cursor-abs(1st-precurrsor) 

[pr_sampled, t_sampled]=find_sampled_pulse(t_in,p_in,steps_per_UI,opt);
[cursor,id]=max(pr_sampled);
cursor_idx=find(t_in==t_sampled(id));
left_idx=cursor_idx-floor(steps_per_UI/2);
right_idx=left_idx+steps_per_UI;
pre_cursorUI=length(left_idx:-steps_per_UI:1);
%cursorUI=pre_cursorUI;
dp1=min(left_idx:-steps_per_UI:1);
dp2=max(right_idx:steps_per_UI:length(p_in(:,1)))-1;
t_out=t_in(dp1:dp2);
p_out=p_in(dp1:dp2);
p_out_mat=reshape(p_out,steps_per_UI,[]);
for i=pre_cursorUI-1:pre_cursorUI+1
    cursorRow=find(p_out_mat(:,i)==cursor);
    if ~isempty(cursorRow)
        cursorUI=i;
        break
    end 
end
%{
% Find pulse peak
[curval curind] = max(p_in);

% Move one UI left from the peak
if curind-steps_per_UI >= 1
    dp1=curind-steps_per_UI;
else
    dp1=1;
end

% sweep from dp1 to the pulse peak such that p(left)=p(left+steps_per_UI) 
indx=dp1:curind;
dv=abs(p_in(indx)-p_in(indx+steps_per_UI));
vm=p_in(indx+steps_per_UI/2);
[dv_min,dv_ind]=min(dv);
left=indx(dv_ind);
right=left+steps_per_UI;
curval_new=vm(dv_ind);
curind_new=indx(dv_ind)+steps_per_UI/2;

% Determine starting point relative to cursor UI
dp1 = mod(left,steps_per_UI); %data point 1 - ie first data point
if dp1==0,
    left = left + 1;
    dp1 = mod(left,steps_per_UI);
end
  
% Cursor UI
cursorUI = (left-dp1)/steps_per_UI + 1; %UI # of cursor

% Post-cursor UI count
N_post_UI = floor((length(p_in)-right)/steps_per_UI); %# of post cursor UIs

% Total UI
nUI = (N_post_UI+cursorUI);    

% End point
dp2 = right + steps_per_UI*N_post_UI; %data point 2 - last data point

% Output 
p_out=p_in(dp1:dp2-1);
mUI=floor(length(p_out)/steps_per_UI);
pr_mat=reshape(p_out(1:mUI*steps_per_UI),steps_per_UI,[]);
%}
