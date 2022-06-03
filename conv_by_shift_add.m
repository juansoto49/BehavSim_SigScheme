function [t,v]=conv_by_shift_add(dt,v0,pattern,steps_per_UI)

% Contruct waveform based on pattern and pulse response
%
% Input
%      dt = time step
%      v0 = unequalized pulse response
%      pattern = data pattern
%      steps_per_UI = steps per UI
%
% Output
%      t = time
%      v = resulting waveform

% Make sure the input is a row vector
if size(v0,1) ~= 1
    v_tmp=v0';
else
    v_tmp=v0;
end

% Shift the pulse response
n_sft=length(pattern);
v=0;
for i=1:n_sft
    if pattern(i) ~= 0
        v=v+pattern(i)*[zeros(1,(i-1)*steps_per_UI) v_tmp zeros(1,(n_sft-i)*steps_per_UI)];
    end
end
t=[0:length(v)-1]*dt;
if size(v0,1) ~= 1
    t=t'; v=v';
end


