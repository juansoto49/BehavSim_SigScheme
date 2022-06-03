function [pr_sampled, t_sampled]=find_sampled_pulse(t,pr,steps_per_UI,opt)

% opt = 1: use the pulse peak as the sampling point
% opt = 2: zero 1st precursor
% opt = 3: use maximize cursor-abs(1st precursor) 

% Find pulse peak
[curval curind] = max(pr);

% Check if opt 2 is valid. If not, set opt=3
if opt == 2
    ii=find(pr(1:curind)<0,1,'last');    
    if isempty(ii) || pr(ii+steps_per_UI) ~= max(pr(ii:steps_per_UI:end))        
        disp('Failed to find the zero-1st precursor; swich to option 1')
    else
        if abs(pr(ii)) > abs(pr(ii+1)) % pick the smaller one as the 1st pre-cursor
            ii=ii+1;
        end
        curind=ii+steps_per_UI;        
    end        
end

if opt == 3
    % Search from left to right such that abs(cursor-1st precursor) is maximum
    indx=(curind-steps_per_UI):(curind+steps_per_UI/2);    
    cur_precur=pr(indx)-abs(pr(indx-steps_per_UI));
    [~,ind]=max(cur_precur);
    curind=indx(ind); 
    curval=pr(curind);
end

% Sampled pulse
preind=fliplr(curind-steps_per_UI:-steps_per_UI:1);
postind=curind+steps_per_UI:steps_per_UI:length(pr);
smpind=[preind curind postind];
pr_sampled=pr(smpind);
t_sampled=t(smpind);
