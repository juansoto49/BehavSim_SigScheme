% NRZ, PAM4, QPSK and QAM16 Behavioral Simulation Top
clear all
close all
clc

%addpath('C:\User\path_to_function\Functions\');

%Inputs for S-Parameter Pulse Response Extraction
%------------------------------------------------------------------------------------------------------------------------------
res.UI = 31.25e-12;             %Select Bit Unit Interval (depending on your signaling scheme
res.Tr = 3.125e-12;             %Select Rise Time
Vsw=1;                          %Select voltage swing for the ideal tx mmodel
fpad=200e9;                     %Define padding frequency for ifft
xtlk_scale=2;                   %Define crosstalk scaling factor
f=[0:10e6:100e9]';              %Range of frequencies for input S-Parameter Model
Zref=42.5;                      %Define reference impedance for S-Parameter Model
file_name='C:\User\File_to_sp\spfile.snp';  %Input S-Parameter File
%-------------------------------------------------------------------------------------------------------------------------------
%Inputs for EQ Optimization Algorithm and PDA Analysis
%-------------------------------------------------------------------------------------------------------------------------------
UI=res.UI;                      %Reusing pulse response Bit Unit Interval (do not modify); 
steps_per_UI=64;                %Steps for each Bit Invertal
noisefloor=0;                   %Noise floor required for PDA
Npre=5;                         %Number of pre-cursor ISI to include on analysis
Npost=20;                       %Number of post-cursor ISI to include on analysis included
add_FFE=1;                      %Add FFE Equalization
add_CTLE=1;                     %Add CTLE Equalization
add_DFE=1;                      %Add DFE Equalization
txffe_table=[                   %Table for possible taps for FFE Equalization (User can increase or reduce this table and amount of taps. This demo has two-precursors and one-postcursor
0	0	1	0
0	-0.083	0.917	0
0	-0.167	0.833	0
0	0	0.917	-0.083
0	0	0.833	-0.167
0.042	-0.208	0.75	0
0.042	-0.125	0.708	-0.125
0.083	-0.208	0.709	0
0.083	-0.250	0.667	0
0.083	-0.250	0.625	-0.042
];
% CTLE - pole-zero model. This is the pulse-zero model  2-pole 1-zero.
% These parameters can be modified. To modify the algorithm. Modify the
% CTLE function source code
fz1 = 250e6;
g_dc=[-15:1:-5];
fz3 = 7.7e9;
fp1 = 1.3*fz1;
fp2 = 7.7e9;
fp3 = 22e9;
fp4 = 28e9;
fp5 = 32e9;
fp6 = 32e9;
% % CTLE - pole-zero model This can be selected using the gen_CTLE_ir
% function instead of the current implemented
% fp1=5e9;
% fp2=10e9;
% fz=3.55e9;
% g_dc=[-9:1:0];
%This is the possible DFE taps available. User can select amount of taps.
%Algorithm will only use the desired taps
Ndfe=3;
sig_target=0.1; %This is the AGC Magnitude Target
dfe_bound=[
-0.030 0.030
-0.020 0.020
-0.020 0.020
-0.015 0.015
-0.010 0.010
-0.005 0.005
];
%----------------------------------------------------------------------------------------------------------------------------------------
%Inputs for Bit-by-Bit Simulation
%----------------------------------------------------------------------------------------------------------------------------------------
Modulation_Selection = 'NRZ';   %Select Between NRZ, QPSK, QAM16 or PAM4
N = 10000;                      %Number of bits to simulate


%----------------------------------------------------------------------------------------------------------------------------------------
%Pulse Response Extraction from S-Parameter
%--------------------------------------------------------------------------
S=read_SParam(file_name,f,Zref); % load TS file, interpolate, extrapolate and renormalization
sdd=s2sdd(S.Parameters,1);  % Convert to diff S-parameters
sdd21=squeeze(sdd(7,2,:));  % Diff IL
sdd31=squeeze(sdd(7,6,:));  % Diff NEXT
% Pad to fpad
%--------------------------------------------------------------------------
tmp(1,1,:)=sdd21;
[Snew,fext]=S_decayExtrap(tmp,f,fpad,0);
sdd21_pad=squeeze(Snew);
tmp(1,1,:)=sdd31;
[Snew,fext]=S_decayExtrap(tmp,f,fpad,0);
sdd31_pad=squeeze(Snew);
% Input pulse spectrum
%--------------------------------------------------------------------------
f1=fext; 
if f1(1)==0; f1(1)=1e-30; end; % avoid f(1)=0.
s=complex(0,2*pi*f1); 
Vin=1./(res.Tr*s.^2).*(1.-exp(-s*res.Tr)).*(1.-exp(-s*res.UI));
% Method 2 - get pulse response using impulse response
%--------------------------------------------------------------------------
% IL impulse
P=sdd21_pad;
h_tmp = real(ifft([P; flipud(conj(P(2:end)))])); % impulse response
h = h_tmp(1:length(fext));
h_il=h;
% Input pulse
t0=[0; res.Tr; res.UI; res.UI+res.Tr];
v0=[0;  1;  1; 0];
% Interpolate input pulse
dt=0.5/fext(end); 
t=[0:length(fext)-1]'*dt;
t_int=[0:dt:t0(end)];
v0_int=interp1(t0,v0,t_int);
% IL pulse and step
p2(:,1)=filter(v0_int,1,h)*Vsw;
%------------------------------------------------
% Crosstalk pulse
P=sdd31_pad;
h_tmp = real(ifft([P; flipud(conj(P(2:end)))])); % impulse response
h = h_tmp(1:length(fext));
p2(:,2) = filter(v0_int,1,h)*xtlk_scale*Vsw;
figure
title('Crosstalk pulse response')
plot(t,p2(:,2),'LineWidth',2)
xlabel('Time,s')
ylabel('Pulse Response, V')
legend('Method 2')
grid on
figure
title('IL pulse response')
plot(t,p2(:,1),'LineWidth',2)
xlabel('Time,s')
ylabel('Pulse Response, V')
legend('Method 2')
grid on
res.t = t;
res.p = p2;

p=res.p;                        %Pulse response results from S-Parameter Extraction
t=res.t;                        %Time axis for pulse response result.
tstep=UI/steps_per_UI;

% Resample pulse response
t_new=[0:tstep:t(end)]';
for i=1:size(p,2)
    p_new(:,i)=interp1(t,p(:,i),t_new,'pchip');
end

% Trunction
[~,peakind]=max(p_new(:,1));
ID=find((t_new>=t_new(peakind)-Npre*UI) & (t_new<=t_new(peakind)+Npost*UI));
if ~isempty(ID)
    t_uneq=t_new(ID); 
    t_uneq=t_uneq-t_uneq(1);
    p_uneq(:,1)=p_new(ID,1);
end
for i=2:size(p_new,2)
    p_uneq(:,i)=p_new(1:length(t_uneq),i);
end
figure
plot(t_uneq,p_uneq(:,1))
xlabel('Time,s')
ylabel('Pulse Response, V')


%-------------------------------------------------------------------------

if add_FFE
    Nffe=size(txffe_table,1);
else
    Nffe=1;
end

if add_CTLE
    Nctle=length(g_dc);
else
    Nctle=1;
end

% Adaptive to sweep all possible EQ settings
eye_height_matrix = [];
eye_height_best=-1.0;
for i_ffe=1:Nffe

    % Apply TXFFE
    coeff=txffe_table(i_ffe,:);
    for i=1:size(p_uneq,2)
        [t_eq,p_ffe(:,i)]=conv_by_shift_add(tstep,p_uneq(:,i),coeff,steps_per_UI);
    end
    
    for i_ctle=1:Nctle
        
        % Apply CTLE
        if add_CTLE
            ir=get_CTLE_gen6_ir_basic(fz1,fz3,fp1,fp2,fp3,fp4,fp5,fp6,g_dc(i_ctle),t_eq);
            %ir=get_CTLE_ir_basic(fz,fp1,fp2,g_dc(i_ctle),t_eq);
            for i=1:size(p_ffe,2)
                p_ctle(:,i)=filter(ir,1,p_ffe(:,i));
            end
        else
            p_ctle=p_ffe;
        end
        
        % Find the sampling point or cursor
        pr=p_ctle(:,1);
        %[pr_sampled, sample_idx, curind]=find_sampled_pulse(pr,steps_per_UI,2);
        [pr_sampled, t_sampled]=find_sampled_pulse(t_eq,pr,steps_per_UI,2);
        %{
        figure
        plot(t_eq,pr)
        hold all
        plot(t_sampled,pr_sampled,'O')
        %}
        % Cursor
        pmax=max(pr);
        
        % AGC
        [cursor,curind]=max(pr_sampled);
        AGC=sig_target/pmax;
        p_ctle=p_ctle*AGC;
        pr_sampled=pr_sampled*AGC;
        cursor=AGC*cursor;
        
        % Add DFE
        if add_DFE
            for i=1:Ndfe
                postcur=pr_sampled(curind+i);
                if postcur > dfe_bound(i,2)
                    dfe_taps(i)=dfe_bound(i,2);
                    pr_sampled(curind+i)= postcur - dfe_bound(i,2);
                elseif postcur < dfe_bound(i,1)
                    dfe_taps(i)=dfe_bound(i,1);
                    pr_sampled(curind+i)= postcur - dfe_bound(i,1);
                else
                    dfe_taps(i)=postcur;
                    pr_sampled(curind+i)= 0;
                end
            end
        else
            dfe_taps=[];
        end
        
        pr_sampled(curind)=0;
        if curind+Npost+1 < length(pr_sampled)
            pr_sampled(curind+Npost+1:end)=0;
        end
        if curind-Npre-1 > 1
            pr_sampled(1:curind-Npre-1)=0;
        end
        
        % Total ISI
        ISI_total=sum(abs(pr_sampled));
        
        % Crosstalk
        for i=1:size(p_ctle,2)-1
            nUI=floor(size(p_ctle,1)/steps_per_UI);
            pr_xtlk=p_ctle(1:nUI*steps_per_UI,i+1);
            pr_xtlk_mat=reshape(pr_xtlk,steps_per_UI,[]);
            XTK(i)=max(sum(abs(pr_xtlk_mat),2));
        end
        
        % Total crosstalk
        XTK_total=sum(XTK);
       
        eye_height=cursor-(ISI_total+XTK_total);
        %for debug_start
        eye_height_matrix(i_ffe) = eye_height; 
        %for debug_end
        if eye_height > eye_height_best
            eye_height_best=eye_height;
            cursor_best=cursor;
            ISI_total_best=ISI_total;
            XTK_total_best=XTK_total;
            coeff_best=coeff;
            g_dc_best=g_dc(i_ctle);
            p_ctle_best=p_ctle;   
            dfe_taps_best=dfe_taps;
        end
        
    end
end

figure
plot(t_uneq,p_uneq(:,1),t_eq,p_ctle_best(:,1),'LineWidth',2)
xlabel('Time,s')
ylabel('Pulse Response, V')
grid on
legend('Unequalize','FFE+CTLE Equalized')

%--------------------------------------------------------------------------

% Center the eye
opt=2;
[pr, pr_mat, cursorUI]=pulse_centering(t_eq,p_ctle_best(:,1),steps_per_UI,opt);
% pr_mat: Re-arrange pulse response such that each row is a sampled response and
% pr_mat(:,cursorUI) represent the possible cursors. 
%pr_mat=reshape(pr,steps_per_UI,[]);

% Apply DFE
if add_DFE    
    for i=1:Ndfe        
        pr_mat(:,cursorUI+i)=pr_mat(:,cursorUI+i)-dfe_taps_best(i);
    end
end
pr_dfe=reshape(pr_mat,[],1);
t_dfe=[0:length(pr_dfe)-1]'*tstep;

[pr_dfe_max,pr_dfe_max_idx] = max(pr_dfe);
dfe_marker_idxs = [pr_dfe_max_idx-2*steps_per_UI-1:-steps_per_UI:1 pr_dfe_max_idx-(steps_per_UI)-1 pr_dfe_max_idx pr_dfe_max_idx+steps_per_UI-1 pr_dfe_max_idx+2*steps_per_UI-1:steps_per_UI:length(t_dfe)];

figure
[p_uneq_max,p_uneq_max_idx] = max(p_uneq(:,1));
p_uneq_marker_idxs = [p_uneq_max_idx-steps_per_UI:-steps_per_UI:1 p_uneq_max_idx p_uneq_max_idx+steps_per_UI-1 p_uneq_max_idx+2*steps_per_UI-1:steps_per_UI:length(t_uneq)];
[ctle_max,ctle_max_idx] = max(p_ctle_best(:,1));
ctle_marker_idxs = [ctle_max_idx-2*steps_per_UI-1:-steps_per_UI:1 ctle_max_idx-(steps_per_UI)-1 ctle_max_idx ctle_max_idx+steps_per_UI-1 ctle_max_idx+2*steps_per_UI-1:steps_per_UI:length(t_eq)];
[pr_dfe_max,pr_dfe_max_idx] = max(pr_dfe);
dfe_marker_idxs = [pr_dfe_max_idx-2*steps_per_UI-1:-steps_per_UI:1 pr_dfe_max_idx-(steps_per_UI)-1 pr_dfe_max_idx pr_dfe_max_idx+steps_per_UI-1 pr_dfe_max_idx+2*steps_per_UI-1:steps_per_UI:length(t_dfe)];

plot(t_uneq./UI,p_uneq(:,1),'.-','MarkerIndices',p_uneq_marker_idxs,'MarkerSize',25,'LineWidth',2);
hold on
plot(t_eq./UI-1.75,p_ctle_best(:,1),'.-','MarkerIndices',ctle_marker_idxs,'MarkerSize',25,'LineWidth',2); %Adjusted move to plot the aligned pulse response. The UI movement can be manually adjusted
plot(t_dfe./UI-1.5,pr_dfe,'.-','MarkerIndices',dfe_marker_idxs,'MarkerSize',25,'LineWidth',2);
xlabel('Time, UI')
ylabel('Pulse Response, V')
grid on
legend('Unequalize','FFE+CTLE','FFE+CTLE+DFE')
xlim([0 25]);
hold off

%--------------------------------------------------------------------------
% % Plot eye diagram
% 
if strcmp(Modulation_Selection,'NRZ')

    N=10000;
    pattern=(randi([0 1],N,1)-0.5)/0.5; % generate N-number of random 0's or 1's bipolar (remove -0.5 to unipolar) 
    span = 3;
    sps = 8;
    % Waveform by shifting and adding
    [t_w,v_w]=conv_by_shift_add(tstep,pr_dfe,pattern,steps_per_UI);
    figure
    plot(t_w,v_w,'LineWidth',2)
    grid on
    xlabel('Time,s')
    ylabel('Waveform, V')

    % Chop the waveform with one-UI segements
    mUI=floor(length(v_w)/steps_per_UI);
    v_w=v_w(1:mUI*steps_per_UI);
    wave_mat=reshape(v_w,steps_per_UI,[]);
    x=[0:steps_per_UI-1]'*tstep;
    figure
    plot(x,wave_mat(:,[sps*span+1:end-sps*span]))
    xlabel('Time,s')
    ylabel('NRZ Eye-diagram, V')
    
    %PSD Results
    [pxx,f1] = pwelch(v_w,[],[],[],'centered','power');
    figure
    plot(f1./pi.*1024,pow2db(pxx),'LineWidth',2);
    xlim([-100 100]);
    xlabel('Frequency (GHz)','Fontsize',12);
    ylabel('Power (dB)','Fontsize',12);
    legend('Bipolar NRZ','Fontsize',10);
    %--------------------------------------------------------------------------
    % PDA

    % Center the eye
    [pr, pr_mat, cursorUI]=pulse_centering(t_dfe,pr_dfe,steps_per_UI,2);

    % Cursor
    cursor = pr_mat(:,cursorUI);
    pr_mat(:,cursorUI) = 0;

    % PDA upper and lower curves
    if cursorUI-Npre-1 >= 1
        pr_mat(:,1:cursorUI-Npre-1)=0;
    end
    if cursorUI+Npost+1 < size(pr_mat,2)
        pr_mat(:,cursorUI+Npost+1:end)=0;
    end
    upr=sum(min(pr_mat,-noisefloor),2);
    lwr=sum(max(pr_mat,noisefloor),2);

    x=[0:steps_per_UI-1]'*tstep;
    figure
    plot(x,cursor+upr,'LineWidth',2)
    hold all
    plot(x,lwr,'LineWidth',2)
    xlabel('Time,s')
    ylabel('NRZ PDA Eye, V')
    
elseif strcmp(Modulation_Selection,'PAM4')
    %--------------------------------------------------------------------------
    %PAM4 EYE
    % 
    % PAM4 pattern
    L=4;
    %pattern=randi([0 L-1],N,1)/(L-1); %Unipolar PAM4
    pattern=(randi([0 L-1],N,1)-(L-1)/2)/((L-1)/2); %Bipolar PAM4
    span = 3;
    sps = 2*L;

    % Waveform by shifting and adding
    [t_w,v_w]=conv_by_shift_add(tstep,pr_dfe,pattern,steps_per_UI);

    figure
    plot(t_w,v_w,'LineWidth',2)
    grid on
    xlabel('Time,s')
    ylabel('Waveform, V')

    % Chop the waveform with one-UI segements
    wave_mat=reshape(v_w,steps_per_UI,[]);
    x=[0:steps_per_UI-1]'*tstep;
    figure
    plot(x,wave_mat(:,[sps*span+1:end-sps*span]));
    xlabel('Time,s');
    ylabel('PAM4 Eye-diagram, V');
    
    %PSD Results
    [pxx,f1] = pwelch(v_w,[],[],[],'centered','power');
    figure
    plot(f1./pi.*1024,pow2db(pxx),'LineWidth',2);
    xlim([-100 100]);
    xlabel('Frequency (GHz)','Fontsize',12);
    ylabel('Power (dB)','Fontsize',12);
    legend('PAM4','Fontsize',10);
  
    % Add crosstalk
    for i=1:size(p_new,2)-1
        p_xt=p_new(:,i+1);
        nUI=floor(length(p_xt)/steps_per_UI);
        p_xt=p_xt(1:nUI*steps_per_UI);
        p_xt_mat=reshape(p_xt,steps_per_UI,[]);
        upr_xt(i)=min(sum(min(p_xt_mat,-noisefloor),2));
        lwr_xt(i)=max(sum(max(p_xt_mat,noisefloor),2));
    end
    
    %--------------------------------------------------------------------------
    % PDA

    % Center the eye
    [pr, pr_mat, cursorUI]=pulse_centering(t_dfe,pr_dfe,steps_per_UI,2);

    % Cursor
    cursor = pr_mat(:,cursorUI);
    pr_mat(:,cursorUI) = 0;

    % PDA upper and lower curves
    if cursorUI-Npre-1 >= 1
        pr_mat(:,1:cursorUI-Npre-1)=0;
    end
    if cursorUI+Npost+1 < size(pr_mat,2)
        pr_mat(:,cursorUI+Npost+1:end)=0;
    end
    upr=sum(min(pr_mat,-noisefloor),2);
    lwr=sum(max(pr_mat,noisefloor),2);

    % Add crosstalk
    for i=1:size(p_new,2)-1
        p_xt=p_new(:,i+1);
        nUI=floor(length(p_xt)/steps_per_UI);
        p_xt=p_xt(1:nUI*steps_per_UI);
        p_xt_mat=reshape(p_xt,steps_per_UI,[]);
        upr_xt(i)=min(sum(min(p_xt_mat,-noisefloor),2));
        lwr_xt(i)=max(sum(max(p_xt_mat,noisefloor),2));
    end

    % Sum all crosstalks
    upr_xt_sum=sum(upr_xt);
    lwr_xt_sum=sum(lwr_xt);

    % PDA with crosstalk
    upr_wxt=upr+upr_xt_sum;
    lwr_wxt=lwr+lwr_xt_sum;
    %plot(x,cursor+upr_wxt)
    %plot(x,lwr_wxt)
    %Here starts PAM4 Simulations
    % PAM4 PDA
    %--------------------------------------------------------------------------
    L=4;            % signal level
    a=0:1/(L-1):1;  % bit value
    for i=1:L-1
        pda(i).ISIupr(:,1) = upr(:,1) + a(i+1)*cursor;
        pda(i).ISIlwr(:,1) = lwr(:,1) + a(i)*cursor;
    end

    time_ui=[0:tstep:UI-tstep]';
    figure
    hold all
    plot(time_ui/UI,pda(1).ISIupr)    
    plot(time_ui/UI,pda(1).ISIlwr)
    plot(time_ui/UI,pda(2).ISIupr)    
    plot(time_ui/UI,pda(2).ISIlwr)
    plot(time_ui/UI,pda(3).ISIupr)    
    plot(time_ui/UI,pda(3).ISIlwr)

    % Add crosstalk
    %for i=1:L-1
    %pda(i).ISIupr_wxt(:,1) = upr_wxt(:,1) + a(i+1)*cursor;
    %pda(i).ISIlwr_wxt(:,1) = lwr_wxt(:,1) + a(i)*cursor;
    %end

elseif strcmp(Modulation_Selection,'QPSK')
    %--------------------------------------------------------------------------
    % 
    % %Plot eye diagram of PSK
    M = 4;
    span = 3;
    sps = 2*M;
    pattern=randi([0 M-1],N,1); %generate N-number of random 0's or 1;s
    patternmod= pskmod(pattern,M,pi/4);
    %patternmod = qammod(pattern, M);
    % Waveform by shifting and adding
    [t_w,v_w] = conv_by_shift_add(tstep,pr_dfe,patternmod,steps_per_UI);
    plot(t_w,v_w,'LineWidth',2)
    grid on
    xlabel('Time,s')
    ylabel('PSK Waveform, V')

    txsig = v_w;

    % Chop the waveform with one-UI segements
    mUI=floor(length(v_w)/steps_per_UI);
    v_w=v_w(1:mUI*steps_per_UI);
    wave_mat=reshape(v_w,steps_per_UI,[]);
    x=[0:steps_per_UI-1]'*tstep;
    figure
    subplot(1,2,1)
    plot(x,real(wave_mat(:,[sps*span+1:end-sps*span])));
    xlabel('Time,s');
    ylabel('PSK Eye-diagram In Phase, V');
    subplot(1,2,2);
    plot(x,imag(wave_mat(:,[sps*span+1:end-sps*span])));

    xlabel('Time,s')
    ylabel('PSK Eye-diagram Quadrature, V')

    scatterplot(wave_mat(steps_per_UI/2,[sps*span+1:end-sps*span]));
    
    %PSD Results
    [pxx,f1] = pwelch(v_w,[],[],[],'centered','power');
    figure
    plot(f1./pi.*1024,pow2db(pxx),'LineWidth',2);
    xlim([-100 100]);
    xlabel('Frequency (GHz)','Fontsize',12);
    ylabel('Power (dB)','Fontsize',12);
    legend('QPSK','Fontsize',10);
    
elseif strcmp(Modulation_Selection,'QAM16')
    
    % %Plot eye diagram of QAM
    % %Pattern
    M = 16;
    span = 3;
    sps = 2*M;
    pattern=randi([0 M-1],N,1); %generate N-number of random 4-bit sequences
    patternmod = qammod(pattern, M,'UnitAveragePower',true);

    % Waveform by shifting and adding
    [t_w,v_w] = conv_by_shift_add(tstep,pr_dfe,patternmod,steps_per_UI);
    plot(t_w,v_w,'LineWidth',2)
    grid on
    xlabel('Time,s')
    ylabel('QAM Waveform, V')

    txsig = v_w;

    % Chop the waveform with one-UI segements
    mUI=floor(length(v_w)/steps_per_UI);
    v_w=v_w(1:mUI*steps_per_UI);
    wave_mat=reshape(v_w,steps_per_UI,[]);
    x=[0:steps_per_UI-1]'*tstep;
    figure
    subplot(1,2,1)
    plot(x,real(wave_mat(:,[sps*span+1:end-sps*span])));
    xlabel('Time,s');
    ylabel('QAM Eye-diagram In Phase, V');
    subplot(1,2,2);
    plot(x,imag(wave_mat(:,[sps*span+1:end-sps*span])));

    xlabel('Time,s')
    ylabel('QAM Eye-diagram Quadrature, V')

    scatterplot(wave_mat(steps_per_UI/2,[sps*span+1:end-sps*span]));
    
    %PSD Results
    [pxx,f1] = pwelch(v_w,[],[],[],'centered','power');
    figure
    plot(f1./pi.*1024,pow2db(pxx),'LineWidth',2);
    xlim([-100 100]);
    xlabel('Frequency (GHz)','Fontsize',12);
    ylabel('Power (dB)','Fontsize',12);
    legend('QAM16','Fontsize',10);

    
end




