function S= read_SParam(filename,freq,Zref);
% Read s-parameter and remaps it over frequency grid
%
rf_obj = read(rfckt.datafile,filename);
S=sparameters(rf_obj);
f=S.Frequencies;
sp=S.Parameters;
m=size(sp,1);
n=length(f);

f0=f; S0=sp; 

% Extrapolate toward DC
if f(1) > freq(1)
    f0(1)=freq(1); f0(2:n+1)=f;
    S0(1:m,1:m,1)=sp(1:m,1:m,1);    
    S0(1:m,1:m,2:n+1)=sp;
end

% Extrapolate toward infinity
n=length(f0);
if f(end) < freq(end)
    f0(n+1)=freq(end);
    S0(1:m,1:m,n+1)=sp(1:m,1:m,end);
end

% Interpolate
sp_int=sp_interp(S0,f0,freq);    
S=sparameters(sp_int,freq,S.Impedance);    

% Renorm
if Zref ~= S.Impedance
    S=newref(S,Zref);
end    
