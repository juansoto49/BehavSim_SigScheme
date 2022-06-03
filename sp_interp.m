function S=sp_interp(S0,f0,f)
[m,m,n]=size(S0);
for k=1:m
    for j=1:m
        skj=squeeze(S0(k,j,:));
        skj_mag=abs(skj);
        skj_phs=unwrap(angle(skj));
        mag=interp1(f0,skj_mag,f);
        phs=interp1(f0,skj_phs,f);
        S1(k,j,:)=mag.*exp(sqrt(-1)*phs);
    end
end
S=S1;