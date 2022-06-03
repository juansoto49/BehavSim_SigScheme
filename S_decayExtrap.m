function [Snew,fnew]=S_decayExtrap(S,f,fpad,psv_enf)
% [Snew]=S_decayExtrap(S,freq,fpad)
% This function performs padding for a S matrix

f=f(:);
% form the new frequency array
df=f(2)-f(1); n=length(f);

N=1+floor((fpad-f(1))/df);
fstart=f(1); fend=fstart+(N-1)*df;
fext=[f(n)+5*df:df:fend]';
fnew=[fstart:df:fend]'; 

m=size(S,1); 


for i=1:m
        for j=1:m
            
            H(1:n,1)=S(i,j,:); Hmag=abs(H); Hphi=unwrap(angle(H));
            if Hmag(1)==0; Hmag(1)=1e-30; end;
            Hdb=20*log10(Hmag); Hdba=Hdb(n); Hphia=Hphi(n);            
            %[i,j,Hdb(1:5)']
            
           % if any(Hdb(n-4:n)==-Inf)
           if any(Hdb(n-4:n)<-600)
                Hext=[H;zeros(length(fext),1)];
                tmp=interp1([f;fext],Hext,fnew,'spline');
                Snew(i,j,:)=tmp;
            else            
                % Magnitude, using the last 5 points to extrapolate
                p=polyfit(f(n-4:n)*1e-6,Hdb(n-4:n),1);
                if p(1) < 0
                    Hdbext=p(1)*fext*1e-6+p(2);   
                else
                    Hdbb=-100;
                    k=(Hdba-Hdbb)/(fend-f(n));
                    Hdbext=Hdba-(fext-f(n))*k;
                end
                Hdbext=[Hdb;Hdbext];
                
                idx_small=find(Hdbext<-600); 
                Hdbext(idx_small)=-600; 

                % Phase using the last 5 points to extrapolate
                p=polyfit(f(n-4:n)*1e-6,Hphi(n-4:n),1);
                if p(1) < 1e-5 % 2*pi*10e8*1e-6*1e-5=0.02 pi; also most constant
                    Hphiext=p(1)*fext*1e-6+p(2);
                else
                    Hphib=Hphia-10*pi;
                    k=(Hphia-Hphib)/(fend-f(n));
                    Hphiext=Hphia-(fext-f(n))*k;    
                end
                Hphiext=[Hphi;Hphiext];

                Hmagext=10.^(Hdbext/20);
                Hext=Hmagext.*exp(sqrt(-1)*Hphiext);
                
                % now apply the hilbert transform to make the paded part
                % causal. 
                
                
             tmp=interp1([f;fext],Hext,fnew,'spline');
             
% %below line is using the hilbert transform to real and imag       
%             tmp_full=[tmp; flipud(tmp(2:end))];
%             tmp_hlbt=hilbert(real(tmp_full)); 
%             
%             mid_idx=round(N/2); qt_idx=round(N*0.75);
%             % if orignal sequence (n) is less than half of the resulted sequency (N), 
%             % change the phase at 0.75*N. 
%             if n<mid_idx  
%                 tmp(qt_idx+1:N)=conj(tmp_hlbt(qt_idx+1:N));
%             else
%             % if original sequence (n) is more than half of the resulted sequence (N), 
%             % change the phase immediately after n; .
%             tmp(n+1:N)=conj(tmp_hlbt(n+1:N));
%             end
                
                
            Snew(i,j,:)=tmp(1:N);    
                
                
            end
            % now apply the hilbert transform to make the paded part
            % causal.



% % below line is using the hilbert transform to mag and phase (minimum
% % phase)
% % do not use the minimum phase since most systems are not minimum phase system; 
% 
% 
%              
%              logHmag=log(abs(tmp));
%              idx=find(~isfinite(logHmag)); 
%              if ~isempty(idx)
%                  logHmag(idx)=logHmag(idx(1)-1);
%              end
%              H1=[logHmag; flipud(logHmag(2:end-1))];
%              phase_hlbt=-imag(hilbert(H1));
%              Hmin=abs(tmp).*exp(sqrt(-1)*phase_hlbt(1:N));
%              
%                mid_idx=round(N/2); qt_idx=round(N*0.75);
%                 if n<mid_idx
%                     tmp(qt_idx+1:N)=Hmin(qt_idx+1:N);
%                 else
%                     tmp(n+1:N)=Hmin(n+1:N);
%                 end 
%                  
%                 %tmp(n+1:N)=Hmin(n+1:N);
%             
%             %pause;
             
              
             
             
            %Snew(i,j,:)=interp1([f;fext],Hext,fnew,'spline');
            
        end
end

if psv_enf
    if ~ispassive(Snew)
        %Snew=makepassive(Snew);
        eigenmax=0.99999;
        Snew=passivityscalingfix(Snew, eigenmax); 
    end
end

% if ~ispassive(Snew) && psv_enf
%     Snew=makepassive(Snew); 
% end



