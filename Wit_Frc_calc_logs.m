%*****      Witness Well Fractal Parameter       *****%

function [Hws1]=Wit_Frc_calc_logs(Signal,N,depth,winlength);

%% properties of the signal

signum=length(Signal);                                    % number of data in the signal
sigdis=0.152344;                                                % distance between each data (meter)
siglength=signum*sigdis;                               % lenght of the well that well log run (meter)

%% properties of the window

%winlength=2^(N+2)*sigdis;                               % lenght of the window (meter)
windis=sigdis;                                                         % length of the sliding window movement

%% depth of the boundary

winboundary=ceil(depth/windis);             % number of the window that has been placed lower the boundary

%% wavelet decomposition   

[C,L]=wavedec(Signal,N,'haar');                      % 5 level decomposition of Signal by using 'db1' wavelet
% C is the coefficient matrix that contains the coefficient of details of
% each level and approximation of the last level decomposition
% L is the bookkeeping vector

% approximation and details formation
for j=1:N;
    D{N-j+1}=C(sum(L(1:j,1))+1:sum(L(1:j+1,1)));
end
A=C(1:L(1,1),1);

%% fractal parameter estimation

for jj=1:N
         winnum=winlength/(sigdis*2^jj);                           % number of data in a window
         disnum=windis/((2^jj)*sigdis);                                % number of data between each movement
         strt=ceil((depth/windis)*disnum);                          % the boundary location
         V(1,jj)=var(D{jj}(strt-winnum:strt-1,1));                 % variance of the pervious window
         V(2,jj)=var(D{jj}(strt:strt+winnum-1,1));                % variance of the next window
         S(1,jj)=std(D{jj}(strt-winnum:strt-1,1));                                                            % standard deviation of the pervious window
         S(2,jj)=std(D{jj}(strt:strt+winnum-1,1));                                                           % standard deviation of the next window
         FFTR1{jj}=fft(D{jj}(strt-winnum:strt-1,1));                                                      % upper window fast fourier transform of details
         FFTR2{jj}=fft(D{jj}(strt:strt+winnum-1,1));                                                     % lower window fast fourier transform of details
         
         power1{jj}=FFTR1{jj}.*conj(FFTR1{jj})/(D{jj}(strt-winnum:strt-1,1));                                        % power matrix calculation
         power2{jj}=FFTR2{jj}.*conj(FFTR2{jj})/(D{jj}(strt:strt+winnum-1,1));
            
         Energy(1,jj)= norm(power1{jj}(1:winnum,1:end));                                                                                                            % energy matrix calculation
         Energy(2,jj)= norm(power2{jj}(1:winnum,1:end));                                                                                                            % energy matrix calculation
        
end

for jj=1:N
   levlog(jj)=log2(2^jj);                                                                                   % log2 of the decomposition level
end

enlog=log2(Energy);                                                                                     % log2 of energy of data in each window
varlog=log2(V) ;                                                                                             %  log2 of variance of data in each window
stdlog=log2(S);                                                                                              %  log2 of standard deviation of data in each window
[r2(1,1),m2(1,1),b2(1,1)]=regression(levlog,stdlog(1,1:N));                      % slope= Hws parameter
[r2(2,1),m2(2,1),b2(2,1)]=regression(levlog,stdlog(2,1:N));                      % slope= Hws parameter

Hws1=m2;
end



