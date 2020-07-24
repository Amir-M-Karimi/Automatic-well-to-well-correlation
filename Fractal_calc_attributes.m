%*****       Observation Well Fractal Parameter Calculation       *****%

function [Hws]= Fractal_calc_attributes(Signal,N,winlength);

%% properties of the signal

signum=length(Signal);                            % number of data in the signal
sigdis=0.152344;                                            % distance between each data (meter)
siglength=signum*sigdis;                       % lenght of the well that well log has run (meter)

%% properties of the sliding window

% if N is the numebr of the decomposition level
% winlength=2^(N+2)*sigdis;                        % lenght of the window (meter)
windis=sigdis;                                        % length of the sliding window movement

%% wavelet decomposition                                              

[C,L]=wavedec(Signal,N,'haar');           % N level decomposition of Signal by using 'db1' wavelet
% C is the coefficient matrix that contains the coefficient of details of
% each level and approximation of the last level decomposition
% L is the bookkeeping vector
% approximation and details formation

for j=1:N;
    D{N-j+1}=C(sum(L(1:j,1))+1:sum(L(1:j+1,1)));            % wavelet detalis and approximation matrix 
end
A=C(1:L(1,1),1);

%% fractal parameter estimation

%win=floor((length(Signal)-2^(N+2))/(windis/sigdis)+1);               % number of windows
win=floor((length(Signal)-(winlength/sigdis))/(windis/sigdis)+1);               % number of windows
for jj=1:N
    
    detnum=length(D{jj});                                              % number of data of details
    %winnum=2^(N+2-jj);                                                % number of data in window
    winnum=winlength/(sigdis*2^jj);                            % number of data in a window
    disnum=windis/((2^jj)*sigdis);                                % number of data between each movement
    %FFTResults{jj}= fft(D{jj});                                                                             % fast fourier transform of details
    %Power{jj} = FFTResults{jj}.*conj(FFTResults{jj})/D{jj};                        % power matrix calculation
   
    z=1;
    
     FFTR{jj}=ones(winnum,win);
     
    for ii=1:win
        
         FFTR{jj}(1:end,ii)=fft(D{jj}(z:z+winnum-1,1));
         Power{jj}=FFTR{jj}(1:end,ii).*conj(FFTR{jj}(1:end,ii))/(D{jj}(z:z+winnum-1,1));
         Energy(ii,jj)= norm(Power{jj});
         
         V(ii,jj)=var(D{jj}(z:z+winnum-1,1));                                                         % variance calculation in each window
         S(ii,jj)=std(D{jj}(z:z+winnum-1,1));                                                         % standard deviation in each wnidow
         
       %  Energy(ii,jj)= norm(Power{jj}(z:z+winnum-1,1:end));                           % energy matrix calculation
       
         

  if mod(ii,2^jj)==0  
            z=z+1;
          end
        
    end
end

for jj=1:N
   levlog(jj)=log2(2^jj);                                                                                        % log2 of the decomposition level
end

    enlog=log2(Energy);                                                                                      % log2 of the energy of data in each window
    varlog=log2(V);                                                                                               %  log2 of the variance of data in each window
    stdlog=log2(S);                                                                                               %  log2 of the standard deviation of data in each window
    
 for kk=1:win             
    [r2(kk,1),m2(kk,1),b2(kk,1)]=regression(levlog,stdlog(kk,1:N));               % slope= Hws parameter
 end
 
 
 Hws=m2;
end