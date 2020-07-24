%*****     Witness Well Statistical Parameter Calculation      *****%

function [ Avg1, Cv1, R1, theta1]=Wit_Sta_calc_attributes (Signal,N,winlength);

%% properties of the signal
signum=length(Signal);                            % number of data in the signal
sigdis=0.152344;                                            % distance between each data (meter)
siglength=signum*sigdis;                       % lenght of the well that loged (meter)

%% properties of the sliding window

% winlength=(winlength/sigdis)*sigdis;                                     % lenght of the window (meter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%windis=2*sigdis;                                                     % length of the sliding window movement
windis=sigdis;                                                         % length of the sliding window movement
winnum=winlength/(sigdis);                         % number of data in window
disnum=windis/sigdis;                                           % number of data between each movement
signum=length(Signal);                                          % number of data in the signal
 
     %%%    if mod(signum,2)==0
        %%% win=floor((signum-winnum)/disnum+1);             % number of window   
      %%% else
         %%% win=floor((signum-winnum)/disnum+2);             % number of window
       %%% end
       win=floor((length(Signal)-(winlength/sigdis))/(windis/sigdis)+1);               % number of windows
       
      %% Definition of features
  z=1;
  
  for jj=1:winnum
      regx(jj,1)=jj;
  end
  
  for ii=1:win
       Avg1(ii,1)=mean(Signal(z:z+winnum-1,1));                                                 % The first feature, average value
       Cv1(ii,1)=((var(Signal(z:z+winnum-1,1)))^0.5)/Avg1(ii,1);                          % The second feature, the coefficient of variation
       Max(ii,1)=max(Signal(z:z+winnum-1,1));
       Min(ii,1)=min(Signal(z:z+winnum-1,1));
       R1(ii,1)=Min(ii,1)/Max(ii,1);                                                                           % The min/max ratio, R
       m(ii,:)=polyfit(regx,Signal(z:z+winnum-1,1),1);
    
        if m(ii,1)<=0
           theta1(ii,1)=-atan(m(ii,1))/pi;
       else
           theta1(ii,1)=1-(atan(m(ii,1))/pi);
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% z=z+2;
       z=z+1;
  end

  end