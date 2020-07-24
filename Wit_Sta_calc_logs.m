%*****     Witness Well Statistical Parameter Calculation      *****%

function [ Avg1, Cv1, R1, theta1]=Wit_Sta_calc_logs(Signal,N,depth,winlength);

%% properties of the signal

signum=length(Signal);                            % number of data in the signal
sigdis=0.152344;                                            % distance between each data (meter)
siglength=signum*sigdis;                       % lenght of the well that logged (meter)

%% properties of the window

% winlength=2^(N+2)*sigdis;                                     % lenght of the window (meter)

windis=sigdis;                                                     % length of the sliding window movement
%winnum=2^(N+2);                                                    % number of data in window
winnum=winlength/sigdis;                      % number of data in window
disnum=windis/sigdis;                                           % number of data between each movement

%% depth of the boundary

winboundary=ceil(depth/windis);       % number of the window that has been placed lower the boundary

%% satistical pararmeter calculation


 strt=ceil((depth/windis)*disnum);                      % the boundary location
 z=1;
  
  for jj=1:winnum
      regx(jj,1)=jj;
  end

       Avg1(1,1)=mean(Signal(strt-winnum:strt-1,1));                                                  % The first feature, average value of  lower window
       Avg1(2,1)=mean(Signal(strt:strt+winnum-1,1));                                                 % The first feature, average value of upper window
       
       Cv1(1,1)=((var(Signal(strt-winnum:strt-1,1)))^0.5)/Avg1(1,1);                          % The second feature, the coefficient of variation of  lower window
       Cv1(2,1)=((var(Signal(strt:strt+winnum-1,1)))^0.5)/Avg1(2,1);                          % The second feature, the coefficient of variation of upper window
       
       Max(1,1)=max(Signal(strt-winnum:strt-1,1));
       Max(2,1)=max(Signal(strt:strt+winnum-1,1));
       
       Min(1,1)=min(Signal(strt-winnum:strt-1,1));
       Min(2,1)=min(Signal(strt:strt+winnum-1,1));
       
       R1(1,1)=Min(1,1)/Max(1,1);                                                                   % The min/max ratio, R
       R1(2,1)=Min(2,1)/Max(2,1);
       
       m(1,:)=polyfit(regx,Signal(strt-winnum:strt-1,1),1);
       m(2,:)=polyfit(regx,Signal(strt:strt+winnum-1,1),1);
       
    for ii=1:2
        if m(ii,1)<=0
           theta1(ii,1)=-atan(m(ii,1))/pi;
       else
           theta1(ii,1)=1-(atan(m(ii,1))/pi);
        end
    end
        z=z+2;
  end