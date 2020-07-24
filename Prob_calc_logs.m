%***********      Probability Calculation       ************%

function [Prob,Probability]=Prob_calc_logs(Signal,N,WFeatures,OFeatures,StdFeatures,Num,win,winlength,sigdis);

%% Error Calculation

for i=1:win-(winlength/sigdis);
Upper_win_Err (i,:)=WFeatures(1,:)-OFeatures(i,:);                                         % Error calculation for the upper window
Lower_win_Err (i,:)=WFeatures(2,:)-OFeatures(i+(winlength/sigdis),:);                      % Error calculation for the lower window
end

for j=1:win-(winlength/sigdis);
Up_F_Prob(j,:)= (1./(sqrt(2*pi.*StdFeatures))).*exp(-(Upper_win_Err(j,:).^2)./(2.*StdFeatures));                            %upper window probability of each feature
L_F_Prob(j,:)= (1./(sqrt(2*pi.*StdFeatures))).*exp(-(Lower_win_Err (j,:).^2)./(2.*StdFeatures));                              %lower window probability of each feature
end

%% normalizing

for ii=1:win-(winlength/sigdis);
    for jj=1:Num
norm_Up_F_Prob(ii,jj)=Up_F_Prob(ii,jj)/(max(Up_F_Prob(:,jj)));          % normalized probability of each features in upper window
norm_L_F_Prob(ii,jj)=L_F_Prob(ii,jj)/(max(L_F_Prob(:,jj)));                   % normalized probability of each features in lower window
    end
end


   for jj=1:win-(winlength/sigdis);
            probab(jj,:)=Up_F_Prob(jj,:).*L_F_Prob(jj,:);
        end
        
      
for ii=1:win-(winlength/sigdis);
    for jj=1:Num
Prob(ii,jj)=probab(ii,jj)/(max(probab(:,jj)));          % normalized probability of each features 
    end
end
                                            

Up_Prob=ones(win-(winlength/sigdis),1);
L_Prob=ones(win-(winlength/sigdis),1);

for ii=1:win-(winlength/sigdis);
    for jj=1:Num
 Up_Prob(ii,1)=Up_Prob(ii,1)*norm_Up_F_Prob(ii,jj);                    % Probability of upper window          
 L_Prob(ii,1)=L_Prob(ii,1)*norm_L_F_Prob(ii,jj);                             % Probability of lower window
    end
end

%% independent assumption

Probability=Up_Prob.*L_Prob;                                               % total probability ( if the features be independent)
end

















