%*************     main program    *************%
%%%%%%%%%    PCA witness   %%%%%%%%

ww1 = input ('Enter the witness well log 1:');
ww2 = input ('Enter the witness well log 2:');
ww3 = input ('Enter the witness well log 3:');
ww4 = input ('Enter the witness well log 4:');
Topw=input('enter depth of the top of the witness well (m): ');
bound=input('enter the depth of the boundary in witness well (m): '); 
[wwnew,wwnew1,wwnew2] = PCA_witness_logs(ww1,ww2,ww3,ww4);


%%%%%%%%%    PCA observation   %%%%%%%%

ow1 = input('Enter the observation well log 1:');
ow2 = input('Enter the observation well log 2:');
ow3 = input('Enter the observation well log 3:');
ow4 = input('Enter the observation well log 4:');
Top=input(' enter depth of the top of the observation well (m): ');
realdepth=input('enter the real depth of selected boundary(m): ');


[ownew] = PCA_observation_logs(ow1,ow2,ow3,ow4);

%% window length

disp('choose one of the following amount for the length of the window (m) : (1.218752 , 2.437504 , 4.875008 , 9.750016 , 19.500032)');
winlength=input('length of the window: ');

%% 1st calculation of the parameter in witness well

x1=wwnew';
depth=bound-Topw;
N=input('numebr of decomposition level: ');


%% fractal and statistical features calculation

[Hws1]=Wit_Frc_calc_logs(x1,N,depth,winlength);
[ Avg1, Cv1, R1, theta1]=Wit_Sta_calc_logs(x1,N,depth,winlength);

%% 2nd calculation of the parameter in observation well

x2=ownew';

%% sliding window

sigdis=0.152344;                                                                          % distance between each data (meter)
windis=sigdis;                                                                             % length of the sliding window movement
%winlength=2^(N+2)*sigdis;                                                       % lenght of the window (meter)  
win=floor((length(x2)-(winlength/sigdis))/(windis/sigdis)+1);               % number of windows

%% calculation of the features in sliding window along well log

[ Avg,Cv,R, theta]=Statis_calc_logs(x2,N,winlength);
[Hws]= Fractal_calc_logs(x2,N,winlength);

%% witness and observation features matrix

WFeatures=[Hws1,Avg1, Cv1, R1, theta1];       % witness  features matrix
OFeatures=[Hws,Avg,Cv,R, theta];                        % observation features matrix
StdFeatures=std(OFeatures);                                                    % standard deviation of the features along the observation features matrix
Num=length(WFeatures);                                                          % number of features

[Prob,Probability]=Prob_calc_logs(x2,N,WFeatures,OFeatures,StdFeatures,Num,win,winlength,sigdis);

%% probability plots

for i=1:win-(winlength/sigdis)
Depth(i,1)=Top+windis*(i-1)+winlength;
end

plot(Probability,Depth);
title(' total probability (independence assumption) ', 'color' , 'r')
xlabel('total probability');
ylabel('depth (m)');

figure
subplot(2,4,1)
plot(Prob(1:end,1),Depth)                                          % Hws probability vs. depth
title('subplot 1: Hws probability')
subplot(2,4,2)
plot(Prob(1:end,2),Depth)                                          % Average probability vs. depth
title('subplot 2: average probability')
subplot(2,4,3)
plot(Prob(1:end,3),Depth)                                          % coefficient of variation probability vs. depth
title('subplot 3: coefficient of variation probability')
subplot(2,4,4)
plot(Prob(1:end,4),Depth)                                          % min/max ratio probability vs. depth
title('subplot 4: min/max ratio probability')
subplot(2,4,5)
plot(Prob(1:end,5),Depth)                                          % theta probability vs. depth
title('subplot 5: theta ratio probability')

%% outputs

aim=find(Probability==(max(max(Probability))));
aimdepth=Top+windis*(aim-1)+winlength;             % depth estimated
difdepth=aimdepth-realdepth;                                    % difference between estimated depth from the real depth
Real=ceil((realdepth-Top-winlength)/sigdis);
for j=Real-5:Real+5
Realprob(j,1:5)=Prob(j,1:5);                                            % feature probability of the real depth
PROB=max(Realprob);
end






