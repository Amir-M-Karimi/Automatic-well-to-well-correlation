%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%main program%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1st calculation witness well
disp('choose one of the following amount for the length of the window (m) : (1.218752 , 2.437504 , 4.875008 , 9.750016 , 19.500032)');
winlength=input('length of the window: ');
N=input('numebr of decomposition level: ');
x1= input ('Enter the witness well log:');
Topw=input('enter depth of the top of the witness well (m): ');
bound=input('enter the depth of the boundary in witness well (m): '); 
depth=bound-Topw;
[Hws1]=Wit_Frc_calc_attributes (x1,N,winlength);
[ Avg1, Cv1, R1, theta1]=Wit_Sta_calc_attributes (x1,N,winlength);
[wwnew,wwnew1,wwnew2] = PCA_witness_attributes(Avg1,Cv1,R1,theta1);
Witstat=wwnew';
sigdis=0.152344;
winnum=winlength/sigdis; 
windis=sigdis;
winboundary=ceil(depth/windis);
wup=winboundary-(winnum/2);
wdown=winboundary+(winnum/2);
wstd(1,1)= Witstat(wup,1);
wstd(2,1)= Witstat(wdown,1);
wf(1,1)=Hws1(wup,1);
wf(2,1)=Hws1(wdown,1);
WFeatures=[wf,wstd];       % witness  features matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2st calculation observation
x2 = input('Enter the observation well log:');
Top=input(' enter depth of the top of the observation well (m): ');
realdepth=input('enter the real depth of selected boundary(m): ');
[Hws]= Fractal_calc_attributes(x2,N,winlength);
[ Avg,Cv,R, theta]=Statis_calc_attributes(x2,N,winlength);
[ownew] = PCA_observation_attributes(Avg,Cv,R,theta);
Obstat=ownew';
OFeatures=[Hws,Obstat];
StdFeatures=std(OFeatures);                                                    
Num=length(WFeatures);                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Proability calculation
win=floor((length(x2)-(winlength/sigdis))/(windis/sigdis)+1);
[Prob,Probability]=Prob_calc_attributes(x2,N,WFeatures,OFeatures,StdFeatures,Num,win,winlength,sigdis);
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
plot(Prob(1:end,2),Depth)                                          % PCA  statical  probability vs. depth
title('subplot 2: Statical Parameter proability')
%% outputs
aim=find(Probability==(max(max(Probability))));
aimdepth=Top+windis*(aim-1)+winlength;             % depth estimated
difdepth=aimdepth-realdepth;                                    % difference between estimated depth from the real depth
Real=ceil((realdepth-Top-winlength)/sigdis);
for j=Real-5:Real+5
Realprob(j,1:2)=Prob(j,1:2);                                            % feature probability of the real depth
PROB=max(Realprob);
end