 %% The script is used to compare the floating implementation of algorithm with MATLAB and the fixed-point implementation of ASIC architecture. We use a real dataset names R_Awaken,
 %%% to perform this test (see README file for the content of the dataset) extracting a frame of 1s. In article we performed several times these steps varying SNR.
 
 
 
 x=load('DatasetR_Awaken.mat');
 sf=20e+3;
 sampleframe=1*sf;
 Raw=x.AwR(:,1:sampleframe);
 

 windowL=64;
 %% Double-precision ASO-WA
load('ButtII20kHz.mat') %Q1.8 to Q9.0
a=Ainteger./2^8;
b=Binteger./2^8;
[Bf,Af]=butter(1,[300 3000]./(sf/2),'bandpass');
EBPF=[Bf,Af]-[b,a];
 RMSE_BPF=sqrt(immse([Bf,Af],[b,a]));

 filtered_data=filter(b,a,Raw,[],2);
 time=size(Raw,2);
yt=zeros(size(Raw,1),size(Raw,2)+2);

 for t=3:time
        
yt(:,t)=-a(2)*yt(:,t-1)-a(3)*yt(:,t-2)+...
    +b(1)*Raw(:,t)+b(3)*Raw(:,t-2);
yt(:,t)=(yt(:,t));
end

Filtered_Double=yt(:,3:end);


%%Spatial mean 4
 
MEAN_DOUBLE=(mean(Filtered_Double,1));

%%kASO
k=4;
load('hamminginteger17.mat')
yaso=[];
yk=aso(MEAN_DOUBLE,k);
YASO_DOUBLE=yk(k+1:end-k);

s=[];

for i=1:(size(YASO_DOUBLE,2)-17)
S_DOUBLE(i)=dot(YASO_DOUBLE(i:i+16),Binteger./2^9);   
end

%%BATCH
n=size(MEAN_DOUBLE,2);
nbatch=ceil(size(MEAN_DOUBLE,2)/windowL);
sigma0=((sum(abs(MEAN_DOUBLE(1:windowL)))*1.25));

for k=1:nbatch
                jstrt=(k-1)*windowL+1;
                if k ~= nbatch
                    jend=jstrt+windowL-1;
                else
                    jend=n;
                end
                %AA
                y_AA_tmp(jstrt:jend)=((sum(abs(MEAN_DOUBLE(jstrt:jend)))*1.25));
                y_MAD_tmp(jstrt:jend)=(((median(abs(MEAN_DOUBLE(jstrt:jend))./0.6745))));
                if k>1
                signal_1=abs(MEAN_DOUBLE(jstrt:jend));
                signal_2=min( signal_1, sigma0);
                sigma_2a=floor((sum(signal_2)*1.58));
                sigma0=sigma_2a;
                end
                y_WA_tmp(jstrt:jend) = sigma0;

end
y_AA_DOUBLE=y_AA_tmp;
y_MAD_DOUBLE=y_MAD_tmp;
y_WA_DOUBLE=y_WA_tmp;
var_AA_DOUBLE=y_AA_DOUBLE.^2.*32;
var_MAD_DOUBLE=(y_MAD_DOUBLE.^2).*32;
var_WA_DOUBLE=(y_WA_DOUBLE.^2).*32;
%%Thresholding
latency_block=43:size(S_DOUBLE,2); %da aggiungere i 4 colpi per l'aso e 17 dello smothing
AP_BIT_AA_DOUBLE=S_DOUBLE(latency_block)>var_AA_DOUBLE(latency_block);
AP_BIT_WA_DOUBLE=S_DOUBLE(latency_block)>var_WA_DOUBLE(latency_block);
AP_BIT_MAD_DOUBLE=S_DOUBLE(latency_block)>var_MAD_DOUBLE(latency_block);

 
 
 


 %% Fixed-point ASO-WA
quantizedRaw=num2fixpt(Raw,sfix(10),2^0); %Quantization process

%%Butterworth Filter BPF II order 
load('ButtII20kHz.mat') %Q1.8 to Q9.0
A=Ainteger;
B=Binteger;

yt=zeros(size(Raw,1),size(Raw,2)+2);
time=size(Raw,2);
row_n=4;

for t=3:time
        
yt(:,t)=-floor(A(2)*yt(:,t-1))-floor(A(3)*yt(:,t-2))+...
    +floor(B(1)*quantizedRaw(:,t))+floor(B(3)*quantizedRaw(:,t-2));
yt(:,t)=floor((yt(:,t)/2^8));
end

st=yt(:,3:end);


%%Spatial mean 4
 
ymean=floor(mean(st,1));

%%kASO
k=4;

yaso=[];
yk=aso(ymean,k);
yaso=yk(k+1:end-k);
load('hamminginteger17.mat')
b=Binteger;
a=Ainteger;
s=[];
EHAMM=hamming(17).'-b./2^9;
 RMSE_FIR=sqrt(immse(hamming(17).',b./2^9));

for i=1:(size(yaso,2)-17)
s(i)=floor(dot(yaso(i:i+16),b)./2^9);   
end

%%BATCH
n=size(ymean,2);
nbatch=ceil(size(ymean,2)/windowL);
sigma0=floor((sum(abs(ymean(1:windowL)))*5)/2^8);
% y_AA=[zeros(size(ymean)) zeros(1,window)];
% y_MAD=[zeros(size(ymean)) zeros(1,window)];
% y_WA=[zeros(size(ymean)) zeros(1,window)];
for k=1:nbatch
                jstrt=(k-1)*windowL+1;
                if k ~= nbatch
                    jend=jstrt+windowL-1;
                else
                    jend=n;
                end
                %AA
                y_AA_tmp(jstrt:jend)=floor((sum(abs(ymean(jstrt:jend)))*5)/2^8);
                y_MAD_tmp(jstrt:jend)=floor((floor(median(abs(ymean(jstrt:jend))))*11)/16);
                if k>1
                signal_1=abs(ymean(jstrt:jend));
                signal_2=min( signal_1, sigma0);
                sigma_2a=floor((sum(signal_2)*6)./2^8);
                sigma0=sigma_2a;
                end
                y_WA_tmp(jstrt:jend) = sigma0;

end
y_AA=y_AA_tmp;
y_MAD=y_MAD_tmp;
y_WA=y_WA_tmp;
var_AA=y_AA.^2.*32;
var_MAD=(y_MAD.^2).*32;
var_WA=(y_WA.^2).*32;
%%Thresholding
latency_block=43:size(s,2); %da aggiungere i 4 colpi per l'aso e 17 dello smothing
AP_BIT_AA=s(latency_block)>var_AA(latency_block);
AP_BIT_WA=s(latency_block)>var_WA(latency_block);
AP_BIT_MAD=s(latency_block)>var_MAD(latency_block);

%% Parameters detection
number_spike_D=nnz(diff(x.Pattern_TP(latency_block))>0);
detect_spike_D=(AP_BIT_WA_DOUBLE.*x.Pattern_TP(latency_block));
fault_spike_D=(AP_BIT_WA_DOUBLE.*x.Pattern_FP(latency_block));
TP_D=nnz(diff(detect_spike_D)>0);



FP_D=nnz(diff(fault_spike_D)>0);
FN_D=number_spike_D-TP_D;
Accuracy_D=(TP_D/(TP_D+FN_D+FP_D));


TPR_D=TP_D/(TP_D+FN_D);
FAR_D=FP_D/(TP_D+FP_D);




number_spike_F=nnz(diff(x.Pattern_TP(latency_block))>0);
detect_spike_F=(AP_BIT_WA.*x.Pattern_TP(latency_block));
fault_spike_F=(AP_BIT_WA.*x.Pattern_FP(latency_block));
TP_F=nnz(diff(detect_spike_F)>0);



FP_F=nnz(diff(fault_spike_F)>0);
FN_F=number_spike_F-TP_F;
Accuracy_F=(TP_F/(TP_F+FN_F+FP_F));


TPR_F=TP_F/(TP_F+FN_F);
FAR_F=FP_F/(TP_F+FP_F);

%% Error Measure
 E=S_DOUBLE-s;
 ED=abs(E);
 MED=mean(ED);
 ErrorCnt=nnz(E);
 RMSE_ACC=sqrt(immse(Accuracy_D,Accuracy_F));
 
 RMSE=sqrt(immse(S_DOUBLE,s));
 
 E=AP_BIT_WA_DOUBLE-AP_BIT_AA;
 
 ED_BPF=abs(EBPF);
 ED_BPF=abs(EHAMM);
 
fprintf('The Root Mean Square Error Accuracy is 2%f2.\n',RMSE_ACC);
fprintf('The Mean Error Distance up to Smoothing is 2%f2.\n',MED);
fprintf('The Root Mean Square Error up to Smoothing is 2%f2.\n',RMSE);
fprintf('The RMSE of the BPF is 2%f2.\n',RMSE_BPF);
fprintf('The RMSE of the Hamming is 2%f2.\n',RMSE_FIR);
