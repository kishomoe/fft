clc
clear all;
format long;
Ns=1001;		%采样点
order=11;		%谐波数
filename='Bx.csv';		%导出的气隙磁密数据集
data=csvread(filename,1,0);     %从行偏移量 R1 和列偏移量 C1 开始读取文件中的数据
x=data(:,1);        %取第一列
y=data(:,2);        %取第二列
figure;
plot(x,y,'r');      %plot origional waveform
hold on;
grid on;
fft1=fft(y,Ns);
j=0;
amp_har=zeros(1,(order));
%%
for m=1:1:order
j=j+1;
fft1=fft(y,Ns);
fund_ele_front=fft1(m+1);
fund_ele_back=fft1(Ns+1-m);
amp_har(j)=(abs(fund_ele_front))/Ns*2;
fft1=0*fft1;
fft1(m+1)=fund_ele_front;
fft1(Ns+1-m)=fund_ele_back;
fft1=ifft(fft1,Ns);
fft1=real(fft1);
plot(x,fft1);
title('各次谐波分解示意图')
xlabel('长度(mm)')
ylabel('各次谐波幅值大小')
hold on;
end
%%
k=(1:1:order);
figure;
bar(k,amp_har);
title('各次谐波FFT分解柱状图')
xlabel('谐波次数')
ylabel('各次谐波幅值大小')
grid on;
peak_b=max(fft1);
rms_b=0.707*peak_b;
s=amp_har.^2;
THD = sqrt(sum(s(2:(order+1)/2)))/(amp_har(1));
THD_pct = 100*THD;