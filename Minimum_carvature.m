function [H,f]=Minimum_carvature()
load 'track.mat';
load 'body.mat';
[c,~]=size(track_x);
width = track_w / 2;  % һ��·��
%% ��ֵ
[insert_track_x,insert_track_y,insert_c]=insert_ref(track_x,track_y,c);  %���Բ�ֵ
[spline_track_x,spline_track_y,spline_c]=spline_ref(insert_track_x,insert_track_y,insert_c);%����������ֵ

%% ���ݲ�ֵ֮��ĵ㣬���¼���һ�׵���
[Ds_x,Ds_y,P_x,P_y]=D_track(spline_track_x,spline_track_y,spline_c);

%% ���㵥λ����
[vector]=get_vec(Ds_x,Ds_y,spline_c);

%% ����·������
x_b=zeros(spline_c,2);
y_b=zeros(spline_c,2);
for n=1:spline_c   
    x_b(n,1) = spline_track_x(n) + width* vector(n,2);%���߽�
    y_b(n,1) = spline_track_y(n) - width * vector(n,1);
    x_b(n,2) = spline_track_x(n)- width* vector(n,2);%�Ҳ�߽�
    y_b(n,2) = spline_track_y(n) + width * vector(n,1);
end

%% ����A���� 
[A]=buid_A(spline_c); 

%% ��� Mx  My
[Mx,My]=get_M(vector);

a=[0 0 1 0]; % 1��λ��ȡ������ȡ�ĸ�ϵ��
[A_exc]=get_Aexc(spline_c,a);
%% ����Pxx Pxy Pyy
[Pxx,Pxy,Pyy]=get_P(Ds_x,Ds_y,spline_c);

%% ����Tc Tnx Tny
[Tc,Tnx,Tny]=get_T(A_exc,A,Mx,My);

%% ���� Hx Hxy Hy  �� fx fxy fy
[H,f]=calculate_Hf(Tnx,Tny,Pxx,Pxy,Pyy,Tc,P_x,P_y);

%% ����ʽԼ��
[Ek,k]=inequation(Ds_y,Ds_x,Tc,Tny,Tnx,P_x,P_y,spline_c);

delta_k=1;
while delta_k>0.05
[alpha]=optimize(f,H,Ek,k,body,spline_c);
QP_track_x=track_x+vector(:,2).*alpha;
QP_track_y=track_y-vector(:,1).*alpha;
[k_last_max,loca_last]=Get_Max_K(track_x,track_y,spline_c,Tc);
[H,f,Ek,k]=upd_cons(QP_track_x,QP_track_y,spline_c,A_exc,A); %״̬����
[k_current_max,loca_current]=Get_Max_K(QP_track_x,QP_track_y,spline_c,Tc);
delta_k=abs(k_current_max-k_last_max);
track_x=QP_track_x;
track_y=QP_track_y;
plot(track_x,track_y);
hold on
end


%% ��ͼ
figure(1)
plot(spline_track_x,spline_track_y,'--r')
hold on
plot(x_b(:,1),y_b(:,1),'k')
hold on;
plot(x_b(:,2),y_b(:,2),'k')
hold on;
plot(QP_track_x,QP_track_y);


end



%% +++++++++++++++++++�Ӻ���++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [insert_track_x,insert_track_y,insert_c]=insert_ref(track_x,track_y,Dimension)
len=100;
insert_track_x=[];
insert_track_y=[];
for i=1:Dimension-1
    Distance=sqrt((track_x(i+1)-track_x(i))^2+(track_y(i+1)-track_y(i))^2);
    if Distance>len
        delta_x=track_x(i+1)-track_x(i);
        delta_y=track_y(i+1)-track_y(i);
        num=ceil(Distance/len); %��ɼ���
        temp_track_x=zeros(num,1);
        temp_track_y=zeros(num,1);
        
        for n=1:num
            add=(n-1)/num;
            temp_track_x(n,1)=track_x(i)+delta_x*add;
            temp_track_y(n,1)=track_y(i)+delta_y*add;
        end
        insert_track_x=[insert_track_x;temp_track_x];
        insert_track_y=[insert_track_y;temp_track_y];
    else
        insert_track_x=[insert_track_x;track_x(i)];
        insert_track_y=[insert_track_y;track_y(i)];
    end
end
insert_track_x=[insert_track_x;track_x(i+1)];
insert_track_y=[insert_track_y;track_y(i+1)];
[insert_c,~]=size(insert_track_x);
end

%% ����������ֵ
function [spline_track_x,spline_track_y,spline_c]=spline_ref(insert_track_x,insert_track_y,insert_c)
len=100;
spline_track_x=[];
spline_track_y=[];
[A]=buid_A(insert_c);
[P_x,P_y]=get_P_(insert_track_x,insert_track_y,insert_c);
para_x=A\P_x;
para_y=A\P_y;
loca=0;
for i=1:insert_c-1
    
    Distance=sqrt((insert_track_x(i+1)-insert_track_x(i))^2+(insert_track_y(i+1)-insert_track_y(i))^2);
    if Distance > len
        num=ceil(Distance/len);
        temp_x=zeros(num,1);
        temp_y=zeros(num,1);
        for n=1:num
            t=(1/num)*(n-1);
            temp_x(n,1)=para_x(loca+1)+para_x(loca+2)*t+para_x(loca+3)*t^2+para_x(loca+4)*t^3;
            temp_y(n,1)=para_y(loca+1)+para_y(loca+2)*t+para_y(loca+3)*t^2+para_x(loca+4)*t^3;
        end
        spline_track_x=[spline_track_x;temp_x];
        spline_track_y=[spline_track_y;temp_y];
    else
        spline_track_x=[spline_track_x;insert_track_x(i)];
        spline_track_y=[spline_track_y;insert_track_y(i)];
    end
    loca=loca+4;
end
spline_track_x=[spline_track_x;insert_track_x(i+1)];
spline_track_y=[spline_track_y;insert_track_y(i+1)];
[spline_c,~]=size(spline_track_x);
end

%% ����һ�׵���
function [Ds_x,Ds_y,P_x,P_y]=D_track(track_x,track_y,c)
a=[0 1 0 0];
[A_exc]=get_Aexc(c,a);

[A]=buid_A(c);
[P_x,P_y]=get_P_(track_x,track_y,c);
Ds_x=A_exc*(A \ P_x);
Ds_y=A_exc*(A \ P_y);
end

%% ���㵥λ����
function [vector]=get_vec(Ds_x,Ds_y,c)
   vector=zeros(c,2);
   for i=1:c
       mold=sqrt(Ds_x(i)^2+Ds_y(i)^2);
       vector(i,1)=Ds_x(i)/mold;
       vector(i,2)=Ds_y(i)/mold;
   end
end
 
 %% ����A����
 function [A]=buid_A(c)
A_cell=cell(c,c);
B_cell=cell(c-1,c-1);
C_cell=cell(c-1,1);
D_cell=cell(1,c);
a=[1 0 0 0;
   1 1 1 1;
   0 1 2 3;
   0 0 2 6];

b=[0 0 0 0;
   0 0 0 0;
   0 -1 0 0;
   0 0 -2 0];
for i=1:c %A matrix
    for j=1:c
        if i==j
            
            A_cell{i,j}=a;
        else
             A_cell{i,j}=zeros(4,4);
         end
    end
end
A=cell2mat(A_cell);

for i=1:c-1 % C matrix
    C_cell{i,1}=zeros(4,4);
end

for i=1:c % D matrix
    if i==1
        D_cell{1,i}=b;
    else
        D_cell{1,i}=zeros(4,4);
    end
end  

for i=1:c-1  %B matrix
    for j=1:c-1
        if i==j
            B_cell{i,j}=b;
        else
             B_cell{i,j}=zeros(4,4);
         end
     end
end

B_cell=[C_cell,B_cell;D_cell];
B=cell2mat(B_cell);
A=A+B;
 end
 %% ��þ���P_x P_y
 function[P_x,P_y]=get_P_(track_x,track_y,c)
 P_x_cell=cell(c,1);
 
for a=1:c
    if a==c
     p_x=[track_x(a); track_x(1); 0; 0];
   
    else
     p_x=[track_x(a); track_x(a+1); 0; 0];
    end
    P_x_cell{a,1}=p_x;
end
 P_x=cell2mat(P_x_cell);


P_y_cell=cell(c,1);
for n=1:c
    if n==c
     p_y=[track_y(n); track_y(1); 0; 0];
   
    else
     p_y=[track_y(n); track_y(n+1); 0; 0];
    end
    P_y_cell{n,1}=p_y;
end
 P_y=cell2mat(P_y_cell);
 end

 %%  ���M����
 function [Tc,Tnx,Tny]=get_T(A_exc,A,Mx,My)
Tc=2*A_exc/A;
Tnx=(2*A_exc/A)*Mx;
Tny=(2*A_exc/A)*My;
end

%% ��ȡ Mx My ����
function [Mx,My]=get_M(vector)
[c,~]=size(vector);
        Mx_cell=cell(c,1);
         
for i=1:c 
    A=zeros(4,c);%����Mx����
    if i<c-1
        A(1,i)=vector(i,2);
        A(2,i+1)=vector(i+1,2);
    else
        A(2,1)=vector(1,2);
        A(1,i)=vector(i,2);
    end
    Mx_cell{i,1}=A;
end
Mx=cell2mat(Mx_cell);

My_cell=cell(c,1); %����My����
for i=1:c
    A=zeros(4,c);
    if i<c-1
        A(1,i)=-vector(i,1);
        A(2,i+1)=-vector(i+1,1);
    else
        A(2,1)=-vector(1,1);
        A(1,i)=-vector(i,1);
    end
    My_cell{i,1}=A;
end
My=cell2mat(My_cell);
end

%% ��ȡ����
function [A_exc]=get_Aexc(c,a)
A_exc_cell=cell(c,c);

for i=1:c
    for j=1:c
        if i==j
            A_exc_cell{i,j}=a;
        else
            A_exc_cell{i,j}=zeros(1,4);
        end
    end
end

A_exc=cell2mat(A_exc_cell);
end

%% P����
function [Pxx,Pxy,Pyy]=get_P(Ds_x,Ds_y,c)
Pxx=zeros(c,c);
Pxy=zeros(c,c);
Pyy=zeros(c,c);


for i=1:c
    for j=1:c
        if i==j
            Pxx(i,j)=Ds_y(i)^2/(Ds_x(i)^2+Ds_y(i)^2)^3;
            Pxy(i,j)=-2*Ds_x(i)*Ds_y(i)/(Ds_x(i)^2+Ds_y(i)^2)^3;
            Pyy(i,j)=Ds_x(i)^2/(Ds_x(i)^2+Ds_y(i)^2)^3;
        else
            Pxx(i,j)=0;
            Pxy(i,j)=0;
            Pyy(i,j)=0;
        end
    end
end
end


 %% H f
function [H,f]=calculate_Hf(Tnx,Tny,Pxx,Pxy,Pyy,Tc,P_x,P_y)
Hx=Tnx'*Pxx*Tnx;
Hxy=Tny'*Pxy*Tnx;
Hy=Tny'*Pyy*Tny;

fx=2*Tnx'*Pxx'*Tc*P_x;
fxy=Tny'*Pxy'*Tc*P_x+Tnx'*Pxy'*Tc*P_y;
fy=2*Tny'*Pyy'*Tc*P_y;

H=(Hx+Hxy+Hy)+(Hx+Hxy+Hy)';
f=(fx+fxy+fy)';
end
 
%% 
function [Q_y,Q_x]=get_Qy(Ds_y,Ds_x,c)
Q_x=zeros(c,c);
Q_y=zeros(c,c);

for a=1:c
    for b=1:c
        if a==b
            Q_x(a,b)=Ds_y(a)/(Ds_x(a)^2+Ds_y(a)^2)^(3/2);
            Q_y(a,b)=Ds_x(a)/(Ds_x(a)^2+Ds_y(a)^2)^(3/2);
        end
    end
end
end

 %% ����ʽԼ��
function [Ek,k]=inequation(Ds_y,Ds_x,Tc,Tny,Tnx,P_x,P_y,c)
const=ones(c,1);
K_bound=0.12*const;
[Q_y,Q_x]=get_Qy(Ds_y,Ds_x,c);
Kref=Q_y*Tc*P_y-Q_x*Tc*P_x;
K_upper_bound=K_bound-Kref;
K_lower_bound=K_bound+Kref;

%Az<=b �� A=Ek,b=k;
k=[K_upper_bound;K_lower_bound];
Ek=[Q_y*Tny-Q_y*Tnx;-Q_y*Tny+Q_y*Tnx];
end

%% ״̬����
function [H,f,Ek,k]=upd_cons(track_x,track_y,c,A_exc,A)
[Ds_x,Ds_y,P_x,P_y]=D_track(track_x,track_y,c);
[vector]=get_vec(Ds_x,Ds_y,c);
[Mx,My]=get_M(vector);
[Pxx,Pxy,Pyy]=get_P(Ds_x,Ds_y,c);
[Tc,Tnx,Tny]=get_T(A_exc,A,Mx,My);
[H,f]=calculate_Hf(Tnx,Tny,Pxx,Pxy,Pyy,Tc,P_x,P_y);
[Ek,k]=inequation(Ds_y,Ds_x,Tc,Tny,Tnx,P_x,P_y,c);
end

%% �Ż�
 function [alpha]=optimize(f,H,Ek,k,body,c)
 lb=(-3+body.Wb/2000)*ones(c,1);
 ub=(3-body.Wb/2000)*ones(c,1);
 alpha=quadprog(H,f,Ek,k,[],[],lb,ub);
 end
 
 
 %% �����������
function [K_QP_MAX,location]=Get_Max_K(track_x,track_y,c,Tc)  
[Ds_x,Ds_y,P_x,P_y]=D_track(track_x,track_y,c);
[Q_y,Q_x]=get_Qy(Ds_y,Ds_x,c);
k_QP=Q_y*Tc*P_y-Q_x*Tc*P_x;
[K_QP_MAX,location]=max(abs(k_QP));
end