function [f3]=buildF2Fig(configs,U_data)

axisSize=configs.axisSize;
w=configs.w;
N=configs.N;
x=U_data(1:N,2); y=x;

ucaxisMax=ceil(max(U_data(:,5)));
F2labels=["F2aaaa", "F2aaab", "F2aabb", "F2abaa", "F2abab", "F2abbb", "F2bbaa", "F2bbab", "F2bbbb"];
f3=figure('Name','F2 Simulations');
f3.Position=[1201 925 1000 1000];
pos=zeros(9,4);
pos(:,3)=0.27;
pos(:,4)=0.25;
    
for i=1:3
    pos(3*i,1)=0.69;
    pos(3*i-1,1)=0.37;
    pos(3*i-2,1)=0.05;

    pos(i,2)=0.69;
    pos(i+3,2)=0.37;
    pos(i+6,2)=0.05;
end


figure(f3)
for i=1:9
    subplot('Position',pos(i,:))
    U_F2=reshape(U_data(:,i+5),N,N);
    surf(x,y,U_F2,"EdgeColor","none")
    view(2)
    colorbar
    clim([0 ucaxisMax])
    axis([0 axisSize*w 0 axisSize*w])
    title(F2labels(i))
end
    
    
    