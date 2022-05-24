function [Val_min,Val_moy,Val_max] = calcul_int(y,percent)

for k=1:size(y,2)
    if min(y(:,k))==max(y(:,k))
        Val_moy(1,k)=y(1,k);
        Val_min(1,k)=y(1,k);
        Val_max(1,k)=y(1,k);
        
    else
        F=[1/size(y,1):(1-2/size(y,1))/(size(y,1)-1):1-1/size(y,1)]';
        ys=sort(y(:,k));
        Val_moy(1,k) = mean(y);%interp1(F,ys,0.5);%
        Val_min(1,k) = interp1(F,ys,percent/2);
        Val_max(1,k) = interp1(F,ys,1-percent/2);
    end
end
if nargout<3
    figure;plot(Val_moy,'*')
    hold on;plot(Val_min,'--r')
    plot(Val_max,'--r')
end