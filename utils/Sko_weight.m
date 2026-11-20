function wp=Sko_weight(M,wb,A,n)
m=1/(M^(1/n));
a=A^(1/n);
den=[];
num=[];
for i=0:n
        C=nchoosek(n,i);
        numi=C*(m^(n-i))*(wb^(i));
        deni=C*(wb*a)^(i);
        num=[num,numi];
        den=[den,deni];
end
wp=tf(num,den);