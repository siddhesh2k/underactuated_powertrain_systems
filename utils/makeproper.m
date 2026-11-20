function Gpp=makeproper(Gp,tau)
Gp=tf(Gp);
s=tf('s');
tg=isproper(Gp);
if tg==1
    Gpp=Gp;
else
    np=size(pole(Gp));
    np=np(1,1);
    nz=size(zero(Gp));
    nz=nz(1,1);
    nr=nz-np;
    Gpp=minreal(Gp*1/(tau*s+1)^(nr),1e-3);
end
