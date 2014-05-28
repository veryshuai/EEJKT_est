function [yy0,yy1,yy2,sp_h,sp_f,sp_p_uni,ud_p,sp_d,Phi_init,drw1,drw2,drw3,drw4,drw5,drw6,drw7,drw8,drw9,drw10,drw11,drw12,drw13,drw14,drw15] = ...
    getshks(j)

persistent shk

if isempty(shk) == 1
    load data/srates.mat;
    load data/eps.mat;
    load data/drw.mat;
    shk.yy0 = yy0; shk.yy1=yy1; shk.yy2=yy2; shk.sp_h=sp_h; shk.sp_f=sp_f; shk.sp_p_uni=sp_p_uni; shk.ud_p=ud_p;
    shk.drw1=drw1; shk.drw2=drw2; shk.drw3=drw3; shk.drw4=drw4; shk.drw5=drw5; shk.drw6=drw6; shk.drw7=drw7; shk.drw8=drw8; shk.drw9=drw9; 
    shk.drw10=drw10; shk.drw11=drw11; shk.drw12=drw12; shk.drw13=drw13; shk.drw14 = drw14; shk.drw15 = drw15; shk.Phi_init = Phi_init;
    shk.sp_d = sp_d;
end    

if j == 0
    yy0 = shk.yy0; yy1=shk.yy1; yy2=shk.yy2; sp_h=shk.sp_h; sp_f=shk.sp_f; sp_p_uni=shk.sp_p_uni; ud_p=shk.ud_p;
    drw1=shk.drw1; drw2=shk.drw2; drw3=shk.drw3; drw4=shk.drw4; drw5=shk.drw5; drw6=shk.drw6; drw7=shk.drw7; drw8=shk.drw8; drw9=shk.drw9; 
    drw10=shk.drw10; drw11=shk.drw11; drw12=shk.drw12; drw13=shk.drw13; drw14 = shk.drw14; drw15 = shk.drw15; Phi_init = shk.Phi_init;
    sp_d = shk.sp_d;
else
    yy0 = shk.yy0; yy1=shk.yy1; yy2=shk.yy2; sp_h=shk.sp_h; sp_f=shk.sp_f; sp_p_uni=shk.sp_p_uni; ud_p=shk.ud_p;
    drw1=shk.drw1(j,:); drw2=shk.drw2(j,:); drw3=shk.drw3(j,:); drw4=shk.drw4(j,:); drw5=shk.drw5(j,:); drw6=shk.drw6(j,:); drw7=shk.drw7(j,:); drw8=shk.drw8(j,:); drw9=shk.drw9(j,:); 
    drw10=shk.drw10(j,:); drw11=shk.drw11(j,:); drw12=shk.drw12(j,:); drw13=shk.drw13(j,:); drw14 = shk.drw14(j,:); drw15 = shk.drw15; Phi_init = shk.Phi_init;
    sp_d = shk.sp_d;
end

end
