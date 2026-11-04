function [verim, gviol, out] = im_core(x)
%IM_CORE  Analytical design model of the induction motor.
%   D=16: Lamda_katsayi and KD_katsayi are now variables. Dout=0.180.
%   cos_phi is computed from the equivalent circuit. 10 constraints (below).
%
%   x ORDER:
%     1 Lamda_katsayi 2 KD_katsayi 3 Bg_katsayi 4 Kst 5 Jcos_katsayi
%     6 bos[mm] 7 hos[mm] 8 hw[mm] 9 Bts[T] 10 Jb[A/mm2] 11 Btr[T]
%     12 hor[mm] 13 bor[mm] 14 rotor_slot[deg] 15 Bcr[T] 16 k
%
%   CONSTRAINTS (gviol, each element >=0; 0 = satisfied):
%     g1 Bcs<1.7 | g2 Dshaft_in>0 | g3 Sn in[0.005,0.05] | g4 tbk in[1.75,3]
%     g5 tLR in[0.70,2.75] | g6 iLR in[6,7] | g7 Fmts<Fmtr1<=1.1*Fmts (band)
%     g8 Fmtr1~=Fmtr2 | g9 Kst2~=Kst | g10 cos_phi_1n~=cos_phi

verim = 0; out = struct();
NC = 10; gviol = 1e6 * ones(1, NC);

try
    % ---- Fixed design inputs ----
    Pn = 5500; V1_line = 400; V1_ph = V1_line/sqrt(3);
    f1 = 50; m = 3; p1_cift = 2; p1_tek = 2*p1_cift;
    pf = power_factor(p1_cift,Pn);
    eff = 0.01*efficiency_lookup("IE3",p1_cift,Pn);

    % ---- Optimization variables (x) ----
    Lamda_katsayi = x(1); KD_katsayi = x(2);
    Bg_katsayi = x(3);    Kst = x(4);  Jcos_katsayi = x(5);
    bos = x(6)/1000; hos = x(7)/1000; hw = x(8)/1000;
    Bts = x(9); Jb = x(10); Btr = x(11);
    hor = x(12)/1000; bor = x(13)/1000; rotor_slot_derece = x(14);
    Bcr = x(15); k = x(16);

    % ---- Fixed discrete choices ----
    Nr=28; y=7; q1=3; a1=1; ap=1;

    %%  MAIN DIMENSIONS
    KE = 0.98-0.005*p1_cift;
    Sgap = KE*Pn/eff/pf;
    if p1_tek == 2,     Lamda = Lamda_katsayi*(0.6+1)/2;
    elseif p1_tek == 4, Lamda = Lamda_katsayi*(1.2+1.8)/2;
    elseif p1_tek == 6, Lamda = Lamda_katsayi*(1.6+2.2)/2;
    elseif p1_tek == 8, Lamda = Lamda_katsayi*(2+3)/2; end
    Co_katsayi = 1;
    Co = Co_katsayi*esson_constant(p1_cift,Pn);
    Dis1 = (p1_tek/pi/Lamda*p1_cift/f1*Sgap/Co)^(1/3);
    L1 = Lamda*pi*Dis1/p1_tek;
    Dout = 0.180;
    if p1_tek == 2,     KD = KD_katsayi*(0.54+0.58)/2;
    elseif p1_tek == 4, KD = KD_katsayi*(0.61+0.63)/2;
    elseif p1_tek == 6, KD = KD_katsayi*(0.68+0.71)/2;
    elseif p1_tek == 8, KD = KD_katsayi*(0.72+0.74)/2; end
    Dis2 = Dout*KD; L2 = Dis1^2*L1/Dis2^2; Dis3 = Dis2;
    L3 = ceil(1000*L2)/1000; Tau3 = pi*Dis3/p1_tek;
    q3 = m*q1; Tau_s = Tau3/q3;
    g = 0.00035;

    %%  STATOR WINDING
    Ns = p1_tek*q1*m;
    Alfa_ec_rad = 2*pi*p1_cift/Ns;
    Kq1 = sin(q1*Alfa_ec_rad/2)/q1/sin(Alfa_ec_rad/2);
    Ky1 = sin(pi/2*y/q3); Kw1 = Kq1*Ky1;
    if p1_tek == 2,     Bg1 = ((Dis1^2*L1)/(Dis3^2*L3))*Bg_katsayi*(0.5+0.75)/2;
    elseif p1_tek == 4, Bg1 = ((Dis1^2*L1)/(Dis3^2*L3))*Bg_katsayi*(0.65+0.78)/2;
    elseif p1_tek == 6, Bg1 = ((Dis1^2*L1)/(Dis3^2*L3))*Bg_katsayi*(0.7+0.82)/2;
    elseif p1_tek == 8, Bg1 = ((Dis1^2*L1)/(Dis3^2*L3))*Bg_katsayi*(0.75+0.85)/2; end
    Kf = Kst_Kf(Kst); Alfa_i = Kst_Alfa_i(Kst);
    Phi = Alfa_i*Tau3*L3*Bg1;
    W1 = KE*V1_ph/4/Kf/Kw1/f1/Phi;
    if mod(a1*(W1-mod(W1,p1_cift*q1))/p1_cift/q1,2) == 0
        ns = a1*(W1-mod(W1,p1_cift*q1))/p1_cift/q1;
    else
        ns = a1*(W1-mod(W1,p1_cift*q1))/p1_cift/q1-1;
    end
    W11 = p1_cift*q1*ns; Bg2 = Bg1*W1/W11; Phi = Alfa_i*Tau3*L3*Bg2;
    I1n = Pn/eff/pf/sqrt(3)/V1_line;
    if p1_tek == 2 || p1_tek == 4, Jcos = Jcos_katsayi*(4+7)/2;
    else, Jcos = Jcos_katsayi*(5+8)/2; end
    A_Co = I1n/Jcos/a1; d_Co = conductor_diameter(sqrt(4*A_Co/pi/ap));

    %%  STATOR SLOT
    if Pn < 10000, Kfill = (0.4/0.375)*(0.35+0.4)/2; else, Kfill = (0.4/0.375)*(0.4+0.44)/2; end
    A_su = pi*d_Co^2*ap*ns/4/Kfill; Kfe = 0.96;
    bts = Bg2*Tau_s/Bts/Kfe; bs1 = pi*(Dis3+2*hos+2*hw)/Ns-bts;
    bs2 = sqrt(bs1^2+4*(10^-6*A_su)*tan(pi/Ns)); hs = 2*(10^-6*A_su)/(bs1+bs2);
    Kc = 1.15; mu_o = 4*pi*10^-7; Fmg = Kc*g*Bg2/mu_o;
    Hts = BmHm(Bts); Fmts = Hts*(hs+hos+hw);
    Fmtr1 = Kst*Fmg-Fmts;                        % related to constraints g7/g8
    hcs = (Dout-(Dis3+2*(hos+hw+hs)))/2;
    Bcs = Phi/2/L3/hcs;                          % constraint g1: Bcs<1.7

    %%  ROTOR SLOT
    KI = (0.8*pf+0.2); Ib = KI*2*m*W11*Kw1*I1n/Nr; Ab = Ib/(Jb*10^6);
    Ier = Ib/2/sin(pi*p1_cift/Nr); Jer = (0.75/0.775)*Jb*(0.75+0.8)/2; Aer = Ier/(Jer*10^6);
    Tau_r = pi*(Dis3-2*g)/Nr; btr = Bg2*Tau_r/Kfe/Btr; Dre = Dis3-2*g;
    d1 = (pi*(Dre-2*hor)-Nr*btr)/(Nr+pi);
    rotor_slot_rad = deg2rad(rotor_slot_derece);
    d2 = sqrt((Ab-(rotor_slot_rad+1/(4*tan(pi/Nr)))*d1^2)/(rotor_slot_rad-1/(4*tan(pi/Nr))));
    hr = (d1-d2)/2/tan(pi/Nr);
    Htr = BmHm(Btr); Fmtr2 = Htr*(hr+hor+(d1+d2)/2);   % constraint g8
    hcr = Phi/2/L3/Bcr;
    Dshaft_in = Dis3-2*g-2*(hor+(d1+d2)/2+hr+hcr);      % constraint g2
    Der = Dre-2*0.003; b = k*(hr+hor+(d1+d2)/2); a = Aer/b;

    %%  MAGNETIZATION
    Gamma1 = bos^2/(5*g+bos); Gamma2 = bor^2/(5*g+bor);
    Kc1 = Tau_s/(Tau_s-Gamma1); Kc2 = Tau_r/(Tau_r-Gamma2); Kc = Kc1*Kc2;
    Hcs = BmHm(Bcs); Hcr = BmHm(Bcr);
    Ccs = 0.88*exp(-0.4*Bcs^2); Ccr = 0.88*exp(-0.4*Bcr^2);
    Fmcs = Ccs*pi*(Dout-hcs)/p1_tek*Hcs; Fmcr = Ccr*pi*(Dshaft_in+hcr)/p1_tek*Hcr;
    F1m = 2*(Fmg+Fmts+Fmtr2+Fmcs+Fmcr); Ks = F1m/2/Fmg-1;
    Imu = pi*p1_cift*(F1m/2)/(3*sqrt(2)*W11*Kw1);
    Kst2 = (Fmts+Fmtr2)/Fmg;                            % constraint g9

    %%  RESISTANCES AND REACTANCES
    Beta = y/q3; y_Beta = Beta*Tau3;
    if p1_tek == 2,     l_end = 2*y_Beta-0.04;
    elseif p1_tek == 4, l_end = 2*y_Beta-0.02;
    elseif p1_tek == 6, l_end = pi/2*y_Beta+0.018;
    elseif p1_tek == 8, l_end = 2.2*y_Beta-0.012; end
    l_c = 2*(L3+l_end);
    Ro_Co_20 = 1.78*10^-8; T = 75; Ro_Co_80 = Ro_Co_20*(1+(T-20)/273);
    Rs = Ro_Co_80*l_c*W11/a1/(A_Co*10^-6);
    l_er = pi*(Der-b)/Nr;
    Ro_Al_20 = 3.1*10^-8; Ro_Al_80 = Ro_Al_20*(1+(T-20)/273);
    Omega1 = 2*pi*f1; Beta_s = sqrt(Omega1*mu_o/2/Ro_Al_20);
    S = 1; Xi = Beta_s*hr*sqrt(S);
    KR1 = Xi*((sinh(2*Xi)+sin(2*Xi))/(cosh(2*Xi)-cos(2*Xi)));
    Rbe_80_S1 = Ro_Al_80*(L3/Ab*KR1+l_er/2/Aer/(sin(pi*p1_cift/Nr))^2);
    Rr_S1 = 4*m/Nr*(W11*Kw1)^2*Rbe_80_S1;
    Delta_s = (2/3*hs/(bs1+bs2)+2*hw/(bos+bs1)+hos/bos)*((1+3*Beta)/4);
    Cs = 1-0.033*bos^2/g/Tau_s; Phi1 = pi*(6*Beta-5.5);
    Gamma_ds = (0.18*sin(Phi1)+1.24)*0.01;
    Delta_ds = 0.9*Tau_s*q1^2*Kw1^2*Cs*Gamma_ds/Kc/g/(1+Kst);
    Delta_ec = 0.34*q1/L3*(l_end-0.64*Beta*Tau3);
    Xsl = 2*mu_o*Omega1*L3*W11^2/p1_cift/q1*(Delta_s+Delta_ds+Delta_ec);
    Delta_r = 0.66+2*hr/3/(d1+d2)+hor/bor;
    Gamma_dr = 9*(6*p1_cift/Nr)^2*10^-2;
    Delta_dr = 0.9*Tau_r*Gamma_dr/Kc/g*(Nr/6/p1_cift)^2;
    Delta_er = 2.3*(Der-b)/Nr/L3/4/(sin(pi*p1_cift/Nr))^2*log10(4.7*(Der-b)/(b+2*a));
    Kx1 = 3/2/Xi*(sinh(2*Xi)-sin(2*Xi))/(cosh(2*Xi)-cos(2*Xi));
    Xbe_S1 = Omega1*mu_o*L3*(Delta_r*Kx1+Delta_dr+Delta_er);
    Xrl_S1 = 4*m*(W11*Kw1)^2/Nr*Xbe_S1;
    KR2 = 1; Rbe_80_Sn = Ro_Al_80*(L3/Ab*KR2+l_er/2/Aer/(sin(pi*p1_cift/Nr)^2));
    Rr_Sn = Rr_S1*Rbe_80_Sn/Rbe_80_S1;
    Kx2 = KR2; Xbe_Sn = Omega1*mu_o*L3*(Delta_r*Kx2+Delta_dr+Delta_er);
    Xrl_Sn = 4*m*(W11*Kw1)^2/Nr*Xbe_Sn;
    Xsl_sat_S1 = (0.7+0.8)/2*Xsl; Xrl_sat_S1 = (0.6+0.7)/2*Xrl_Sn;
    Xm = sqrt((V1_ph/Imu)^2-Rs^2)-Xsl;
    c = 1; Kskew = sin(pi/2*c*Tau_s/Tau3)/(pi/2*c*Tau_s/Tau3);
    Xm_skew = Kskew*Xm; Xrl_skew = Xm*(1-Kskew^2);
    Xrl_skew_S1 = Xrl_sat_S1+Xrl_skew; Xrl_skew_Sn = Xrl_S1+Xrl_skew;

    %%  LOSSES AND EFFICIENCY
    P_Co = 3*Rs*I1n^2; I_rn = KI*I1n; P_Al = 3*Rr_Sn*I_rn^2;
    P_mv = 0.012*Pn; P_stray = 0.01*Pn;
    p10 = 1.38; Kt = (1.6+1.8)/2; Gi = 7800;
    G_tl = Gi*Ns*bts*(hs+hw+hos)*L3*Kfe; P_tl = Kt*p10*(f1/50)^1.3*Bts^1.7*G_tl;
    Ky = (1.6+1.9)/2; G_ts = Gi*pi/4*(Dout^2-(Dout-2*hcs)^2)*L3*Kfe;
    P_yl = Ky*p10*(f1/50)^1.3*Bcs^1.7*G_ts;
    Kps = 1/(2.2-Bts); Kpr = 1/(2.2-Btr); Bps = (Kc1-1)*Bg2; Bpr = (Kc2-1)*Bg2;
    G_tr = Gi*L3*Kfe*Nr*(hr+(d1+d2)/2)*btr;
    Ps_iron = 0.5*10^-4*((Nr*f1/p1_cift*Kps*Bps)^2*G_ts+(Ns*f1/p1_cift*Kpr*Bpr)^2*G_tr);
    P_iron = P_tl+P_yl+Ps_iron;
    Rc = (0.97*V1_ph)^2/(P_iron/3);                     % for cos_phi
    Losses = P_Co+P_Al+P_iron+P_mv+P_stray;
    eff_h = Pn/(Pn+Losses);

    %%  OPERATING CHARACTERISTICS
    Sn = P_Al/(Pn+P_Al+P_mv+P_stray);                   % constraint g3
    Tn = Pn/(2*pi*f1/p1_cift*(1-Sn));
    Cm = 1+Xsl/Xm_skew;
    Tbk = 3*p1_cift/2/Omega1*V1_ph^2/(Rs+sqrt(Rs^2+(Xsl+Cm*Xrl_skew_Sn)^2));
    ILR = V1_ph/sqrt((Rs+Rr_S1)^2+(Xsl_sat_S1+Xrl_skew_S1)^2);
    TLR = 3*Rr_S1*ILR^2/Omega1*p1_cift;
    tbk = Tbk/Tn; tLR = TLR/Tn; iLR = ILR/I1n;          % constraints g4,g5,g6
    eff_h = 100*eff_h;
    cos_phi_1n = Pn/(3*V1_ph*I1n*(eff_h/100));          % constraint g10
    cos_phi = cos(angle(Rs+Xsl*1i+1/((1/(Xm*1i))+(1/Rc)+(1/(Rr_Sn/Sn+Xrl_skew_Sn*1i)))));

    %%  TEMPERATURE (report)
    if p1_tek == 2, Alfa_conv=60; elseif p1_tek==4, Alfa_conv=50; elseif p1_tek==6, Alfa_conv=40; else, Alfa_conv=32; end
    Alfa_cond = 0.25/0.0003; A_ls = (2*hs+bs2)*L3*Ns; A_frame = pi*Dout*(L3+Tau3)*3;
    Teta_Co = 20 + P_Co/Alfa_cond/A_ls + Losses/Alfa_conv/A_frame;

    %%  ---- OBJECTIVE ----
    if ~isfinite(eff_h) || ~isreal(eff_h) || eff_h <= 0 || eff_h > 100
        error('imDesign:core:badEff','Invalid efficiency: %g', eff_h);
    end
    verim = eff_h;

    %%  ---- CONSTRAINT VIOLATIONS (10) ----
    g1 = ihl_ust(Bcs, 1.7);                 % Bcs < 1.7
    g2 = ihl_alt(Dshaft_in, 0);             % Dshaft_in > 0
    g3 = ihl_aralik(Sn,  0.005, 0.05);
    g4 = ihl_aralik(tbk, 1.75, 3.0);
    g5 = ihl_aralik(tLR, 0.70, 2.75);
    g6 = ihl_aralik(iLR, 6.0, 7.0);
    % g7: Fmts < Fmtr1 <= 1.1*Fmts  (band)
    base = max(1, abs(Fmts));
    g7 = ( max(0, Fmts - Fmtr1) + max(0, Fmtr1 - 1.1*Fmts) ) * 100 / base;
    % g8: Fmtr1 ~= Fmtr2 (equality, +-4% band)
    tol_F12 = 0.04*max(1,abs(Fmtr2));   % g8 band: +-4%
    g8 = max(0, abs(Fmtr1 - Fmtr2) - tol_F12) * 100 / max(1,abs(Fmtr2));
    % g9: Kst2 ~= Kst (equality, +-0.02 band)
    g9 = max(0, abs(Kst2 - Kst) - 0.02) * 100 / max(1e-3,abs(Kst));
    % g10: cos_phi_1n ~= cos_phi (equality, +-0.03 band)
    g10 = max(0, abs(cos_phi_1n - cos_phi) - 0.03) * 100 / max(1e-3,abs(cos_phi));

    gviol = [g1 g2 g3 g4 g5 g6 g7 g8 g9 g10];
    if any(~isfinite(gviol)) || any(~isreal(gviol))
        error('imDesign:core:badConstraint','Invalid constraint value.');
    end
    gviol = max(gviol, 0);

    out.eff_h=eff_h; out.Losses=Losses; out.P_Co=P_Co; out.P_Al=P_Al; out.P_iron=P_iron;
    out.Bcs=Bcs; out.Bg2=Bg2; out.Sn=Sn; out.tbk=tbk; out.tLR=tLR; out.iLR=iLR;
    out.Kst_in=Kst; out.Kst2=Kst2; out.Fmtr1=Fmtr1; out.Fmtr2=Fmtr2; out.Fmts=Fmts;
    out.cos_phi_1n=cos_phi_1n; out.cos_phi=cos_phi;
    out.Dshaft_in=Dshaft_in; out.Teta_Co=Teta_Co; out.W11=W11; out.ns=ns; out.L3=L3; out.p1_tek=p1_tek;

catch ME
    verim = 0; gviol = 1e6 * ones(1, NC); out.err = ME.message;
end
end

% --- Constraint-violation helpers (penalty measure, 0 = satisfied) ---
function v = ihl_ust(deger, ust)    % upper-bound violation
    v = max(0, deger - ust) * 100 / max(1e-9, abs(ust));
end
function v = ihl_alt(deger, alt)    % lower-bound violation
    if alt == 0, v = max(0, alt - deger) * 1000;
    else, v = max(0, alt - deger) * 100 / max(1e-9, abs(alt)); end
end
function v = ihl_aralik(deger, alt, ust)   % range (lower + upper) violation
    v = ihl_alt(deger, alt) + ihl_ust(deger, ust);
end
