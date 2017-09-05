model IkedaEq
  // B1
  constant Real RTOT = 20, KR = 0.3, RTOP = 3, KL = 0.2, DEN = 1;
  Real QCO0, QCO, PAS, PVS, PAP, PVP;
  //B2
  constant Real XHB = 15, FO2I = 0.21, FCOI = 0, PBA = 760, PBL = PBA - 47, VAL = 3, MRCO = 0.2318, MRO2 = 0.2591, VI0 = 5;
  Real FCOA, FO2A, UCOV, UO2V, UCOA, UO2A, VI, XCO3;
  Real k1, k6, k3, k4, k2, k5, H, VR, UHBO, DCLA, UHB, PCOA, PO2A, PHA;
  // B3
  // constant parameters
  constant Real QIN = 0.001, QVIN = 0, QIWL = 0.0005, QMWP = 0.0005, QLF0 = 0.002, VIF0 = 8.8, CFC = 0.007, VRBC = 1.8, CRAV = 5.93;
  Real VB, HT, VEC, PC, XPP, YPLC, XPIF, YPLF, YPLV, YPLG, YPG, VP_plus, VP_minus, QCFR, Pf, x, PIF, QLF, QPLC, PPCO, PICO;
  Real VIN, VIF, ZPP, ZPIF, ZPLG, ZPG, VP;
  // B4
  constant Real CSM = 0.0003, YNIN = 0.12, CKEI = 0.001, YKIN = 0.047, XGLO = 108, YINS = 0, YGLI = 0, YMNI = 0, YURI = 0.15, CGL1 = 1, CGL2 = 1, CGL3 = 0.03, CHEI = 5, CBFI = 10 ^ (-9);
  Real ZNE, ZKE, ZKI, ZHI, YINT, ZGLE, ZMNE, ZURE, VIC;
  Real y4, F41, YGLU, XNE, XKI, XKE, YHI, PHI, YGLS, z4, XGLE, OSMP, XMNE, XURE, VTW, YMNU, YURU, QIC;
  // B5
  constant Real YCAI = 0.007, YCLI = 0.1328, YMGI = 0.008, YOGI = 0.01, YPOI = 0.025, YSOI = 0.02;
  Real F50, YCO3, YCA, YMG, YSO4, YPO4, YORG, XCAE, XMGE, XSO4, XPO4, XOGE, XCLE, XCLA, STBC, YCLU;
  Real ZCAE, ZMGE, ZSO4, ZPO4, ZOGE, ZCLE;
  // B6
  constant Real CPRX = 0.2, YNH0 = 0.024, YTA0 = 0.0068;
  Real PHU2, PHA1, PHU1;
  Real YNH4, F62, YTA1, STPO, PHU, STPG, YTA, YNOD, YNU, GFR1_6, YND, YKU, YKD, OSMU, QWD, QWU, YNH;
  // B7
  Real CKAL = 0.5, CNAL = 0.1, COAD = 0.5, CPAD = 1, CPAL = 0.01, CPVL = 0.1, GFR0 = 0.1, ACTH = 1, VEC0 = 11, TADH = 30, TALD = 30;
  Real ALD0, ADH0;
  Real ADH, ALD, THDF, GFR1, GFR, z7, ALD1;
initial equation
  // B2
  FCOA = 0.0561;
  //
  FO2A = 0.1473;
  //
  UCOV = 0.6075;
  // l(STPD)/l.blood
  UO2V = 0.1515;
  // l(STPD)/l.blood
  VI = 5;
  // l/min
  XCO3 = 24;
  // mEq/l
  // B3
  VIN = 0.01;
  // l
  VIF = VIF0;
  // l
  ZPP = 154;
  // g
  ZPIF = 176;
  // g
  ZPLG = 70;
  // g
  ZPG = 20;
  // g
  VP = 2.2;
  // l
  // B4
  ZNE = 1540;
  // mEq
  ZKE = 49.5;
  // mEq
  ZKI = 2800;
  // mEq
  ZHI = 100;
  YINT = 0;
  ZGLE = 66;
  // mg
  ZMNE = 0;
  // mM
  ZURE = 77.5;
  // mM
  VIC = 20;
  // l
  // B5
  ZCAE = 55;
  // mEq
  ZMGE = 33;
  // mEq
  ZSO4 = 11;
  // mEq
  ZPO4 = 12.1;
  // mM
  ZOGE = 66;
  // mM
  ZCLE = 1144;
  // mEq
  // B6
  PHU2 = 6;
  PHA1 = 7.4;
  PHU1 = 6;
  // B7
  ALD0 = 0;
  ADH0 = 0;
equation
  //*****************************************************************
  // BLOCK 1 - CARDIOVASCULAR SYSTEM
  // Input = VB
  // Output = QCO, PAS, PVS, PVP, PAP
  //*****************************************************************
  // constant parameters
  QCO0 = VB * DEN;
  // l/min
  QCO = QCO0 + 1;
  // l/min
  PAS = 20 + RTOT * QCO0;
  // mmHg
  PVS = max(0, (-10.33) + QCO0 / KR);
  // mmHg
  PAP = 8 + RTOP * QCO0;
  // mmHg
  PVP = max((-16) + QCO0 / KL, 0);
  // mmHg
  //*****************************************************************
  // BLOCK 2 - RESPIRATION
  // Input = STBC, QCO, VTW
  // Output = DCLA, PCOA, XCO3, PHA
  //*****************************************************************
  // equations (3) to (8)
  der(FCOA) = (VI * (FCOI - FCOA) + 863 / (PBA - 47) * QCO * (UCOV - UCOA)) / VAL;
  der(FO2A) = (VI * (FO2I - FO2A) + 863 / (PBA - 47) * QCO * (UO2V - UO2A)) / VAL;
  der(UCOV) = (MRCO + QCO * (UCOA - UCOV)) / VTW;
  der(UO2V) = ((-MRO2) + QCO * (UO2A - UO2V)) / VTW;
  UCOA = 6.732 * 10 ^ (-4) * PCOA + 0.02226 * XCO3;
  // l(STPD)/l.blood
  UO2A = 3.168 * 10 ^ (-5) * PO2A + UHBO;
  // l(STPD)/l.blood
  // equation (1)
  der(VI) = (VR * VI0 - VI) / 2;
  // equation (2) and F21
  k1 = if PHA <= 7.4 then 0.22 else 0.0258;
  k6 = if PHA <= 7.4 then -12.734 else -5.003;
  k3 = 0.58;
  k4 = 3.496;
  k2 = if PCOA > 40 then 1 else 0.0396;
  k5 = if PCOA > 40 then -32.08 else 160.11;
  H = 10 ^ (9 - PHA);
  VR = k1 * H + k2 * (k3 + k4 / (PO2A - 32)) * (PCOA + k5) + k6;
  assert(VR >= 0, "VR out of bounds! Original LIMIT VR >= 0; ");
  // F23
  /*
    f=(1-exp(-PO2A*g))^2
    g=0.0066815*PHA^3-0.10098*PHA^2+0.44921*PHA-0.454
    */
  UHBO = UHB * (1 - exp((-PO2A * 0.0066815 * PHA ^ 3) - 0.10098 * PHA ^ 2 + 0.44921 * PHA - 0.454)) ^ 2;
  // l.O2(STPD)/l.blood
  DCLA = XCO3 - STBC;
  // mEq/l
  UHB = XHB / 75;
  // l.O2(STPD)/l.blood
  PCOA = FCOA * (PBA - 47);
  // mmHg
  PO2A = FO2A * (PBA - 47);
  // mmHg
  // equation (9) and F22
  PHA = 6.1 + log10(XCO3 / (0.03 * PCOA));
  //
  // F24
  der(XCO3) = STBC - (0.527 * XHB + 3.7) * (PHA - 7.4) + 0.375 * (UHB - UHBO) / 0.02226 - XCO3;
  //*****************************************************************
  // BLOCK 3 - EXTRACELLULAR SPACE
  // Input = PAS, PVS, QIC, QWU
  // Output = PPCO, VB// VEC, XPP
  //*****************************************************************
  VB = VRBC + VP;
  HT = VRBC / VB;
  //
  VEC = VP + VIF;
  // l
  der(VIN) = QIN - VIN / 10;
  der(VIF) = QCFR - QLF - QIC;
  PC = (CRAV * PVS + PAS) / (1 + CRAV);
  // mmHg
  der(ZPP) = YPLF - YPLG - YPLV - YPLC;
  XPP = ZPP / VP;
  // g/l
  YPLC = QPLC * (XPP - XPIF);
  // g/min
  XPIF = ZPIF / VIF;
  // g/l
  der(ZPIF) = YPLC - YPG - YPLF;
  YPLF = XPIF * QLF;
  // g/min
  YPLV = XPP * 0.00047 - 0.0329;
  // g/min
  YPLG = 0.00023 * (XPP - ZPLG);
  // g/min
  der(ZPLG) = (XPP - ZPLG) / 24;
  YPG = 0.0057 * (XPIF - ZPG);
  //g/min
  der(ZPG) = (XPIF - ZPG) / 150;
  // equation (10)
  VP_plus = VIN / 10 + QVIN + QMWP + QLF;
  // equation (11)
  VP_minus = QIWL + QWU + QCFR;
  // equation (12)
  der(VP) = VP_plus - VP_minus;
  // equation (13)
  QCFR = CFC * Pf;
  // l/min
  // equation (14)
  Pf = PC - PPCO - PIF + PICO;
  // mmHg
  // F31;  // NB: no line breaks in the equation for the actual program
  x = VIF / VIF0;
  //
  PIF = if x <= 0.9 then -15 elseif x > 0.9 and x <= 1 then 87 * x - 93.3
   elseif x > 1 and x <= 2 then -6.3 * (2 - x) ^ 10 else x - 2;
  // mmHg
  // F32
  QLF = QLF0 * (24 / (1 + exp(-0.4977 * PIF)));
  // l/min
  // F33
  QPLC = 2.768 * 10 ^ (-6) * PC ^ 2;
  // l/min
  // F34
  PPCO = 0.4 * XPP;
  // mmHg
  // F35
  PICO = 0.25 * XPIF;
  // mmHg
  //*****************************************************************
  // BLOCK 4 - INTRACELLULAR SPACE and ELECTROLYTES
  // Input = GFR, PHA, VEC, YKU, YNU
  // Output = QSMP, QIC, VTW, XKE, XNE, YGLU, YHI, YMNU, YNU
  //*****************************************************************
  // equation (15)
  der(ZNE) = YNIN - YNU + YHI;
  // equation (16)
  der(ZKE) = YKIN - YKU - y4;
  der(ZKI) = y4;
  y4 = z4 + CKEI * (2800 * F41 - ZKI);
  // F41
  F41 = 1 + 0.5 * log10(XKE / (56.744 - 7.06 * PHA));
  // F42
  YGLU = if XGLE * GFR < 0.65 then 0 else XGLE * GFR - 0.65;
  // mg/min
  XNE = ZNE / VEC;
  // mEq/l
  XKI = ZKI / VIC;
  // mEq/l
  XKE = ZKE / VEC;
  // mEq/l
  der(ZHI) = YHI;
  YHI = CHEI * (0.4 - PHA + PHI);
  PHI = -log10(CBFI * ZHI);
  der(YINT) = 1 / 1.50 * (XGLE - XGLO / 18 - YINT);
  YGLS = CGL1 * YINT + CGL2 * YINS;
  z4 = CGL3 * YGLS;
  der(ZGLE) = YGLI / 180 - YGLS - YGLU;
  XGLE = ZGLE / VEC;
  // mg/l
  OSMP = (XNE + XKE) * 1.86 + XGLE + XURE + XMNE + 9.73;
  // mOsm/l
  der(ZMNE) = YMNI - YMNU;
  der(ZURE) = YURI - YURU;
  XMNE = ZMNE / VEC;
  // mM/l
  XURE = ZURE / VTW;
  // mM/l
  VTW = VEC + VIC;
  // l
  YMNU = 1 * GFR * XMNE;
  // mM/min
  YURU = XURE * GFR * 0.6;
  // mM/min
  der(VIC) = QIC;
  QIC = CSM * ((-XNE) - XKE - XGLE + 10.5 + XKI);
  // l/min
  //*****************************************************************
  // BLOCK 5 - KIDNEY 2
  // Input = DCLA,GFR,PCOA,VEC,XCO3,XKE,XNE,XPP,YHI,YKU,YNH4,YNU, STPG
  // Output = STBC, YCO3, YORG, YPO4
  //*****************************************************************
  // F50
  F50 = (-PCOA / 120) + 4 / 3;
  // F51;  // NB: no line breaks in the equation for the actual program
  YCO3 = if XCO3 * GFR * F50 <= 2 then 0 elseif XCO3 * GFR * F50 > 2 and XCO3 * GFR * F50 <= 4 then 0.1638 * (XCO3 * GFR * F50 - 2) ^ 2.61 else XCO3 * GFR * F50 - 3;
  // mEq/min
  // F52
  YCA = if XCAE * GFR < 0.493 then 0 else XCAE * GFR - 0.493;
  // mEq/min
  // F53
  YMG = if XMGE * GFR < 0.292 then 0 else XMGE * GFR - 0.292;
  // mEq/min
  // F54
  YSO4 = if XSO4 * GFR < 0.08 then 0 else XSO4 * GFR - 0.08;
  // mEq/min
  // F55
  YPO4 = if XPO4 * GFR <= 0.11 then 5 / 22 * XPO4 * GFR else XPO4 * GFR - 0.085;
  //mM/min
  // F56
  YORG = if XOGE * GFR <= 0.6 then XOGE * GFR / 60 else XOGE * GFR / 3 - 0.19;
  // mM/min
  der(ZCAE) = YCAI - YCA;
  der(ZMGE) = YMGI - YMG;
  der(ZSO4) = YSOI - YSO4;
  der(ZPO4) = YPOI - YPO4;
  der(ZOGE) = YOGI - YORG;
  der(ZCLE) = YCLI - YCLU;
  XCAE = ZCAE / VEC;
  // mEq/l
  XMGE = ZMGE / VEC;
  // mEq/l
  XSO4 = ZSO4 / VEC;
  // mEq/l
  XPO4 = ZPO4 / VEC;
  // mM/l
  XOGE = ZOGE / VEC;
  // mM/l
  XCLE = ZCLE / VEC;
  // mEq/l
  XCLA = XCLE - DCLA;
  // mEq/l
  STBC = XCAE + XMGE - XSO4 - 1.8 * XPO4 - XOGE - XCLE + XNE + XKE - 0.2214 * XPP;
  // mEq/l
  YCLU = max(0, YNU + YKU - STPG + YNH4 - YCO3 + YCA + YMG - YSO4);
  // mEq/min
  //*****************************************************************
  // BLOCK 6 - KIDNEY 1
  // Input = ADH, ALD, GFR, OSMP, PHA, THDF, XKE, XNE, YCO3, YGLU, YORG,
  // YPO4, YMNU, YURU
  // Output = STPG, QWU, YKU, YNH, YNH4,YNU
  //*****************************************************************
  YNH4 = YNH0 * ((-0.5 * PHU1) + 4);
  // mEq/min
  // F62
  F62 = YTA0 * ((-2.5 * PHA1) + 19.5);
  // F63
  YTA1 = if PHU2 <= 4 then 0 elseif PHU2 > 4 and PHU2 <= 5 then F62 * (PHU2 - 4) else F62;
  // mEq/min
  der(PHU2) = PHU - PHU2;
  // F64
  STPO = YPO4 * (1 + 1 / (1 + 10 ^ (6.8 - PHA)));
  // mM/min
  // F65;  // NB: no line breaks in the equation for the actual program
  PHU = -log10(((-((10 ^ (-4.3) + 10 ^ (-6.8)) * (STPG - YPO4 - 1 / (1 + 10 ^ (PHA - 4.3)) * YORG) - 10 ^ (-6.8) * YPO4 - 10 ^ (-4.3) * YORG)) + (((10 ^ (-4.3) + 10 ^ (-6.8)) * (STPG - YPO4 - 1 / (1 + 10 ^ (PHA - 4.3)) * YORG) - 10 ^ (-6.8) * YPO4 - 10 ^ (-4.3) * YORG) ^ 2 - 4 * (STPG - YPO4 - 1 / (1 + 10 ^ (PHA - 4.3)) * YORG) * (10 ^ (-4.3) * 10 ^ (-6.8) * (STPG - YPO4 - 1 / (1 + 10 ^ (PHA - 4.3)) * YORG - YPO4 - YORG))) ^ 0.5) / 2 / (STPG - YPO4 - 1 / (1 + 10 ^ (PHA - 4.3)) * YORG));
  STPG = max(0, STPO + YORG - YTA);
  YTA = YTA1 + 0.009 + ALD * 0.001;
  der(PHA1) = (PHA - PHA1) / 200;
  der(PHU1) = (PHU - PHU1) / 300;
  YNOD = max(0, YTA1 + YNH4 - YCO3);
  // mEq/min
  YNU = max(YND * 0.116 - YNOD, 0);
  // mEq/min
  GFR1_6 = THDF * GFR * CPRX;
  // l/min
  YND = XNE * GFR1_6 * 0.5 * 0.9 - ALD * 0.09;
  // mEq/min
  YKU = 0.39 * YKD;
  // mEq/min
  YKD = ALD * 0.018 * XKE + 0.9 * 0.5 * GFR1_6 * XKE;
  // mEq/min
  OSMU = (YGLU + YURU + YMNU + 1.86 * (YKU + YNU)) / QWU;
  // mOsm/l
  QWD = (YGLU + YURU + YMNU + (YND + YKD) * 1.86 + 0.32) / OSMP;
  // l/min
  QWU = QWD - QWD * 0.9 * ADH;
  // l/min
  YNH = 0.5 * GFR1_6 * XNE;
  // mEq/min
  //*****************************************************************
  // BLOCK 7 - CONTROLLER of RENAL FUNCTION
  // Input = OSMP, PAS, PPCO, PVP, VEC, XKE, XNE,YNH
  // Output = ADH, ALD, GFR, THDF
  //*****************************************************************
  // F71
  ADH = 1.1 / (1 + exp(-0.5 * (ADH0 + 4.605)));
  // F72
  ALD = 10 / (1 + exp(-0.4394 * (ALD0 - 5)));
  // F73
  THDF = if PPCO <= 28 then (-5 * (PPCO / 28 - 1)) + 1 else 1;
  // l
  // F74;  // NB: no line breaks in the equation for the actual program
  GFR1 = if PAS < 40 then 0 elseif PAS >= 40 and PAS < 80 then 0.02 * PAS - 0.8
   elseif PAS >= 80 and PAS < 100 then (-0.0005 * (PAS - 100) ^ 2) + 1 else 1;
  GFR = GFR1 * GFR0 * VEC / VEC0;
  // l/min
  der(ALD0) = (ALD1 - ALD0) / TALD;
  der(ADH0) = (z7 - ADH0) / TADH;
  z7 = (OSMP - 287) * COAD - (PVP - 4) * CPAD;
  ALD1 = (ACTH - 1) * 1 + (XKE - 4.5) * CKAL - (PVP - 4) * CPVL - (YNH - 1.4) * CNAL - (PAS - 100) * CPAL;
  annotation(uses(Modelica(version = "3.2.1")), experiment(StartTime = 0, StopTime = 180, Tolerance = 1e-06, Interval = 0.36));
end IkedaEq;