function [dmin, dmax, isim] = im_bounds()
%IM_BOUNDS  Lower/upper bounds and names of the 16 continuous design variables.
%   Lamda_katsayi and KD_katsayi are now variables (previously fixed) -> D = 16.
%   Lengths are in mm; im_core converts them to meters via /1000.
%   KD_katsayi range given as 0.62/0.61 .. 0.62/0.63; numerically
%   [0.62/0.63, 0.62/0.61] = [0.98413, 1.01639] (ordered so lower < upper).
%   rotor_slot: 0 < deg < 360; since 0 and 360 are singular, bounded to [1,359].
%        1      2          3     4    5      6     7    8     9    10   11   12   13  14    15   16
%       Lamk   KDk        Bgk   Kst  Jcosk  bos   hos  hw    Bts  Jb   Btr  hor  bor deg   Bcr  k
    dmin=[0.80  0.62/0.63  0.85  0.10 0.75   1.75  0.5  0.75  1.20 1.5  1.20 0.3  1.0 1     1.20 1.0];
    dmax=[1.20  0.62/0.61  1.20  1.00 1.25   3.25  1.0  3.75  1.70 4.0  1.70 0.8  2.0 359   1.50 1.2];
    isim = {'Lamda_katsayi','KD_katsayi','Bg_katsayi','Kst','Jcos_katsayi', ...
            'bos[mm]','hos[mm]','hw[mm]','Bts[T]','Jb[A/mm2]','Btr[T]','hor[mm]', ...
            'bor[mm]','rotor_slot[deg]','Bcr[T]','k'};
end
