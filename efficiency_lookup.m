function [eff] = efficiency_lookup(eff_sinifi,cift_kutup_sayisi,guc)
%EFF_HESAPLAMA  Target efficiency (%) for an IE class from a cubic fit.
%   Returns the IEC IE1/IE2/IE3 target efficiency as a 3rd-order polynomial
%   in log10(power[kW]). Column index = pole-pair count (cift_kutup_sayisi);
%   rows are the cubic coefficients (highest order first).

IE1_katsayilari = [0.5234 0.5234 0.0786;
                   -5.0499 -5.0499 -3.5838;
                   17.4180	17.4180	17.2918;
                   74.3171 74.3171 72.2383];
IE2_katsayilari = [0.2972 0.0278 0.0148;
                   -3.3454 -1.9247 -2.4978;
                   13.0651 10.4395 13.2470;
                   79.0770 80.9761 77.5603];
IE3_katsayilari = [0.3569 0.0773 0.1252;
                   -3.3076 -1.8951 -2.6130;
                   11.6108 9.2984 11.9963;
                   82.2503 83.7025 80.4769];

if eff_sinifi == "IE1"
    eff = IE1_katsayilari(1,cift_kutup_sayisi)*(log10(guc/1000))^3+...
        IE1_katsayilari(2,cift_kutup_sayisi)*(log10(guc/1000))^2+...
        IE1_katsayilari(3,cift_kutup_sayisi)*(log10(guc/1000))^1+...
        IE1_katsayilari(4,cift_kutup_sayisi);
elseif eff_sinifi == "IE2"
    eff = IE2_katsayilari(1,cift_kutup_sayisi)*(log10(guc/1000))^3+...
        IE2_katsayilari(2,cift_kutup_sayisi)*(log10(guc/1000))^2+...
        IE2_katsayilari(3,cift_kutup_sayisi)*(log10(guc/1000))^1+...
        IE2_katsayilari(4,cift_kutup_sayisi);
elseif eff_sinifi == "IE3"
    eff = IE3_katsayilari(1,cift_kutup_sayisi)*(log10(guc/1000))^3+...
        IE3_katsayilari(2,cift_kutup_sayisi)*(log10(guc/1000))^2+...
        IE3_katsayilari(3,cift_kutup_sayisi)*(log10(guc/1000))^1+...
        IE3_katsayilari(4,cift_kutup_sayisi);
end
