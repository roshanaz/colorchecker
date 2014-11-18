function [XYZim]=RGBtoXYZconver(RGBim)
var_R =  RGBim(:,1) ; 
var_G =  RGBim(:,2) ;
var_B =  RGBim(:,3) ;

if ( var_R > 0.04045 ) var_R = ( ( var_R + 0.055 ) / 1.055 ) ^ 2.4;
else                   var_R = var_R / 12.92;
end
if ( var_G > 0.04045 ) var_G = ( ( var_G + 0.055 ) / 1.055 ) ^ 2.4;
else                   var_G = var_G / 12.92;
end

if ( var_B > 0.04045 ) var_B = ( ( var_B + 0.055 ) / 1.055 ) ^ 2.4;
else                   var_B = var_B / 12.92;
end

var_R = var_R * 100;
var_G = var_G * 100;
var_B = var_B * 100;


X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;

XYZim = [X;Y;Z];
end