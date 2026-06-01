%In this version of the slope-finding function, fit a line to the data
%using least squares and recover the slope (and intercept). Note that this
%is a renamed version of terminal_slope_ols.

function [slopels,b_ls,Vmid_index_ls,timevec_ls,Vmax_index]=terminal_slope_ls(temp_T,temp_VL)

%Get index for max VL value 
[maxv,Vmax_index] = max(temp_VL);
% Vmax_index = find(temp_VL == maxv);

T = temp_T;
viral_load = temp_VL;

zero_cross = viral_load(2:end).*viral_load(1:end-1) < 0;
zero_cross_ind = find(zero_cross);

if length(zero_cross_ind) >1
    T = T(zero_cross_ind(1): zero_cross_ind(end));
    viral_load = viral_load(zero_cross_ind(1): zero_cross_ind(end));
end

%Cut off extra data from >10 days after peak 
if (T(end) - T(Vmax_index)) > 10
    viral_load(Vmax_index + 11:end) = [];
    T(Vmax_index + 11:end) = [];
end

%Find index for cutoff value for the data input for LS, determined by when VL values start increasing
%more than 70% of the way into the time vector
%Note: Depending on the data, another option to add in is to cut off the
%data 60% of the way in if there are multiple increases in a row 
for i= 2:(length(viral_load))
     if viral_load(i-1)<viral_load(i) && (i-1)/length(viral_load)>=0.7% && temp_VL(i + 1)<temp_VL(i+2)&& temp_VL(i + 2)<temp_VL(i+3)
         midv = viral_load(i-1);
         break
     else
         midv = viral_load(end);
     end
end
Vmid_index_ls = min(find(viral_load == midv));

%Update T, temp_VL to only include max point to cutoff point for fitting 
T_tofit=T(Vmax_index:Vmid_index_ls);
temp_VL_tofit = temp_VL(Vmax_index:Vmid_index_ls);

%Make T, temp_VL, and the truncated versions into column vectors if they are not already 
if iscolumn(T) == 0
    T=T';
end
if iscolumn(temp_VL_tofit) == 0
    temp_VL_tofit = temp_VL_tofit';
end
if iscolumn(T_tofit) == 0
    T_tofit=T_tofit';
end
if iscolumn(temp_VL) == 0
    temp_VL = temp_VL';
end


%Fit a line to the truncated data with the condition that the line go through the max value 
%https://www.mathworks.com/matlabcentral/answers/94272-how-do-i-constrain-a-fitted-curve-through-specific-points-like-the-origin-in-matlab
%^Source for fitting with this specific condition 
n=1;
V(:,n+1) = ones(length(temp_VL_tofit),1,class(temp_VL_tofit));
for j=n:-1:1
    V(:,j) = T_tofit.*V(:,j+1);
end
C=V;
d=temp_VL_tofit;
A=[];
b=[];

Aeq=T(Vmax_index).^(n:-1:0); %Condition to go through max value 
beq=temp_VL(Vmax_index);

p=lsqlin(C,d,A,b,Aeq,beq);
yhat = polyval(p,T);
slopels = p(1);
b_ls = p(2);
timevec_ls = T_tofit;

end