function [Tbreak, Yield, Purity, Rs]  = CalArea2(A,Bot,Pu)

%Bot: a peak baseline (lowcut) for fast calculation of peak width
%Recommend Bot value as 0.5% of peak height
%Pu: Purity cut

pur1_p1 = 0;
pur2_p2 = 0;
Rs=0;

T = A(:,1);
y1 = A(:,2);
y2 = A(:,3);

y1area = trapz(T,y1);
y2area = trapz(T,y2);

%finds the time when the peak is at max height
tPeak1 = mean(T(find(y1 == max(y1))));
tPeak2 = mean(T(find(y2 == max(y2))));

%finds the time when the signal is beyond the start cutoff (bot*peakheight)
indBase1 = find(y1 >= Bot*max(y1));
cut1 = min(indBase1); %cut1
cut2 = max(indBase1); %cut2

indBase2 = find(y2 >= Bot*max(y2));
cut3 = min(indBase2);
cut4 = max(indBase2);


%we need to include the case where the peaks overlap -> if cut3 < cut2
if cut3 < cut2
    %we want to make both cuts equal, and at the point where the signal
    %values are equal
%        cut = [cut1 cut2 cut3 cut4];

    intersect = round(mean(find(abs(y1(cut3:cut2)-y2(cut3:cut2)) < 0.1))); %loose tolerance
    %Not sure if this works yet!!!!! may need to use a version of isclose(), could do >= and <= with a small tol
    cut2 = cut3+intersect;
    cut3 = cut2;
    %adjust indBases 

%     cut_r = [cut1 cut2 cut3 cut4]'

end
% a = T(cut1)
% d = T(cut4)
% b  = T(cut2)
% c = T(cut3)

Wb1 = T(cut2)-T(cut1);
Wb2 = T(cut4)-T(cut3);
% if(Wb1 < 0 || Wb2 < 0)
%     fprintf('overlap is too great to find meaningful intersection point')
%     Tbreak=[T(cut1) T(cut2); T(cut3) T(cut4)];
%     Purity=0;
%     Yield=0;
%     return
% end
% Rs = (tPeak2-tPeak1)/(0.5*(Wb1+Wb2));

%And now we need to calculate the actual purity and yield values
while pur1_p1 < Pu || pur2_p2 < Pu
    col1_p1 = trapz(T(cut1:cut2),y1(cut1:cut2));
    col2_p1 = trapz(T(cut1:cut2),y2(cut1:cut2));
    col1_p2 = trapz(T(cut3:cut4),y1(cut3:cut4));
    col2_p2 = trapz(T(cut3:cut4),y2(cut3:cut4));

    pur1_p1 = col1_p1/(col1_p1+col2_p1);
    pur2_p2 = col2_p2/(col1_p2+col2_p2);
    %pur2_p1 = 1-pur1_p1;
    %pur1_p2 = 1-pur2_p2;
    if pur1_p1 < Pu
        cut2 = cut2-1;
    elseif pur2_p2 < Pu
        cut3 = cut3+1;
    end
end
Purity = [pur1_p1 pur2_p2];

yield1_p1 = col1_p1/y1area;
yield2_p2 = col2_p2/y2area;
Yield = [yield1_p1 yield2_p2];
%yield2_p1 = col2_p1/y2area;
%yield1_p2 = col1_p2/y1area;

Tbreak = [T(cut1)  T(cut2);
            T(cut3) T(cut4)];


end

% for comp=1:numComp
%                     Tbreak(1,comp) = T(cut(1,comp));
%                     Tbreak(2,comp) = T(cut(2,comp));
%                 end
%      