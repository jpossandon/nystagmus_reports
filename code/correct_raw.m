function [xgaz,ygaz] = correct_raw(xraw,yraw,caldata)

xgazaux = caldata.ux'*[ones(1,size(xraw'-caldata.rawCenter(1),2));xraw'-caldata.rawCenter(1);yraw'-caldata.rawCenter(2);(xraw'-caldata.rawCenter(1)).^2;(yraw'-caldata.rawCenter(2)).^2];
ygazaux = caldata.uy'*[ones(1,size(yraw'-caldata.rawCenter(2),2));xraw'-caldata.rawCenter(1);yraw'-caldata.rawCenter(2);(xraw'-caldata.rawCenter(1)).^2;(yraw'-caldata.rawCenter(2)).^2];

if strcmp(caldata.calibType,'HV9')
    xgaz    = nan(1,length(xgazaux));
    ygaz    = nan(1,length(ygazaux));
    for ii = 1:4
        switch ii   %cuadrants
            case 1
                auxindx = find(xgazaux<0 & ygazaux<0);
            case 2
                auxindx = find(xgazaux>0 & ygazaux<0);
            case 3
                auxindx = find(xgazaux<0 & ygazaux>0);
            case 4
                auxindx = find(xgazaux>0 & ygazaux>0);
        end
        xgaz(auxindx) = xgazaux(auxindx)+caldata.m(ii).*xgazaux(auxindx).*ygazaux(auxindx);
        ygaz(auxindx) = ygazaux(auxindx)+caldata.n(ii).*xgazaux(auxindx).*ygazaux(auxindx);
    end
else
    xgaz = xgazaux;
    ygaz = ygazaux;
end
xgaz    = xgaz+caldata.rect(3)/2;
ygaz    = (ygaz+caldata.rect(4)/2);

%%
% remove extreme values around NaNs
xnan = find(isnan(xgaz));
xnan([xnan==1 | xnan==2 | xnan==3 | xnan==4 | xnan==5 | xnan==6 | xnan==7 | xnan==8| xnan==9 | xnan==10 | xnan==length(xgaz) | xnan==length(xgaz)-1 | xnan==length(xgaz)-2  | xnan==length(xgaz)-3 | xnan==length(xgaz)-4 | xnan==length(xgaz)-5 | xnan==length(xgaz)-6 | xnan==length(xgaz)-7 | xnan==length(xgaz)-8 | xnan==length(xgaz)-9 | xnan==length(xgaz)-10]) = [];
remx = union([xnan-1,xnan-2,xnan-3,xnan-4,xnan-5,xnan-6,xnan-7,xnan-8,xnan-9,xnan-10],[xnan+1,xnan+2,xnan+3,xnan+4,xnan+5,xnan+6,xnan+7,xnan+8,xnan+9,xnan+10]);
xgaz(remx) = NaN;

ynan = find(isnan(ygaz));
ynan([ynan==1 | ynan==2 | ynan==3 | ynan==4 | ynan==5 | ynan==6 | ynan==7 | ynan==8| ynan==9 | ynan==10 | ynan==length(ygaz) | ynan==length(ygaz)-1 | ynan==length(ygaz)-2  | ynan==length(ygaz)-3 | ynan==length(ygaz)-4 | ynan==length(ygaz)-5 | ynan==length(ygaz)-6 | ynan==length(ygaz)-7 | ynan==length(ygaz)-8 | ynan==length(ygaz)-9 | ynan==length(ygaz)-10]) = [];
remy = union([ynan-1,ynan-2,ynan-3,ynan-4,ynan-5,ynan-6,ynan-7,ynan-8,ynan-9,ynan-10],[ynan+1,ynan+2,ynan+3,ynan+4,ynan+5,ynan+6,ynan+7,ynan+8,ynan+9,ynan+10]);
ygaz(remy) = NaN;
    