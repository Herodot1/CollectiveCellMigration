function [pMeanAll,pStdAll] = MSDScalingCoeff(VelFieldData)

Offset = 2;
pMeanAll = NaN(1700,length(VelFieldData));
pStdAll  = NaN(1700,length(VelFieldData));
tic
for i = 1:length(VelFieldData)
    pGroup = NaN(size(VelFieldData(i).tau,1) ,size(VelFieldData(i).MSDAll,2));
    for j = 1:size(VelFieldData(i).MSDAll,2)
        p = [];
        for k = Offset:size(VelFieldData(i).tau,1) 
            if length(VelFieldData(i).MSDAll(1:k,j)) - sum(isnan(VelFieldData(i).MSDAll(1:k,j))) > Offset
                pTemp = polyfit(log10(VelFieldData(i).tau(1:k)),log10(VelFieldData(i).MSDAll(1:k,j)),1);     
                p(k) = pTemp(1);
            end
        end
        p(p==0) = NaN;
        pGroup(1:length(p),j) = p;
    end
    pMeanAll(1:size(pGroup,1),i) = nanmean(pGroup,2);
    pStdAll(1:size(pGroup,1),i) = nanstd(pGroup,[],2);
end
toc