function [pMeanAllTemp,pStdAllTemp] = MSDScalingCoeffTemp(VelFieldData,Sampling,dt)

% Calculate scaling coefficient of the MSD for different time windows:

Times = [1:Sampling:size(VelFieldData(1).tau,1)];
tau = dt*[1:size(VelFieldData(1).MSDTempAll,1)];

Offset = 2;
pMeanAllTemp = NaN(size(VelFieldData(1).MSDTempAll,1),length(VelFieldData),length(Times));
pStdAllTemp  = NaN(size(VelFieldData(1).MSDTempAll,1),length(VelFieldData),length(Times));

for m = 1:length(Times)
    for i = 1:length(VelFieldData)
        pGroup = NaN(size(VelFieldData(i).MSDTempAll,1) ,size(VelFieldData(i).MSDTempAll,3));
        for j = 1:size(VelFieldData(i).MSDTempAll,3)
            p = [];
            for k = Offset:size(VelFieldData(i).MSDTempAll,1)
                if length(VelFieldData(i).MSDTempAll(1:k,Times(m),j)) - sum(isnan(VelFieldData(i).MSDTempAll(1:k,Times(m),j))) > Offset
                    %pTemp = polyfit(log10(VelFieldData(i).tau(1:k)),log10(VelFieldData(i).MSDTempAll(1:k,Times(m),j)),1);
                    pTemp = polyfit(log10(tau(1:k)),log10(VelFieldData(i).MSDTempAll(1:k,Times(m),j)),1);
                    p(k) = pTemp(1);
                end
            end
            p(p==0) = NaN;
            pGroup(1:length(p),j) = p;
        end
        pMeanAllTemp(1:size(pGroup,1),i,m) = nanmean(pGroup,2);
        pStdAllTemp(1:size(pGroup,1),i,m) = nanstd(pGroup,[],2);
    end
end
