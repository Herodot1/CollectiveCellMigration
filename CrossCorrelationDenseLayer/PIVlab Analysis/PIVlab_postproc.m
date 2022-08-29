% Vector map post processing in PIVlab
function [u_out,v_out] = PIVlab_postproc (u,v,calu,calv, valid_vel, do_stdev_check,stdthresh, do_local_median,neigh_thresh)
%% velocity limits
if numel(valid_vel)>0 %velocity limits were activated
    umin=valid_vel(1);
    umax=valid_vel(2);
    vmin=valid_vel(3);
    vmax=valid_vel(4);
    u(u*calu<umin)=NaN;
    u(u*calu>umax)=NaN;
    v(u*calu<umin)=NaN;
    v(u*calu>umax)=NaN;
    v(v*calv<vmin)=NaN;
    v(v*calv>vmax)=NaN;
    u(v*calv<vmin)=NaN;
    u(v*calv>vmax)=NaN;
end

%% stddev check
if do_stdev_check==1
    meanu=nanmean(u(:));
    meanv=nanmean(v(:));
    std2u=nanstd(reshape(u,size(u,1)*size(u,2),1));
    std2v=nanstd(reshape(v,size(v,1)*size(v,2),1));
    minvalu=meanu-stdthresh*std2u;
    maxvalu=meanu+stdthresh*std2u;
    minvalv=meanv-stdthresh*std2v;
    maxvalv=meanv+stdthresh*std2v;
    u(u<minvalu)=NaN;
    u(u>maxvalu)=NaN;
    v(v<minvalv)=NaN;
    v(v>maxvalv)=NaN;
end

%% local median check
if do_local_median==1    
    % original
    Kernelsize = 3;
    neigh_filt=medfilt2(u,[Kernelsize,Kernelsize],'symmetric');
    neigh_filt=inpaint_nans(neigh_filt);
    neigh_filt=abs(neigh_filt-u);
    u(neigh_filt>neigh_thresh)=nan;
    
    neigh_filt=medfilt2(v,[Kernelsize,Kernelsize],'symmetric');
    neigh_filt=inpaint_nans(neigh_filt);
    neigh_filt=abs(neigh_filt-v);
    v(neigh_filt>neigh_thresh)=nan;
    
    
%     % matPIV version
%     m = 3; % kernelsize
%     nu=zeros(size(u)+2*floor(m/2))*nan;
%     nv=zeros(size(u)+2*floor(m/2))*nan;
%     nu(floor(m/2)+1:end-floor(m/2),floor(m/2)+1:end-floor(m/2))=u;
%     nv(floor(m/2)+1:end-floor(m/2),floor(m/2)+1:end-floor(m/2))=v;
%     
%     INx=zeros(size(nu));
%     INx(floor(m/2)+1:end-floor(m/2),floor(m/2)+1:end-floor(m/2))=zeros(size(u));
%     
%     prev=isnan(nu); previndx=find(prev==1);
%     U2=nu+i*nv; teller=1; [ma,na]=size(U2); histo=zeros(size(nu));
%     histostd=zeros(size(nu));hista=zeros(size(nu));histastd=zeros(size(nu));
%     %fprintf([' Local ',stat,' filter running: '])
%     ff = 1; % use median as reference
%     for ii=m-1:1:na-m+2
%         for jj=m-1:1:ma-m+2
%             if INx(jj,ii)~=1
%                 
%                 tmp=U2(jj-floor(m/2):jj+floor(m/2),ii-floor(m/2):ii+floor(m/2));
%                 tmp(ceil(m/2),ceil(m/2))=NaN;
%                 if ff==1
%                     usum=mnanmedian(tmp(:));
%                 elseif ff==2
%                     usum=mnanmean(tmp(:));
%                 end
%                 histostd(jj,ii)=mnanstd(tmp(:));
%             else
%                 usum=nan; tmp=NaN; histostd(jj,ii)=nan;
%             end
%             %         u1=real(usum).^2 - real(U2(jj,ii)).^2;
%             %         v1=imag(usum).^2 - imag(U2(jj,ii)).^2;
%             %
%             %         histo(jj,ii)=u1+i*v1;
%             histo(jj,ii)=usum;
%             %histostd(jj,ii)=mnanstd(real(tmp(:))) + i*mnanstd(imag(tmp(:)));
%             
%             %th1=angle(usum); th2=angle(U2(jj,ii));
%             %if th1<0, th1=2*pi+th1; end
%             %if th2<0, th2=2*pi+th2; end
%             %hista(jj,ii)=(th1-th2);
%             %if hista(jj,ii)<0, hista(jj,ii)=2*pi+hista(jj,ii); end
%             %histastd(jj,ii)=mnanstd(abs(angle(tmp(:))));
%         end
%         %fprintf('.')
%         
%     end
%     
%     %%%%%%%% Locate gridpoints with a higher value than the threshold
%     
%     %[cy,cx]=find((real(histo)>threshold*real(histostd) | ...
%     %    imag(histo)>threshold*imag(histostd)));
%     [cy,cx]=find( ( real(U2)>real(histo)+neigh_thresh*real(histostd) |...
%         imag(U2)>imag(histo)+neigh_thresh*imag(histostd) |...
%         real(U2)<real(histo)-neigh_thresh*real(histostd) |...
%         imag(U2)<imag(histo)-neigh_thresh*imag(histostd) ) );
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     for jj=1:length(cy)
%         %uv2(jj)=u(cy(jj),cx(jj)); vv2(jj)=v(cy(jj),cx(jj));
%         %xv2(jj)=x(cy(jj),cx(jj)); yv2(jj)=y(cy(jj),cx(jj));
%         % Now we asign NotANumber (NaN) to all the points in the matrix that
%         % exceeds our threshold.
%         nu(cy(jj),cx(jj))=NaN;  nv(cy(jj),cx(jj))=NaN;
%     end
%     u=nu(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
%     v=nv(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
%     
%     rest=length(cy);
%     
%     rest2=sum(isnan(u(:)))-sum(prev(:));


end

%% Gradient filter
%{
if do_gradient==1
    u_filled=inpaint_nans(u);
    v_filled=inpaint_nans(v);
    gradient_filt_x =abs(gradient(u_filled));
    gradient_filt_y =abs(gradient(v_filled));
    u(gradient_filt_x>neigh_thresh)=nan;
    v(gradient_filt_y>neigh_thresh)=nan;
end
%}

u(isnan(v))=NaN;
v(isnan(u))=NaN;
u_out=u;
v_out=v;