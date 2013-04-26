function normScr = objfun(objType,penType,lab,allMasks,p,lambda)

warning off;

% time span of simultaneous streams
st=inf; ed=-inf;
for k=1:length(allMasks)
    if any(any(allMasks{k}.pitch))
        st = min(st, find(sum(allMasks{k}.pitch,1)>0,1,'first'));
        ed = max(ed, find(sum(allMasks{k}.pitch,1)>0,1,'last'));    
    end
end
nFrame = ed-st;

fprintf(' (Obj=%s,Pen=%s,p=%.2f(%d),lambda=%.2f)',objType,penType,p,nFrame,lambda);

N = size(lab,1);
scr = zeros(1,N);
pen = zeros(1,N);
for k=1:N   
    if mod(k,round(N/5))==0, fprintf('.'); end
    
    % hypothesis grouping
    dat1=[]; dat2=[]; pitch1=[]; pitch2=[];
    for i=1:length(lab(k,:))
        if lab(k,i)==1
            dat1 = [dat1 allMasks{i}.dat];
            pitch1 = [pitch1; allMasks{i}.pitch];
        else
            dat2 = [dat2 allMasks{i}.dat];
            pitch2 = [pitch2; allMasks{i}.pitch];
        end
    end    
                
    % penalized objective function
    if sum(lab(k,:))>0 && sum(lab(k,:))<length(allMasks) % no all-0 or all-1 label        
        % objective
        switch objType
            case 'trace'                
                n1=size(dat1,2); n2=size(dat2,2);
                grand_me = mean([dat1 dat2],2);                
                m1=mean(dat1,2); m2=mean(dat2,2);
                Sw = (n1-1)*cov(dat1')+(n2-1)*cov(dat2');
                Sb = n1*(m1-grand_me)*(m1-grand_me)' + n2*(m2-grand_me)*(m2-grand_me)';
                scr(k) = trace(Sb/(Sw+eps));                 
            case 'mahal'
                sigma = cov([dat1 dat2]');
                m1 = mean(dat1,2);
                m2 = mean(dat2,2);
                scr(k) = (m1-m2)'*1/sigma*(m1-m2);  
            case 'T2'
                N1 = size(dat1,2);
                N2 = size(dat2,2);
                sigma = cov([dat1 dat2]');
                m1 = mean(dat1,2);
                m2 = mean(dat2,2);
                scr(k) = N1*N2/(N1+N2)*(m1-m2)'*1/sigma*(m1-m2);
            case 'GLR'
                N1=size(dat1,2); N2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                sigma1=cov(dat1'); sigma2=cov(dat2');
                sigma=cov([dat1 dat2]');
                lamb_mu = ( 1 + N1*N2/(N1+N2)^2*(m1-m2)'*1/(sigma+eps)*(m1-m2) )^( (-(N1+N2))/2 );
                alpha = N1/(N1+N2);
                lamb_sigma = ( (det(sigma1)+eps)^alpha * (det(sigma2)+eps)^(1-alpha) / (det(sigma)+eps) )^((N1+N2)/2);
                scr(k) = -log(lamb_mu*lamb_sigma+eps);
            case 'KL'
                m1=mean(dat1,2); m2=mean(dat2,2);
                sigma1=cov(dat1'); sigma2=cov(dat2');
                sigma1=sigma1+eps; sigma2=sigma2+eps;
                scr(k) = .5*(m1-m2)'*(inv(sigma1)+inv(sigma2))*(m1-m2)+...
                    .5*trace(sigma2/sigma1+sigma1/sigma2-2*eye(size(sigma1,1)));
            case 'Bhatt'
                m1=mean(dat1,2); m2=mean(dat2,2);
                sigma1=cov(dat1'); sigma2=cov(dat2');
                sigma1=sigma1+eps; sigma2=sigma2+eps;
                scr(k) = .25*(m1-m2)'*1/(sigma1+sigma2)*(m1-m2)+...
                    .5*log( det(sigma1+sigma2)/( 2*sqrt(abs(det(sigma1*sigma2))+eps)+eps )+eps );
            case 'L2'
                N1=size(dat1,2); N2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                scr(k) = N1*N2/(N1+N2)*sum((m1-m2).^2);
            case 'L_inf'
                m1=mean(dat1,2); m2=mean(dat2,2);
                scr(k) = max(abs(m1-m2));
             case 'mse'
                scr(k) = sum(sum( (dat1-repmat(mean(dat1,2),1,size(dat1,2))).^2 ))+...
                    sum(sum( (dat2-repmat(mean(dat2,2),1,size(dat2,2))).^2 ));
                scr(k) = -scr(k);
            case 'CH'
                grand_me = mean([dat1 dat2],2);
                n1=size(dat1,2); n2=size(dat2,2);
                Sw = (n1-1)*cov(dat1')+(n2-1)*cov(dat2');
                Sb = sum(sum(dat1)~=0)*(mean(dat1,2)-grand_me)*(mean(dat1,2)-grand_me)'...
                    +sum(sum(dat2)~=0)*(mean(dat2,2)-grand_me)*(mean(dat2,2)-grand_me)';
                scr(k) = trace(Sb)/trace(Sw);       
            case 'DB-euc'
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);                
                d1 = sum(diag( (dat1-repmat(m1,1,n1))'*(dat1-repmat(m1,1,n1)) ));
                d2 = sum(diag( (dat2-repmat(m2,1,n2))'*(dat2-repmat(m2,1,n2)) ));
                d3 = (m1-m2)'*(m1-m2);
                scr(k) = -(1/n1*d1+1/n2*d2)/d3;
            case 'DB-mahal'
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                sigma = cov([dat1 dat2]');
                d1 = sum(diag( (dat1-repmat(m1,1,n1))'*1/sigma*(dat1-repmat(m1,1,n1)) ));
                d2 = sum(diag( (dat2-repmat(m2,1,n2))'*1/sigma*(dat2-repmat(m2,1,n2)) ));
                d3 = (m1-m2)'*1/sigma*(m1-m2);
                scr(k) = -(1/n1*d1+1/n2*d2)/d3;
            case 'MR-euc'
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);                
                d1 = sum(diag( (dat1-repmat(m1,1,n1))'*(dat1-repmat(m1,1,n1)) ));
                d2 = sum(diag( (dat2-repmat(m2,1,n2))'*(dat2-repmat(m2,1,n2)) ));
                d3 = (m1-m2)'*(m1-m2);
                scr(k) = -(d1+d2)/d3;
            case 'MR-mahal'
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                sigma = cov([dat1 dat2]');
                d1 = sum(diag( (dat1-repmat(m1,1,n1))'*1/sigma*(dat1-repmat(m1,1,n1)) ));
                d2 = sum(diag( (dat2-repmat(m2,1,n2))'*1/sigma*(dat2-repmat(m2,1,n2)) ));
                d3 = (m1-m2)'*1/sigma*(m1-m2);
                scr(k) = -(d1+d2)/d3;    
            case 'FR'
                grand_me = mean([dat1 dat2],2);
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                Sw = (n1-1)*cov(dat1')+(n2-1)*cov(dat2');
                Sb = n1*(m1-grand_me)*(m1-grand_me)' + n2*(m2-grand_me)*(m2-grand_me)';
                St = Sw+Sb;
                scr(k) = det(St+eps)/det(Sw+eps); 
            case 'Dunn-euc'
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                d12 = min(sum( ( repmat(dat1,1,n2) - dat2(:,reshape(repmat(1:n2,n1,1),1,n1*n2)) ).^2 , 1 ));
                diam1 = 2*1/n1*sum(sum( (dat1-repmat(m1,1,n1)).^2 ));
                diam2 = 2*1/n2*sum(sum( (dat2-repmat(m2,1,n2)).^2 ));
                scr(k) = d12/max(diam1,diam2);
            case 'Dunn-mahal'
                n1=size(dat1,2); n2=size(dat2,2);
                m1=mean(dat1,2); m2=mean(dat2,2);
                sigma = cov([dat1 dat2]');
                tmp = repmat(dat1,1,n2) - dat2(:,reshape(repmat(1:n2,n1,1),1,n1*n2));
                d12 = min(diag(tmp'*1/sigma*tmp));
                diam1 = 2*1/n1*sum(diag( (dat1-repmat(m1,1,n1))'*1/sigma*(dat1-repmat(m1,1,n1)) ));
                diam2 = 2*1/n2*sum(diag( (dat2-repmat(m2,1,n2))'*1/sigma*(dat2-repmat(m2,1,n2)) ));
                scr(k) = d12/max(diam1,diam2);
        end
        % penalty        
        overlap = sum(sum(sign(pitch1),1)>=2) + sum(sum(sign(pitch2),1)>=2);
        switch penType
            case 'sigmoid1'
                x0=round(p*nFrame); y0 = .99; x1=1/2*round(p*nFrame); y1=.5;
                a=(log(1/y0-1)-log(1/y1-1))/(x0-x1);   
                b=x0-log(1/y0-1)/a;                
                pen(k) = 1./(1+exp(a*(overlap-b)));
            case 'sigmoid2'
                x0=round(p*nFrame); y0 = .99; x1=3/4*round(p*nFrame); y1=.5;
                a=(log(1/y0-1)-log(1/y1-1))/(x0-x1);  b=x0-log(1/y0-1)/a;
                pen(k) = 1./(1+exp(a*(overlap-b)));
            case 'sigmoid3'
                x0=round(p*nFrame); y0 = .99; x1=7/8*round(p*nFrame); y1=.5;
                a=(log(1/y0-1)-log(1/y1-1))/(x0-x1);  b=x0-log(1/y0-1)/a;
                pen(k) = 1./(1+exp(a*(overlap-b)));
            case 'tansig'
                x0=round(p*nFrame); y0=.99;
                a = log(2/(1+y0)-1)/x0;
                pen(k) = 2./(1+exp(a*overlap))-1;
            case 'linear'            
                x0=round(p*nFrame);
                if overlap>=x0
                    pen(k) = 1;
                else
                    pen(k) = 1/x0*overlap;
                end
            case 'step'
                x0=round(p*nFrame);
                if overlap>=x0
                    pen(k) = 1;
                else
                    pen(k) = 0;
                end
        end
    else
        scr(k) = -inf;        
    end
end
normScr = lambda*(scr-min(scr))/(max(scr)-min(scr)+eps)-(1-lambda)*pen;
fprintf('\n');