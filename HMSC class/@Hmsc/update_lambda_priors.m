function [psijh,delta] = update_lambda_priors(nf,nur,a1r,a2r,b1r,b2r,psijh,delta,lambda,factorcov)
nr = length(nf);
for i1=1:nr
   ns = size(lambda{i1},2);
   delta1=delta{i1};
   if factorcov(i1)
      tau = cumprod(delta1, 1);
      lambda1 = permute(lambda{i1},[2,1,3]);
      psijh1 = gamrnd(nur/2 + 0.5,1./(nur/2 + 0.5*bsxfun(@times,lambda1.^2,permute(tau,[3,1,2]))));
      mat = bsxfun(@times,psijh1,lambda1.^2);
      sumMat = permute(sum(mat),[2,3,1]);
      ad = a1r + 0.5*ns*nf(i1);
      tauh = bsxfun(@rdivide, tau, delta1(1,:));
      bd = b1r + 0.5*sum(tauh.*sumMat, 1);
      delta1(1,:) = gamrnd(ad,1./bd);
      for h = 2:nf(i1)
         ad = a2r + 0.5*ns*(nf(i1)-h+1);
         tauh = bsxfun(@rdivide, tau, delta1(h,:));
         bd = b1r + 0.5*sum(tauh(h:end,:).*sumMat(h:end,:), 1);
         delta1(h,:) = gamrnd(ad,1./bd);
         tau = cumprod(delta1, 1);
      end
      % Previos version which is likely not very correct
      %       tauh = cumprod(delta1);
      % 	if factorcov(i1)
      % 		tauh2 = permute(tauh,[1,3,2]);
      % 		lambda1 = permute(lambda{i1},[2,1,3]);
      % 		psijh1 = gamrnd(nur/2 + 0.5,1./(nur/2 + 0.5*bsxfun(@times,lambda1.^2,permute(tauh2,[2,1,3]))));
      % 		mat = bsxfun(@times,psijh1,lambda1.^2);
      % 		ad = a1r + 0.5*ns*nf(i1);
      % 		bd = b1r + 0.5*(1/delta1(1))*sum(tauh2.*permute( sum(mat),[2,1,3] ));
      % 		delta1(1,:) = gamrnd(ad,1/bd);
      % 		for h = 2:nf(i1)
      % 			ad = a2r + 0.5*ns*(nf(i1)-h+1);
      % 			bd = b2r + 0.5*(1/delta1(h))*sum(tauh2(h:end,:,:).*permute( sum(mat(:,h:end,:)),[2,1,3] ));
      % 			delta1(h,:) = gamrnd(ad,1/bd);
      % 			tauh = cumprod(delta1);
      % 			tauh2 = permute(tauh,[1,3,2]);
      % 		end
      % 	else
   else
      tauh = cumprod(delta1);
      lambda1=lambda{i1}';
      psijh1 = gamrnd(nur/2 + 0.5,1./(nur/2 + 0.5*bsxfun(@times,lambda1.^2,tauh')));
      mat = bsxfun(@times,psijh1,lambda1.^2);
      ad = a1r + 0.5*ns*nf(i1);
      bd = b1r + 0.5*(1/delta1(1))*sum(tauh.*sum (mat)');
      delta1(1) = gamrnd(ad,1/bd);
      for h = 2:nf(i1)
         ad = a2r + 0.5*ns*(nf(i1)-h+1);
         bd = b2r + 0.5*(1/delta1(h))*sum(tauh(h:end).*sum (mat(:,h:end))');
         delta1(h) = gamrnd(ad,1/bd);
         tauh = cumprod(delta1);
      end
   end
   psijh{i1}=psijh1;
   delta{i1}=delta1;
end