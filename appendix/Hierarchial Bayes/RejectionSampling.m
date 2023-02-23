function samples = RejectionSampling(f, M, N, a, b)

batchsize = round(2 * N);

samples = zeros(N, 1);
samples_reached = 0;

while samples_reached < N
      samples_needed = N - samples_reached;
      cand = gamrnd(a, b, batchsize, 1);
      Us = rand(batchsize, 1);
      acc = Us < (f(cand)/M);
      samples_new = sum(acc);
      if samples_new == 0
          continue;
      end
      cand_acc = cand(acc);
      if sum(acc) >  samples_needed
         cand_acc = cand_acc(1:(samples_needed)); 
         samples_new = samples_needed;
      end
      samples((samples_reached+1):(samples_reached + samples_new)) = cand_acc; 
      samples_reached = samples_reached + samples_new;
      
    
end






