function jobepilog(pathto)
% Epilog script to compute voltaes and spikes

batch_mv = load_raw_batch_agnostic(pathto);
fprintf('Cell batch loaded\n');
