%% read matrices
qust_B=read_dumped('../data/qust_B.dat', 961, 19);
qust_A=read_dumped('../data/qust_A.dat', 24, 19);
qust_A_plus=read_dumped('../data/qust_A_plus.dat', 19, 24);
qust_interp=read_dumped('../data/qust_interp.dat', 961, 24);

qugrs_B=read_dumped('../data/qugrs_B.dat', 961, 19);
qugrs_A=read_dumped('../data/qugrs_A.dat', 24, 19);
qugrs_A_plus=read_dumped('../data/qugrs_A_plus.dat', 19, 24);
qugrs_interp=read_dumped('../data/qugrs_interp.dat', 961, 24);

% plot differences
figure; imagesc(abs(qugrs_A - qust_A)); colorbar
figure; imagesc(abs(qugrs_B - qust_B)); colorbar
figure; imagesc(abs(qugrs_A_plus - qust_A_plus)); colorbar
figure; imagesc(abs(qugrs_interp - qust_interp)); colorbar

%% compare pseudoinverse to matlab pseudoinverse
figure; imagesc(abs(qust_A_plus - pinv(qust_A))); colorbar
figure; imagesc(abs(qugrs_A_plus - pinv(qugrs_A))); colorbar

%% read SVD matrices
qust_U=read_dumped('../data/qust_U.dat', 24, 19);
qust_S=diag(read_dumped('../data/qust_S.dat', 19, 1));
qust_V=read_dumped('../data/qust_V.dat', 19, 19);
qust_S_plus=read_dumped('../data/qust_S_plus.dat', 19, 19);
qust_U_times_S_plus=read_dumped('../data/qust_U_times_S_plus.dat', 24, 19);


qugrs_U=read_dumped('../data/qugrs_U.dat', 24, 19);
qugrs_S=diag(read_dumped('../data/qugrs_S.dat', 19, 1));
qugrs_V=read_dumped('../data/qugrs_V.dat', 19, 19);
qugrs_S_plus=read_dumped('../data/qugrs_S_plus.dat', 19, 19);
qugrs_U_times_S_plus=read_dumped('../data/qugrs_U_times_S_plus.dat', 24, 19);

figure; imagesc(abs(qugrs_U - qust_U)); colorbar
figure; imagesc(abs(qugrs_S - qust_S)); colorbar
figure; imagesc(abs(qugrs_V - qust_V)); colorbar

%% check multiplications

qugrs_S_plus_m = diag(1./diag(qugrs_S));
qugrs_U_times_S_plus_m = qugrs_U * qugrs_S_plus;
qugrs_A_plus = qugrs_V * qugrs_U_times_S_plus';

figure; imagesc(abs(qugrs_S_plus - qugrs_S_plus_m)); colorbar
figure; imagesc(abs(qugrs_U_times_S_plus - qugrs_U_times_S_plus_m)); colorbar
figure; imagesc(abs(qugrs_A_plus - qugrs_A_plus_m)); colorbar
