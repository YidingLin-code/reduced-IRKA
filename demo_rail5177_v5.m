%Created by Yiding Lin, date: 20260212
%Reference: Yiding Lin and Valeria Simoncini, A reduced-IRKA method for
%large-scale $\mathcal{H}_2$-optimal model order reduction, Journal of Scientific Computing, 2026

%The R-IRKA algorithm frequently makes use of the classical IRKA algorithm.
%We thank Serkan Gugercin and Zlatko Drma\v{c} for making their Matlab
%codes of IRKA available, and Chris Beattie for an insightful conversation on IRKA.

clc
clear
close all
format compact
format short e

load('rail5177.mat');
%load('rail20209.mat');


r=11;%order

tol_IRKA=1e-8; % tolerance for IRKA
outer_tol=1e-8;% tolerance for R-IRKA or TR-IRKA
inner_tol=5e-9;% tolerance for inner IRKA in R-IRKA or TR-IRKA

n=size(A,1);
%%%%%%%%%%%%%
fprintf('~~~~~~~~~R_IRKA~~~~~~~~~\n');
[accum_iter_number,accum_error_vector,accum_shifts,accum_Btd,accum_Ctd,accum_V,accum_W,accum_inner_its]=R_IRKA_MIMO_v3(A,B,C,E,r,outer_tol,inner_tol);

R_IRKA_iteration_number=accum_iter_number

%%%%%%%%%%%%%
fprintf('~~~~~~~~~TR_IRKA~~~~~~~~~\n');
[DRKSM_iter_number,DRKSM_error_vector,DRKSM_shifts,DRKSM_Btd,DRKSM_Ctd,DRKSM_V,DRKSM_W,DRKSM_inner_its]=TR_IRKA_MIMO_v3(A,B,C,E,r,outer_tol,inner_tol);

TR_IRKA_iteration_number=DRKSM_iter_number
%%%%%%%%%%%%%
fprintf('~~~~~~~~~IRKA~~~~~~~~~\n');
[IRKA_iter_num,IRKA_error_vector,IRKA_shifts,IRKA_Btd,IRKA_Ctd,IRKA_V,IRKA_W]=IRKA_MIMO_v5(A,B,C,E,r,tol_IRKA);

IRKA_iteration_number=IRKA_iter_num

%%%%%%%%%%%%% error  figure
gcf_error=figure();
semilogy(IRKA_error_vector,'r*-')
hold on
semilogy(accum_error_vector,'bo-')
hold on
semilogy(DRKSM_error_vector,'kp-')
legend('IRKA','R-IRKA','TR-IRKA')
title('rail5177')
xlabel('Iteration number $$k$$','Interpreter','latex','FontSize',14);
ylabel('$$\chi^{\rm rel}_k$$','Interpreter','latex','FontSize',14);

%%%%%%%%%%%%%%%% shifts figure
gcf_ritz=figure();
plot(real(IRKA_shifts),imag(IRKA_shifts),'r*')
hold on
plot(real(accum_shifts),imag(accum_shifts),'bo')
hold on
plot(real(DRKSM_shifts),imag(DRKSM_shifts),'kp')
legend('IRKA','R-IRKA','TR-IRKA')
title('rail5177')
xlabel('Real part','Interpreter','latex','FontSize',14);
ylabel('Imaginary  part','Interpreter','latex','FontSize',14);

