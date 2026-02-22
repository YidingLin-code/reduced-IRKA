function [i_IRKA,error,new_shifts,Br,Cr,Vr,Wr]=IRKA_MIMO_v5(A,B,C,E,r,tol_IRKA)
%This is the classical IRKA algorithm.
%Reference: S. Gugercin, A. C. Antoulas, and C. A. Beattie, H2 model reduction
%for large-scale linear dynamical systems, SIAM J. Matrix Anal. Appl., 30
%(2008),pp. 609–638


%opts.tol=1e-1;
opts.disp = 0;
rand('state',0);
[Vr,Lr]=eigs(A,E,r,0,opts);

Vr = orth([real(Vr),imag(Vr)]);
Vr=Vr(:,1:r);
Wr = Vr;

old_shifts=ones(r,1);

for i_IRKA=1:300
    
    %i_IRKA   
    Er = Wr'*E*Vr;
    Ar = Wr'*A*Vr; Br = Wr'*B; Cr = C*Vr;
    [X,ritz_values]=eig(full(Ar),full(Er));
    EX = X\(Er*X);
    Br = EX\(X\Br);
    Cr = Cr*X;
    
    tem_shifts = -diag(ritz_values);
    [new_shifts,ord] = sort(tem_shifts);
    Br = Br(ord,:);
    Cr = Cr(:,ord);
    
    error(i_IRKA)=norm(new_shifts-old_shifts)/norm(new_shifts);
    
    IRKA_error=error(i_IRKA);
    
    if error(i_IRKA)<tol_IRKA
        break;
    end
    old_shifts=new_shifts;
    
    [Vr,Wr] =krylov_red_dir_v3(A,B,C,E,new_shifts,Br,Cr);
    
end

[Vr,Wr] =krylov_red_dir_v3(A,B,C,E,new_shifts,Br,Cr);

%final_error=error_shifts(jj);
%semilogy(error)