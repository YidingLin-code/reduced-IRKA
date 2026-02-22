function[new_shifts,Br,Cr,Wr,Vr,its_number]=inner_IRKA_real_shifts_its_it2_v3(A,B,C,E,r,shifts,Br,Cr,inner_tol);
%This is the inner iteration (except for the first iteration) in R-IRKA or TR-IRKA.
%It is accomplished by a classical IRKA algorithm. 
%The initial settings are obtained from the output of the preceding inner iteration.

%Input: r, order; inner_tol, tolerance for stopping IRKA iteration;
%shifts, initial H2 otpimal shifts;Br,Cr, initial tangent directions.

%Output:new_shifts, final H2 optimal shifts;
%its_number, iteration number;
%Wr,Vr, orthonormal basis of H2 optimal subspaces; 
%Br,Cr, tangent directions.

old_shifts=shifts;
[Vr,Wr] =krylov_red_dir_v3(A,B,C,E,shifts,Br,Cr);

for i_IRKA=1:300
    
 %   inner_its=i_IRKA
Er = Wr'*E*Vr;
%size(Wr),size(Vr)
Ar = Wr'*A*Vr; Br = Wr'*B; Cr = C*Vr;
[X,ritz_values]=eig(full(Ar),full(Er));
EX = X\(Er*X);
Br = EX\(X\Br);
Cr = Cr*X;
    
tem_shifts = -diag(ritz_values);
[new_shifts,ord] = sort(tem_shifts);
Br = Br(ord,:);
Cr = Cr(:,ord);

% new_shifts
error(i_IRKA)=norm(new_shifts-old_shifts)/norm(new_shifts);

if error(i_IRKA)<inner_tol
    break;
end
old_shifts=new_shifts;
%inner_error=error(i_IRKA)

[Vr,Wr] =krylov_red_dir_v3(A,B,C,E,new_shifts,Br,Cr);
    
end
its_number=i_IRKA;


%semilogy(error)