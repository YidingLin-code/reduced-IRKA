function[new_shifts,Br,Cr,Wr,Vr,its_number]=inner_IRKA_real_shifts_its_it1_v3(A,B,C,E,r,inner_tol);
%This is the first inner iteration in R-IRKA or TR-IRKA.
%This inner iteration is accomplished by a classical IRKA algorithm.

%Input: r, order; inner_tol, tolerance for stopping IRKA iteration.

%Output:new_shifts, final H2 optimal shifts;
%its_number, iteration number;
%Wr,Vr, orthonormal basis of H2 optimal subspaces; 
%Br,Cr, tangent directions.

opts.tol=1e-1;
opts.disp = 0;
 rand('state',0);
[Vr,Lr]=eigs(A,E,r,0,opts);
Vr = orth([real(Vr),imag(Vr)]);
Vr=Vr(:,1:r);
Wr = Vr;

old_shifts=ones(r,1);

for i_IRKA=1:300
 % inner_its=i_IRKA
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

%new_shifts
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
