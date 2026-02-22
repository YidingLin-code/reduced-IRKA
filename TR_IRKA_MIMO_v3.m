function [iter_number,error_shifts,new_shifts,new_Br,new_Cr,Vr,Wr,inner_its_vector]=TR_IRKA_MIMO_v3(A,B,C,E,r,outer_tol,inner_tol)
%Created by Yiding Lin, date: 20260212
%Reference: A reduced-IRKA method for large-scale $\mathcal{H}_2$-optimal
%model order reduction, Yiding Lin and Valeria Simoncini
%The R-IRKA algorithm frequently makes use of the classical IRKA algorithm.
%We thank Serkan Gugercin and Zlatko Drma\v{c} for making their Matlab
%codes of IRKA available, and Chris Beattie for an insightful conversation on IRKA.

%This is the main function of TR-IRKA(Truncated R-IRKA).

%Input: r, order; 
%outer_tol, tolerance for the outer iteration.
%inner_tol, tolerance for the inner iteration.

%Output:new_shifts, final H2 optimal shifts;
%iter_number, outer iteration number;
%inner_its_vector,  this vector stores every inner iteration number;
%Wr,Vr, orthonormal basis of H2 optimal subspaces; 
%new_Br,new_Cr, tangent directions.
%error_shifts, this vector stores $$\chi_k$$ for every outer iteration.


%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
%FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
%COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
%IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
%CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

r_two = 2*r;
%opts.disp = 1e-1;
%opts.tol = 1e-1;
opts.disp = 0;
rand('state',0);
[Vr,Lr]=eigs(A,E,r_two,0,opts);
Vr = orth([real(Vr),imag(Vr)]);
Vr =Vr (:,1:r_two);
Wr = Vr;

old_shifts=ones(r,1);
shifts=ones(r_two,1);
all_shifts=shifts;

n=size(A,1);
Vr_tem=zeros(n,r);
Wr_tem=zeros(n,r);


for jj=1:30
   % outer_jj=jj
    
    Er = Wr'*E*Vr;
    Ar = Wr'*A*Vr; Br = Wr'*B; Cr = C*Vr;
    %ritz_value=eig(full(Ar),full(Er));
    
    if jj<1.5
        [new_shifts,new_Br,new_Cr,new_Vr,new_Wr,inner_its]=inner_IRKA_real_shifts_its_it1_v3(Ar,Br,Cr,Er,r,inner_tol);
    else
        [new_shifts,new_Br,new_Cr,new_Vr,new_Wr,inner_its]=inner_IRKA_real_shifts_its_it2_v3(Ar,Br,Cr,Er,r,new_shifts,new_Br,new_Cr,inner_tol);
    end
    
    inner_its_vector(jj)=inner_its;
    
    %new_shifts
    %figure()
    %plot(real(new_shifts), imag(new_shifts),'*')
    error_shifts(jj)=norm(new_shifts-old_shifts)/norm(new_shifts);
    %current_outer_error= error_shifts(jj)
    old_shifts=new_shifts;
    all_shifts=[all_shifts;new_shifts];   
    
    if error_shifts(jj)<outer_tol
        break;
    end  
    
    if  jj>1.5
        Vr=Vr(:,1:2*r);%truncate V: Only the front columns are retained.
        Wr=Wr(:,1:2*r);
    end
    
    
    [Vr_tem,Wr_tem] =krylov_red_dir_no_qr_v3(A,B,C,E,new_shifts,new_Br,new_Cr);
    
    %     [Vr,~] = qr(full(Vr_tem),0);
    %     [Wr,~] = qr(full(Wr_tem),0);
    
    [Vr,~] = qr([Vr_tem,Vr],0);%Note that the latest Vr_tem is placed in front of the old Vr.
    [Wr,~] = qr([Wr_tem,Wr],0);
    
    %[Vr,~] = qr([Vr,Vr_tem],0);
    % [Wr,~] = qr([Wr,Wr_tem],0);
    
end
%jj_end=jj

%The orthonormal basis is to be obtained after termination of the iteration process.
%[Vr_tem,Wr_tem] =krylov_red_dir_v3(A,B,C,E,new_shifts,new_Br,new_Cr);
[Vr_tem,Wr_tem] =krylov_red_dir_no_qr_v3(A,B,C,E,new_shifts,new_Br,new_Cr);
[Vr,~] = qr([Vr_tem],0);
[Wr,~] = qr([Wr_tem],0);

ind_real = find(abs(imag(new_shifts)./abs(new_shifts))<outer_tol);
new_shifts(ind_real)=real(new_shifts(ind_real));

iter_number=jj;
%final_error=error_shifts(jj);
%semilogy(error_shifts)
