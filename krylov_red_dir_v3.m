function [Vr,Wr] = krylov_red_dir_v3(AA,BB,CC,EE,S2,Br,Cr)
% Created by S. Gugercin

ind_real = find(abs(imag(S2)./abs(S2))<1e-8);
ind_imag = find(abs(imag(S2)./abs(S2))>=1e-8);
S2real = S2(ind_real);
if ~isempty(S2real)
    S2real = 0.5*(S2real+conj(S2real));
end

S2imag = S2(ind_imag);
Brreal = Br(ind_real,:);
Brreal = 1/2*(Brreal + conj(Brreal));
Brimag = Br(ind_imag,:);

Crreal = Cr(:,ind_real);
Crreal = 1/2*(Crreal + conj(Crreal));
Crimag = Cr(:,ind_imag);

l_real = length(S2real);
l_imag = length(S2imag);
S2imag = S2imag(1:2:l_imag-1);
Brimag = Brimag(1:2:l_imag-1,:);
Crimag = Crimag(:,1:2:l_imag-1);

l_imag_half = length(S2imag);

Vrimag = [];  Wrimag = [];
for i=1:l_imag_half
    Vrimag(:,i) = (S2imag(i)*EE-AA)\(BB*((Brimag(i,:)).'));
    Wrimag(:,i) = (S2imag(i)*EE'-AA')\(CC'*(Crimag(:,i)));
end

Vrreal = [];  Wrreal = [];
for i=1:l_real
    Vrreal(:,i) = (S2real(i)*EE-AA)\(BB*(Brreal(i,:))');
    Wrreal(:,i) = (S2real(i)*EE-AA)'\(CC'*Crreal(:,i));
end
[Vr] = ([Vrreal real(Vrimag) imag(Vrimag)]);
[Wr] = ([Wrreal real(Wrimag) imag(Wrimag)]);
[Vr,Rr] = qr(full(Vr),0);
[Wr,Rr] = qr(full(Wr),0);