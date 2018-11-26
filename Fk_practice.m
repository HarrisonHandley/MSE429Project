% Left hand side (Li^2)
%L_sqr = Lnorm^2; 


syms Px Py Pz Ppx Ppy Pbx Pby alpha_i beta_i gamma_i 

R0_ee=[cosd(alpha_i)*cosd(beta_i), cosd(alpha_i)*sind(beta_i)*sind(gamma_i)-sind(alpha_i)*cosd(gamma_i), cosd(alpha_i)*sind(beta_i)*cosd(gamma_i)+sind(alpha_i)*sind(gamma_i)
           sind(alpha_i)*cosd(beta_i), sind(alpha_i)*sind(beta_i)*sind(gamma_i)+cosd(alpha_i)*cosd(gamma_i), sind(alpha_i)*sind(beta_i)*cosd(gamma_i)-cosd(alpha_i)*sind(gamma_i)
           -sind(beta_i)           , cosd(beta_i)*sind(gamma_i)                                    , cosd(beta_i)*cosd(gamma_i)];
       
R0_p0= R0_ee; 
P_b0_p0= [Px,Py,Pz]';
P_p0_pi= [Ppx,Ppy,0]'; 
P_b0_bi= [Pbx,Pby,0]';
P_bi_pi=simplify(sum(expand((P_b0_p0+R0_p0*P_p0_pi-P_b0_bi).^2)));

Jx = [diff(P_bi_pi,Px);diff(P_bi_pi,Py);diff(P_bi_pi,Pz);diff(P_bi_pi,alpha_i);diff(P_bi_pi,beta_i);diff(P_bi_pi,gamma_i)]';
%diff_px = diff(P_bi_pi,Px);
