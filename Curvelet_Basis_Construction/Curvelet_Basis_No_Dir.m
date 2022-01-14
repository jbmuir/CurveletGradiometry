function [G_mat, scales]=Curvelet_Basis_No_Dir(nx,ny)

  nxnew=nx*2;
  nynew=ny*2;
  ilet=0;
  space_xz0=zeros(ny,nx);
  cc0=fdct_wrapping(space_xz0,1,2, floor(log2(min(nx,ny))-1));
  ns=length(cc0);
  na=length(cc0{2});

  for j=1:1
      for l=1:length(cc0{j})
          [A,B]=size(cc0{j}{l});
          for k1=1:A
              for k2=1:B
                  ilet=ilet+1;
                  cc0_jkl=cc0;
                  cc0_jkl{j}{l}(k1,k2)=1;
                  curvelet_xz0=ifdct_wrapping(cc0_jkl,1,ny,nx);
                  G_mat(:,ilet)=reshape(curvelet_xz0,[],1);
                  scales(ilet) = j;
              end
          end
      end
  end

  space_xz=zeros(nynew,nxnew);

  for j=2:length(cc0)-1
      if(mod(j,2)==0)
          cc=fdct_wrapping(space_xz,1,2,ns+1,na/2);
      else
          cc=fdct_wrapping(space_xz,1,2,ns+1,na);
      end
      for l=1:length(cc{j+1})
          [A,B]=size(cc{j+1}{l});
          cc_jkl=cc;
          cc_jkl{j+1}{l}(floor(A/2),floor(B/2))=1;
          curvelet_xz=ifdct_wrapping(cc_jkl,1,nynew,nxnew);
          temp1(1:ny,1:nx)=curvelet_xz(floor(ny/2):floor(ny/2)+ny-1,floor(nx/2):floor(nx/2)+nx-1);
          maxnorm2=norm(reshape(temp1,[],1))^2;
          for k1=1:A
              for k2=1:B
                  cc_jkl=cc;
                  cc_jkl{j+1}{l}(k1,k2)=1;
                  curvelet_xz=ifdct_wrapping(cc_jkl,1,nynew,nxnew);
                  temp1(1:ny,1:nx)=curvelet_xz(floor(ny/2):floor(ny/2)+ny-1,floor(nx/2):floor(nx/2)+nx-1);
                  nownorm2=norm(reshape(temp1,[],1))^2;
                  if(nownorm2/maxnorm2>0.2)
                      ilet=ilet+1;
                      G_mat(:,ilet)=reshape(temp1,[],1);
                      scales(ilet) = j;
                  end
              end
          end
      end
      clear cc;
      clear cc_jkl;
  end

end
