program reg_simple
implicit real(a-h,o-z)
integer,parameter::nx=1,ny=1,nb=3
integer,parameter::nsmp=10
real(4)::xread(nx)
real(4)::xsmp(nx,nsmp)
real(4)::bsmp(nb,nsmp)
real(4)::bcovinv(nb,nb)
real(4)::ysmp(ny,nsmp)
real(4)::w(ny,nb)

bsmp=0.0

open(11,file="x.dat",convert="big_endian",form="unformatted",access="stream")
open(12,file="y.dat",convert="big_endian",form="unformatted",access="stream")
do ismp=1,nsmp
  read(11,iostat=ios) xread
  read(12,iostat=ios) yread
  write(*,*) ismp, xread, yread
  xsmp(:,ismp)=xread
  ysmp(:,ismp)=yread
end do
close(11)
close(12)

call basisf
bcovinv=matmul(bsmp,transpose(bsmp))

!!! Inverse bcovinv
!call mat_inv(nb, bcovinv, bcovinv)
call inv(nb, bcovinv)

w=matmul(ysmp,matmul(transpose(bsmp),bcovinv))
write(*,*) maxval(xsmp), minval(xsmp)
write(*,*) maxval(bsmp), minval(bsmp)
write(*,*) maxval(w), minval(w)

stop
contains
        subroutine basisf
                do ib=1,nb
                      bsmp(ib,:)=xsmp(1,:)**(ib-1)
                end do
                return
        end subroutine

end program
