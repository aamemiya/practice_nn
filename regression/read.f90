program read

integer,parameter::ndim=1
integer,parameter::nsmp=10
real(4)::xread(ndim)
real(4)::xsmp(ndim,nsmp)

open(11,file="x.dat",convert="big_endian",form="unformatted",access="stream")
!open(11,file="x.dat",convert="big_endian",form="unformatted")
do ismp=1,nsmp
  read(11,iostat=ios) xread
  write(*,*) ismp, xread
  xsmp(:,ismp)=xread
end do 
close(11)

stop
end program
