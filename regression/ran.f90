!============================================
! Box-muller method 
!============================================
real function seiki_ran(indx)
  integer::indx
  real(4)::va,vb
  real(4)::pi

  pi=3.14159263

  call system_clock(count=icount)
  indx_temp=indx+icount

  va = rngu0(indx_temp)
  vb = rngu0(indx_temp)
  seiki_ran = sqrt(-2.0 * log(va)) * sin(2.0*pi*vb)

end function seiki_ran
!============================================
real function range_ran(indx,vmin,vmax)
  integer::indx
  real(4)::vmin,vmax

  call system_clock(count=icount)
  indx_temp=indx+icount
  range_ran = vmin + (vmax-vmin) * rngu0(indx_temp)
  
end function range_ran
!============================================
