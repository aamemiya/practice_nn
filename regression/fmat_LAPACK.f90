subroutine det(ndim,varray,vdet)

integer,intent(in)::ndim
real(4),intent(in)::varray(ndim,ndim)
double precision::darray(ndim,ndim)
integer::ipiv(ndim)
 
darray=dble(varray)
call dgetrf(ndim,ndim,darray,ndim,ipiv,info)
if (info /= 0)then
  write(*,*) "cannot into det"
  vdet=0.0
  return
end if

vdet=1.0
do j=1,ndim
  if (ipiv(j).ne.j)then
    vdet=vdet*(-1.0)*darray(j,j)
  else
    vdet=vdet*darray(j,j)
  end if
end do

return
end subroutine

subroutine inv(ndim,varray)

  integer,intent(in)::ndim
  real(4),intent(inout)::varray(ndim,ndim)
  double precision::darray(ndim,ndim)
  integer::ipiv(ndim)
  integer::lwork
  integer::iwork(64*ndim)

  darray=varray
  
  lwork=64*ndim
  
  call dgetrf(ndim,ndim,darray,ndim,ipiv,info)
if (info.eq.0)then
    call dgetri(ndim,darray,ndim,ipiv,iwork,lwork,info)
 else
    write(*,*) 'cannot into inverse..'
    stop
  end if

  varray=darray
  
  return
end subroutine inv
!=======================================================!

subroutine real_eigen(ndim,vmat_in,varray_out_real,varray_out_aimag) !,vmat_out)
!  use func_sort
  integer,intent(in)::ndim
  real(4),intent(in)::vmat_in(ndim,ndim)
  real(4),intent(out)::varray_out_real(ndim)
  real(4),intent(out)::varray_out_aimag(ndim)
  !  real(4),intent(out)::vmat_out(ndim,ndim)


  complex::cmat_in(ndim,ndim)
  complex::cmat_out(ndim,ndim)
  complex::carray_out(ndim)

  complex::cdummy1(ndim,ndim)
  complex::cdummy2(ndim,ndim)
  
  complex,allocatable::cwork(:)
  complex::cwork_test(1)
  real(4)::rwork(2*ndim)  
  
  complex,parameter::c_ur=(1.0,0.0),c_ui=(0.0,1.0)

  character*1::cn='N'
  character*1::cv='V'

cmat_in = c_ur * vmat_in

!!! initialize =======================!
  call cgeev 	( cn, & !!! JOBVL
		  cv, & !!! JOBVR
		  ndim,  & !!! N
	       cmat_in,   & !!! real array
		  ndim,  & !!! LDA
	       carray_out,  & !!! WORK real
                 cdummy1,  & !!! WORK real
                 ndim,  & !!! LDVL
                 cdummy2,  & !!! WORK real
                 ndim,  & !!! LDVR
                 cwork_test,  & !!! WORK comp
               -1,  & !!! WORK integer dim
               rwork,  & !!! WORK real
		 info    & 
  	) 	

nwork=int(real(cwork_test(1)))
allocate(cwork(nwork))

!!=====================================!
  call cgeev 	( 'N', & !!! JOBVL
		  'V', & !!! JOBVR
		  ndim,  & !!! N
	       cmat_in,   & !!! real array
		  ndim,  & !!! LDA
	       carray_out,  & !!! WORK real
                 cdummy1,  & !!! WORK real
                 ndim,  & !!! LDVL
                 cdummy2,  & !!! WORK real
                 ndim,  & !!! LDVR
                 cwork,  & !!! WORK comp
               nwork,  & !!! WORK integer dim
               rwork,  & !!! WORK real
		 info    & 
  	) 	
  

  varray_out_real=real(carray_out)
  varray_out_aimag=aimag(carray_out)

 call rsort2(varray_out_real,varray_out_aimag,ndim)

 varray_out_real(:)=varray_out_real(ndim:1:-1)
 varray_out_aimag(:)=varray_out_aimag(ndim:1:-1)
 
  return
end subroutine real_eigen
!=======================================================!
subroutine real_eigen_mat_left(ndim,vmat_in,vmat_value,vmat_vec)
  integer,intent(in)::ndim
  real(4),intent(in)::vmat_in(ndim,ndim)
  real(4),intent(out)::vmat_value(ndim,ndim)
  real(4),intent(out)::vmat_vec(ndim,ndim)
  real(4)::varray_out_real(ndim)
  real(4)::varray_out_aimag(ndim)

  complex::cmat_in(ndim,ndim)
  complex::cmat_out(ndim,ndim)
  complex::carray_out(ndim)

  complex::cdummy1(ndim,ndim)
  complex::cdummy2(ndim,ndim)
  
  complex,allocatable::cwork(:)
  complex::cwork_test(1)
  real(4)::rwork(2*ndim)  
  
  complex,parameter::c_ur=(1.0,0.0),c_ui=(0.0,1.0)

  character*1::cn='N'
  character*1::cv='V'

 
cmat_in = c_ur * vmat_in

!!! initialize =======================!
  call cgeev 	( cv, & !!! JOBVL
		  cn, & !!! JOBVR
		  ndim,  & !!! N
	       cmat_in,   & !!! real array
		  ndim,  & !!! LDA
	       carray_out,  & !!! WORK real
                 cdummy1,  & !!! WORK real
                 ndim,  & !!! LDVL
                 cdummy2,  & !!! WORK real
                 ndim,  & !!! LDVR
                 cwork_test,  & !!! WORK comp
               -1,  & !!! WORK integer dim
               rwork,  & !!! WORK real
		 info    & 
  	) 	

nwork=int(real(cwork_test(1)))
allocate(cwork(nwork))

!!=====================================!
  call cgeev 	( 'V', & !!! JOBVL
		  'N', & !!! JOBVR
		  ndim,  & !!! N
	       cmat_in,   & !!! real array
		  ndim,  & !!! LDA
	       carray_out,  & !!! WORK real
               cmat_out,  & !!! WORK real
                 ndim,  & !!! LDVL
                 cdummy2,  & !!! WORK real
                 ndim,  & !!! LDVR
                 cwork,  & !!! WORK comp
               nwork,  & !!! WORK integer dim
               rwork,  & !!! WORK real
		 info    & 
  	) 	
  

  varray_out_real=real(carray_out)
  varray_out_aimag=aimag(carray_out)
  vmat_vec=real(cmat_out)


  do idim=1,ndim
   vmat_value(idim,idim)=varray_out_real(idim)
  end do 

 call rsort2_ex(varray_out_real,varray_out_aimag,vmat_vec,ndim)
  
 varray_out_real(:)=varray_out_real(ndim:1:-1)
 varray_out_aimag(:)=varray_out_aimag(ndim:1:-1)

do idim=1,ndim
   !vmat_vec(1:ndim,idim)= vmat_vec(ndim:1:-1,idim) !!! Error
   vmat_vec(idim,1:ndim)= vmat_vec(idim,ndim:1:-1)
end do

  vmat_value=0.0
 do idim=1,ndim
  vmat_value(idim,idim)=varray_out_real(idim)
 end do 

  return
end subroutine real_eigen_mat_left
!=======================================================!
subroutine real_eigen_mat_right(ndim,vmat_in,vmat_value,vmat_vec)
  integer,intent(in)::ndim
  real(4),intent(in)::vmat_in(ndim,ndim)
  real(4),intent(out)::vmat_value(ndim,ndim)
  real(4),intent(out)::vmat_vec(ndim,ndim)
  real(4)::varray_out_real(ndim)
  real(4)::varray_out_aimag(ndim)

  complex::cmat_in(ndim,ndim)
  complex::cmat_out(ndim,ndim)
  complex::carray_out(ndim)

  complex::cdummy1(ndim,ndim)
  complex::cdummy2(ndim,ndim)
  
  complex,allocatable::cwork(:)
  complex::cwork_test(1)
  real(4)::rwork(2*ndim)  
  
  complex,parameter::c_ur=(1.0,0.0),c_ui=(0.0,1.0)

  character*1::cn='N'
  character*1::cv='V'

 
cmat_in = c_ur * vmat_in

!!! initialize =======================!
  call cgeev 	( cv, & !!! JOBVL
		  cn, & !!! JOBVR
		  ndim,  & !!! N
	       cmat_in,   & !!! real array
		  ndim,  & !!! LDA
	       carray_out,  & !!! WORK real
                 cdummy1,  & !!! WORK real
                 ndim,  & !!! LDVL
                 cdummy2,  & !!! WORK real
                 ndim,  & !!! LDVR
                 cwork_test,  & !!! WORK comp
               -1,  & !!! WORK integer dim
               rwork,  & !!! WORK real
		 info    & 
  	) 	

nwork=int(real(cwork_test(1)))
allocate(cwork(nwork))

!!=====================================!
  call cgeev 	( 'N', & !!! JOBVL
		  'V', & !!! JOBVR
		  ndim,  & !!! N
	       cmat_in,   & !!! real array
		  ndim,  & !!! LDA
	       carray_out,  & !!! WORK real
                 cdummy2,  & !!! WORK real
                 ndim,  & !!! LDVL
               cmat_out,  & !!! WORK real
                 ndim,  & !!! LDVR
                 cwork,  & !!! WORK comp
               nwork,  & !!! WORK integer dim
               rwork,  & !!! WORK real
		 info    & 
  	) 	
  

  varray_out_real=real(carray_out)
  varray_out_aimag=aimag(carray_out)
  vmat_vec=real(cmat_out)

  do idim=1,ndim
   vmat_value(idim,idim)=varray_out_real(idim)
  end do 

 call rsort2_ex(varray_out_real,varray_out_aimag,vmat_vec,ndim)
  
 varray_out_real(:)=varray_out_real(ndim:1:-1)
 varray_out_aimag(:)=varray_out_aimag(ndim:1:-1)

do idim=1,ndim
   !vmat_vec(1:ndim,idim)= vmat_vec(ndim:1:-1,idim) !!! Error
   vmat_vec(idim,1:ndim)= vmat_vec(idim,ndim:1:-1)
end do

  vmat_value=0.0
 do idim=1,ndim
  vmat_value(idim,idim)=varray_out_real(idim)
 end do 

  return
end subroutine real_eigen_mat_right
!=======================================================!
subroutine sqrt_mat(ndim,vmat_in,vmat_out)

integer,intent(in)  :: ndim
real(4),intent(in)  :: vmat_in  (ndim,ndim)
real(4),intent(out) :: vmat_out (ndim,ndim)

real(4) :: vmat_D (ndim,ndim)
real(4) :: vmat_U (ndim,ndim)

call real_eigen_mat_left(ndim,vmat_in,vmat_D,vmat_U)
do idim=1,ndim
 vmat_D(idim,idim)=sqrt(max(vmat_D(idim,idim),0.0))
end do
vmat_out = matmul(matmul(vmat_U,vmat_D),transpose(vmat_U))

return
end subroutine sqrt_mat
!=======================================================!
