  module acc
  contains        
   subroutine compute_acc(eps,sigma,natoms,coord, mass, d, acceleration, f)
     implicit none
     integer, intent(in)          :: natoms
     double precision, intent(in) :: coord(natoms,3)
     double precision, intent(in) :: mass(natoms)
     double precision, intent(in) :: d(natoms,natoms)
     double precision, intent(out) :: acceleration(natoms,3)
     integer                      :: i, j, k
     double precision :: acc
     double precision, intent(in) :: eps, sigma
     double precision, intent(out) :: f(natoms,natoms)    ! force(r_ij)
      
     do i=1,natoms
        do k=1,3
          acc=0.d0 
           do j=1,natoms
             if(j.ne.i) then
             f(i,j)=force(d(i,j),eps,sigma)     ! in this way i store the info in an array
             acc=acc+f(i,j)*((coord(i,k)-coord(j,k))/d(i,j))
             endif
           enddo
          acceleration(i,k)=(-1.d0/mass(i))*acc
        enddo
     enddo

    contains
      double precision function force(r,eps,sigma)
      implicit none
      double precision, intent(in) :: r, eps, sigma
      double precision :: rap
      rap=(sigma/r)**6
      force=24*(eps/r)*(rap-2*rap**2)
      end function force
   end subroutine compute_acc
  end module acc


program main
   use acc     
   implicit none
   integer i, j, k, natoms
   double precision, allocatable :: coord(:,:), mass(:), d(:,:), vel(:,:)
   integer, parameter :: col=3
   double precision sums, eps, sigma, v, Vtot, Ttot, E, T, Etot
   double precision, allocatable :: acceleration(:,:),f(:,:)


   open(10, file='dat.inp')
   open(15, file='dist.out')
   open(20, file='pot.inp')

   read(10,*) natoms
   read(10,*) 
   
   allocate(coord(natoms,col),mass(natoms),d(natoms,natoms),vel(natoms,col))
   allocate(f(natoms,natoms))                      ! array of dV for each d(i,j)
   allocate(acceleration(natoms,col))

   do i=1,natoms
      read(10,*) (coord(i,j),j=1,col), mass(i)
      write(6,*) (coord(i,j),j=1,col), mass(i)
   enddo
   
   do i=1,natoms
     do j=i+1,natoms
       sums=0.0
         do k=1,col
           sums=sums+(coord(i,k)-coord(j,k))**2
         enddo
       d(i,j)=sqrt(sums)
       d(j,i)=d(i,j)
     enddo
       d(i,i)=0
     write(15,*) (d(i,j),j=1,natoms)
   enddo

  read(20,*) eps, sigma 

  Vtot=0.0

  do i=1,natoms
    do j=i+1,natoms
       v=4*eps*((sigma/d(i,j))**12-(sigma/d(i,j))**6)
       write(6,*) i, ",", j, "=", v  
       Vtot=Vtot+v
    enddo
  enddo
  
  write(6,*) "total LJ potential=", Vtot

  vel=0.0
  write(6,*) "initial velocities matrix"
  
  do i=1,natoms
   write(6,*) vel(i,:)
  enddo

  Ttot=T(natoms,mass,vel)
! E=Vtot+Ttot
  Etot=E(Vtot,Ttot)
  write(6,*) "total energy = Vtot + Ttot =", Vtot, "+", Ttot, "=", Etot 


  call compute_acc(eps,sigma,natoms,coord, mass, d, acceleration, f)
  
  do i=1,2
   write(6,*)
  enddo

  write(6,*) "acceleration array"

  do i=1,natoms
     write(6,*) (acceleration(i,k),k=1,col)
  enddo

  end program main

  double precision function T(natoms, mass, vel)
    implicit none
    integer, parameter :: col=3
    integer, intent(in):: natoms
    double precision :: vel(natoms,col)
    double precision, intent(in) :: mass(natoms) 
    integer :: i, j
    double precision :: velocity

    T=0.0
    do i=1,natoms
      velocity=0.0
       do j=1,col
         velocity=velocity+(vel(i,j))**2
       enddo
      T=T+mass(i)*velocity
    enddo
    T=0.5*T
  end function T

  double precision function E(Vtot, Ttot)
    implicit none
    double precision :: Vtot, Ttot
    
    E=Vtot+Ttot
  end function E 


