program main
   implicit none
   integer i, j, k, natoms
   double precision, allocatable :: coord(:,:), mass(:), d(:,:), vel(:,:)
   integer, parameter :: col=3
   double precision sums, eps, sigma, v, Vtot, Ttot, E, T

   open(10, file='dat.inp')
   open(15, file='dist.out')
   open(20, file='pot.inp')

   read(10,*) natoms
   read(10,*) 
   
   allocate(coord(natoms,col),mass(natoms),d(natoms,natoms),vel(natoms,col))

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
  write(6,*) vel
  Ttot=T(natoms,mass,vel)
  E=Vtot+Ttot
  write(6,*) "total energy = Vtot + T =", Vtot, "+", Ttot, "=", E 

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
         vel(i,j)=0.0
         velocity=velocity+(vel(i,j))**2
       enddo
      T=T+mass(i)*velocity
    enddo
    T=0.5*T
  end function T 
