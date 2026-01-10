!---------------------------------------------------------------------------------------!
!-Module to create explicit interface of allocatable dummy arguments:-------------------!
!-it contains the subroutine that reads coordinates and masses of atoms-----------------!                      
                                                                                        ! 
  module reading                                                                        !
  contains                                                                              !
    subroutine atom_coord(natoms, coord, mass, col)                                     !
       implicit none                                                                    !
       integer :: natoms, i, j, col                                                     !
       double precision, allocatable :: mass(:)                                         !
       double precision, allocatable :: coord(:,:)                                      !
                                                                                        !
       do i=1,natoms                                                                    !
         read(10,*) coord(i,1), coord(i,2), coord(i,3) , mass(i)                        !
         write(6,*) (coord(i,j),j=1,col), mass(i)                                       !
       enddo                                                                            !
                                                                                        !
    return                                                                              !
    end subroutine atom_coord                                                           !
  end module reading                                                                    !
                                                                                        !
!---------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------!
!-Module to create explicit interface of allocatable dummy arguments:-------------------!
!-it contains the subroutine that compute distances for each pair of atoms--------------!
                                                                                        !
  module distances                                                                      !
  contains                                                                              !
    subroutine compute_distances(natoms,coord, col, d)                                  !
       implicit none                                                                    !
       integer :: natoms, i, j, k, col                                                  !
       double precision, allocatable :: mass(:)                                         !
       double precision, allocatable :: coord(:,:)                                      !
       double precision, allocatable :: d(:,:)                                          !
       double precision sums                                                            !
                                                                                        !
       do i=1,natoms                                                                    !
         do j=i+1,natoms                                                                !
           sums=0.d0                                                                    !
             do k=1,col                                                                 !
               sums=sums+(coord(i,k)-coord(j,k))**2                                     !
             enddo                                                                      !
           d(i,j)=sqrt(sums)                                                            !
           d(j,i)=d(i,j)                                                                !
         enddo                                                                          !  
           d(i,i)=0.d0                                                                  !
       enddo                                                                            !
                                                                                        !
    return                                                                              !
    end subroutine compute_distances                                                    !
  end module distances                                                                  !
                                                                                        !
!---------------------------------------------------------------------------------------!


!---------------------------------------------------------------------------------------!
!-Module to create explicit interface of allocatable dummy arguments:-------------------!
!-it contains the subroutine that computes the acceleration for each atom---------------!

  module accel
  contains
   subroutine compute_acc(eps,sigma,natoms,coord, mass, d, acceleration,f)
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
     double precision, parameter :: dmin=1.0d-12


     do i=1,natoms
        do k=1,3
          acc=0.d0
           do j=1,natoms
             if(j.ne.i.and.d(i,j).gt.dmin) then
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
  end module accel






!---------------------------------------------------------------------------------------!











!-------------------------------------------------------------------------------------------------!                         
!-Molecular Dynamics Code-------------------------------------------------------------------------!                  
                                                                                                  !                            
  program main                                                                                    !                              
                                                                                                  !                                         
                                                                                                  !                           
                                                                                                  !                                        
!---------------------------------------------------------------------------------------!         !                          
!-declaration block---------------------------------------------------------------------!         !                       
                                                                                        !         !                   
  use reading                                                                           !         !               
  use distances
  use accel 
  
  
  implicit none                                                                         !         !                       
  integer read_natoms                                                                   !         !               
  integer unit, natoms, status                                                          !         !                  
  double precision, allocatable :: mass(:)                                              !         !                    
  double precision, allocatable :: coord(:,:)                                           !         !                
  double precision, allocatable :: d(:,:)
  double precision, allocatable :: vel(:,:)
  double precision, allocatable :: acceleration(:,:), f(:,:)
  double precision :: eps, sigma, V, LJtot, Ttot, T, E, Etot
  integer, parameter :: col=3                                                           !         !               
  integer i, j, k                                                                       !         !              
                                                                                        !         !                
!---------------------------------------------------------------------------------------!         !                       
                                                                                                  !              
!---------------------------------------------------------------------------------------!         !                    
!-opening of input and output file and reading of input variables-----------------------!         !              
                                                                                        !         !            
  open(unit=10, file='dat.inp')                                                         !         !                 
  open(15, file='coord.out')                                                            !         !              
  open(20, file='distances.out')  
               
  unit=10                                                                               !         !              
                                                                                        !         !            
   natoms=read_natoms(unit)                                                             !         !     
   write(6,*) natoms   ! check if natoms is right                                       !         !          
   write(15,*) natoms                                                                   !         !              
                                                                                        !         !         
   read(10,*)          ! read blank line                                                !         !           
   write(6,*)          ! write blank line to reproduce input file                       !         !              
   write(15,*)                                                                          !         !               
                                                                                        !         !              
   allocate(coord(natoms,col),mass(natoms))                                             !         !                
                                                                                        !         !             
   call atom_coord(natoms, coord, mass, col)                                            !         !                  
                                                                                        !         !            
   do i=1,natoms                                                                        !         !              
    write(15,*) (coord(i,j),j=1,col), mass(i)                                           !         !            
   enddo                                                                                !         !            
                                                                                        !         !        
                                                                                        !         !             
!---------------------------------------------------------------------------------------!         !              
                                                                                                  !
                                                                                                  !
                                                                                                  !
!---------------------------------------------------------------------------------------!         !
!-computing internuclear distances between each pair of atoms---------------------------!         !                                        
                                                                                        !         !
   allocate(d(natoms,natoms))                                                           !         !
   call compute_distances(natoms,coord, col, d)                                         !         !
                                                                                        !         !
   do i=1,natoms                                                                        !         ! 
     write(20,*) (d(i,j),j=1,natoms)                                                    !  
   enddo                                                                                !
                                                                                        !
!---------------------------------------------------------------------------------------!   



!---------------------------------------------------------------------------------------!   
!-computing LJ potential energy --------------------------------------------------------!   

  open(25,file='pot.inp')

  read(25,*) eps, sigma
                                                                                                       
  LJtot=V(eps,sigma,natoms,d) 
  
  write(6,*)
  write(6,*) "total LJ potential energy=", LJtot 
                                                                                                       
!---------------------------------------------------------------------------------------!   
                                                                                                       
                                                                                                       
                                                                                                       
                                                                                                       
                                                                                                       
!---------------------------------------------------------------------------------------!   
!-computing kinetic energy--------------------------------------------------------------!   

  allocate(vel(natoms,col))

  vel=0.d0           ! initializing velocities to zero
   
  Ttot=T(natoms,mass,vel)

  write(6,*)
  write(6,*) "total kinetic energy=", Ttot

!---------------------------------------------------------------------------------------!   


!---------------------------------------------------------------------------------------!   
!-computing total energy----------------------------------------------------------------!   

  Etot=E(LJtot,Ttot)
  
  write(6,*)
  write(6,*) "total energy = LJtot + Ttot =", LJtot, "+", Ttot, "=", Etot

!---------------------------------------------------------------------------------------!   


!---------------------------------------------------------------------------------------!   
!-computing acceleration array----------------------------------------------------------!   

  allocate(f(natoms,natoms))                      ! array of dV for each d(i,j)
  allocate(acceleration(natoms,col))

  call compute_acc(eps,sigma,natoms,coord, mass, d, acceleration,f)

  do i=1,2
   write(6,*)
  enddo

  write(6,*) "acceleration array"

  do i=1,natoms
     write(6,*) (acceleration(i,k),k=1,col)
  enddo


!---------------------------------------------------------------------------------------!   


   deallocate(f,acceleration)
   deallocate(coord,mass,d,vel)                                                                   !
                                                                                                  !
   stop                                                                                           !           
   end program main                                                                               !     
!-------------------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------!  
!-function that reads number of atoms from input file-----------------------------------!
                                                                                        !
  integer function read_natoms(unit) result(natoms)                                     ! 
          implicit none                                                                 !
          integer, intent (in) :: unit                                                  !
          read(unit,*) natoms                                                           ! 
  return                                                                                !
  end function read_natoms                                                              !
!---------------------------------------------------------------------------------------!

!---------------------------------------------------------------------------------------!
!-function that calulates total LJ potential of the system------------------------------!
                                                                                        !
  double precision function V(eps,sigma,natoms,d)
         implicit none
         double precision, intent(in) :: eps, sigma
         integer, intent(in) :: natoms
         double precision, intent(in) :: d(natoms,natoms)
         double precision :: Vtot, pot
         integer :: i,j
         
         Vtot=0.d0
         do i=1,natoms
           do j=i+1,natoms
              pot=4*eps*((sigma/d(i,j))**12-(sigma/d(i,j))**6)
              Vtot=Vtot+pot
           enddo
         enddo      
         V=Vtot 
  end function V

!---------------------------------------------------------------------------------------!




!---------------------------------------------------------------------------------------!
!-function that calulates the kinetic energy given the masses and velocities------------!

  double precision function T(natoms,mass,vel)
          implicit none
          integer, parameter :: col=3
          integer, intent(in):: natoms
          double precision, intent(in)   :: vel(natoms,col)
          double precision, intent(in) :: mass(natoms)
          integer :: i,j
          double precision :: velocity    ! velocity vector

          T=0.d0
          do i=1,natoms
             velocity=0.d0
               do j=1,col
                 velocity=velocity+(vel(i,j))**2
               enddo
             T=T+mass(i)*velocity
          enddo   
          T=0.5d0*T

  end function T        

!---------------------------------------------------------------------------------------!



!---------------------------------------------------------------------------------------!
!-function that calulates total energy--------------------------------------------------!

  double precision function E(LJtot,Ttot)
          implicit none
          double precision :: LJtot, Ttot

          E=LJtot+Ttot

  end function E

!---------------------------------------------------------------------------------------!




