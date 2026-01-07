!---------------------------------------------------------------------------------------!
!-Module to create explicit interface of allocatable dummy arguments--------------------!
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
!---------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------!                         
!-Molecular Dynamics Code---------------------------------------------------------------!         !                  
!---------------------------------------------------------------------------------------!         !                            
  program main                                                                                    !                              
                                                                                                  !                                         
                                                                                                  !                           
                                                                                                  !                                        
!---------------------------------------------------------------------------------------!         !                          
!-declaration block---------------------------------------------------------------------!         !                       
                                                                                        !         !                   
  use reading                                                                           !         !               
  implicit none                                                                         !         !                       
  integer read_natoms                                                                   !         !               
  integer unit, natoms, status                                                          !         !                  
  double precision, allocatable :: mass(:)                                              !         !                    
  double precision, allocatable :: coord(:,:)                                           !         !                
  integer, parameter :: col=3                                                           !         !               
  integer i, j                                                                          !         !              
!---------------------------------------------------------------------------------------!         !                
!---------------------------------------------------------------------------------------!         !                       
                                                                                                  !              
!---------------------------------------------------------------------------------------!         !                    
!-opening of input and output file and reading of input variables-----------------------!         !              
                                                                                        !         !            
  open(unit=10, file='dat.inp')                                                         !         !                 
  open(15, file='coord.out')                                                            !         !              
                                                                                        !         !               
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
   deallocate(coord,mass)                                                               !         !                  
!---------------------------------------------------------------------------------------!         !             
!---------------------------------------------------------------------------------------!         !              
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

