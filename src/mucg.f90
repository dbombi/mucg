! Program for prediction of TTT digrams
! Based on MAP programs MUCG46B, MUCG83
    program main_mucg
    use parsubfuns
    implicit none
! Variables
    integer :: step
    real(8), dimension(12) :: c
    real(8), dimension(sz) :: stctemp, stxeq, stxeq2
    real(8) :: mechd, mechsta, gs1
    real(8) :: ws, bs, ms, stored
    logical :: prn !Print steps (.true. or .false.)
! Initialization of variables
    c=0d0
    mechd=0d0
    mechsta=0d0
    gs1=0d0
    stored=0d0 !Stored energy
    step = 20 ! Main loop temperature step
    prn = .true.
! Read input parameters
    call composition(c) ! Reads composition 
    call mech(mechd,mechsta,gs1) ! Reads mechanical properties
!
    call mucg(c,ws,bs,ms,stored,mechd,mechsta,gs1,step,prn,stctemp,stxeq,stxeq2) !Calculate TTT diagram
!
   end program main_mucg