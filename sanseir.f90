!> \file sanseir.f90 
!> \brief Solves a stochastic SEIR model
!> \author S. Scott Collis
!=============================================================================!

!> \brief Useful constants
module constants
  real, parameter :: zero    = 0.0000000000000000000d+0
  real, parameter :: pt25    = 2.5000000000000000000d-1
  real, parameter :: pt33    = 3.3333333333333333333d-1
  real, parameter :: pt5     = 5.0000000000000000000d-1
  real, parameter :: pt66    = 6.6666666666666666666d-1
  real, parameter :: one     = 1.0000000000000000000d+0
  real, parameter :: onept25 = 1.2500000000000000000d+0
  real, parameter :: onept33 = 1.3333333333333333333d+0
  real, parameter :: onept5  = 1.5000000000000000000d+0
  real, parameter :: two     = 2.0000000000000000000d+0
  real, parameter :: three   = 3.0000000000000000000d+0
  real, parameter :: four    = 4.0000000000000000000d+0
  real, parameter :: pi      = 3.1415926535897932385d+0
  real, parameter :: eps     = 1.0000000000000000000d-9
end module constants

!=============================================================================!

!> \brief SEIR Model parmameters
module seir_model
  implicit none
  character(60) :: title="disease progression model"  !< Title of run
  integer, parameter :: mk=10  !< max number of phases
  integer :: nk=1              !< number of phases
  real :: Ro(mk)               !< Reproduction rate per phase
  real :: tk(mk)               !< Time of each phase
  real :: P                    !< Total population
  real :: Io                   !< Initial infected (in Ic)
  real :: alphat               !< Reciprocal of incupation period: \f$\tilde\alpha = \alpha^{-1}\f$
  real :: gammat               !< Reciprocal of infectious period: \f$\tilde\gamma = \gamma^{-1}\f$
  real :: beta(mk)             !< Contact rate (per phase)
  real :: rho(mk)              !< Testing fraction (per phase)
  real :: c                    !< Mortality rate
  real :: cfr                  !< Case fatality rate
  logical :: use_cfr=.false.   !< use CFR
  real :: alpham               !< Mean value for incubation period, \f$\alpha\f$
  real :: gammam               !< Mean value for infectious period, \f$\gamma\f$
  real :: to=0                 !< Starting time (days)
  real :: tf=250               !< End time (days)
  real :: dt                   !< time step
  real :: Fa=0.2               !< Fraction of infected accounted for
  integer :: nt=0              !< number of time steps
  integer :: Erlang_k=2        !< Shape of Erlang distribution
  integer :: Erlang_n=100      !< Number of Erlang samples
  real    :: alpha_min=1       !< Minimum value of \f$\alpha\f$
  real    :: gamma_min=1       !< Minimum value of \f$\gamma\f$
  integer :: ms=50             !< Number of samples in ensemble

  namelist /model/ title, P, Io, alpham, gammam, Erlang_k, Erlang_n, Ro, &
                   rho, c, Fa, alpha_min, gamma_min, cfr, use_cfr, ms
  namelist /time/ to, tf, nt, tk, nk
end module seir_model 

!=============================================================================!
!> \brief SanSEIR solves a stochastic SEIR disease progression model
program sanseir 
!=============================================================================!
!
! This follows the formulation and approach presented in:
!
! Yao Ye Yeo, Yao-Rui Yeo, Wan-Jin Yeo, "A Computational Model for Estimating 
! the Progression of the COVID-19 Cases in the US West and East Coasts",  
! MedRxiv, March 2020, https://doi.org/10.1101/2020.03.24.20043026
!
! Author:  S. Scott Collis
!
! Copyright: (c)2020 S. Scott Collis, All rights reserved
!
! Written: 3-29-2020 
!=============================================================================!
  use constants
  use seir_model 
  implicit none

  integer :: i, n, ns
  external :: seir
  real, external :: erlang_sample

  integer, parameter :: neq=8
  real, allocatable :: U(:,:), V(:), dUdt(:,:), t(:)

  integer, parameter :: mfile=1
  integer :: narg, iarg, nfile=0, ifile(mfile)
  character(80) :: arg, iname, ofile, sfile, rfile

  logical :: use_rk4=.false.   !< Whether to use RK4 (default is RKCK45).

  real :: alphai               !< Sample of \f$\alpha\f$
  real :: gammai               !< Sample of \f$\gamma\f$
  real :: alphas               !< Average \f$\alpha\f$
  real :: gammas               !< Average \f$\famma\f$
  real :: Cr                   !< Cases (reported)
  real :: dCrdt                !< rate of cases reported
!=============================================================================!

! parse arguments

  narg = iargc()
  do iarg = 1, narg
    call getarg(iarg,arg)
    if (arg(1:1) .ne. '-') then
      nfile = nfile + 1
      if (nfile .gt. mfile) then
        write(*,*) '>> Error in argument list, too many file names'
        call exit(1)
      end if
      ifile(nfile) = iarg
    else
      select case (arg(1:3))
      case ('-rk')
        write(*,'("Using RK4 instead of RKCK45")')
        use_rk4 = .true.
      case ('-h')
        write(*,"('--------------------------------------------------')")
        write(*,"('Usage: sanseir [options] [file]')")
        write(*,"('--------------------------------------------------')")
        write(*,"('   -h:  this help                                 ')")
        write(*,"('  -rk:  use RK4 instead of RKCK45                 ')")
        write(*,"('--------------------------------------------------')")
        call exit(0)
      case default
        write(*,"('Argument ',i2,' ignored.')") iarg
      end select
    end if
  end do

! Initialize parameters

  Ro = zero
  beta = zero
  rho = zero

! Read in model parameters

  if (nfile.gt.0) then
    call getarg(ifile(1),iname)
  else
    iname = "seir.inp"
  endif
  open(10,file=iname,status='old')
  read(10,nml=time)
  read(10,nml=model) 
  close(10)

  write(*,'(80("="))')
  write(*,'("SanSEIR: " a)') title
  write(*,'(80("="))')

  write(*,'("Echo of namelist input:")')
  write(*,nml=time)
  if(nk.gt.mk) then
    write(*,*) "nk > mk, terminating"
    call exit(1)
  endif 
  write(*,nml=model)
  
! allocate storage

  allocate ( U(neq,0:nt), V(neq), dUdt(neq,0:nt), t(0:nt) )

! Setup time

  dt = (tf-to)/nt 
  write(*,"('Using timestep, dt = ',1pe12.4e3)") dt

! Setup model parameters

  call random_init(.true.,.true.)
! call random_seed()
  alphai = max(one,erlang_sample(Erlang_k,alpham))
  gammai = max(one,erlang_sample(Erlang_k,gammam))

  alphat = one/alphai
  gammat = one/gammai
  beta   = Ro*gammat

! set the initial condition
! U = {S, E, Ih, Ic, Rh, Rc, Dh, Dc}

#ifdef SANSEIR_ORIGINAL_IC 
  U(:,:) = zero
  U(1,0) = P-Io
  U(2,0) = Io
#else
  U(:,:) = zero
  U(1,0) = P-Io
  U(4,0) = Io
#endif

#ifdef TEST_ERLANG
  alphas = 0
  gammas = 0
  do n = 1, 100000 
    alphai = max(alpha_min,Erlang_sample(Erlang_k,alpham))
    gammai = max(gamma_min,erlang_sample(Erlang_k,gammam))
    alphas = alphas + alphai
    gammas = gammas   + gammai
  end do
  write(*,*) alphas/100000, gammas/100000
  call exit(0)
#endif

! begin the time loop

  do ns = 1, ms
    write(*,'("Computing sample: ",i0)') ns
    t(0) = to
    do i = 1, nt
      U(:,i) = zero
      do n = 1, Erlang_n 
        V(:)   = zero
        alphai = max(alpha_min,Erlang_sample(Erlang_k,alpham))
        gammai = max(gamma_min,Erlang_sample(Erlang_k,gammam))
        alphat = one/alphai
        gammat = one/gammai
        beta   = Ro*gammat
        if (use_rk4) then
          call rk4(neq,U(:,i-1),V,t(i-1),dt,seir)
        else
          call rkck45(neq,U(:,i-1),V,t(i-1),dt,seir)
        endif
        U(:,i) = U(:,i) + V(:)
        call seir(neq, U(:,i-1), t(i-1), V)
        dUdt(:,i) = dUdt(:,i) + V(:)
      end do
      U(:,i) = U(:,i)/Erlang_n
      dUdt(:,i) = dUdt(:,i)/Erlang_n
      t(i) = t(i-1) + dt
    end do
    write(ofile,"('output.',i0)") ns
    write(sfile,"('scaled.',i0)") ns
    write(rfile,"('rates.',i0)") ns
    open(unit=10,file=ofile)
    open(unit=20,file=sfile)
    open(unit=30,file=rfile)
!=============================================================================!
! Column Index:   1  2  3   4   5   6   7   8   9     10    11  12 13
!=============================================================================!
    write(10,'("# t, S, E, Ih, Ic, Rh, Rc, Dh, Dc, Ih+0.2Ic, R, D, Cr")')
    write(20,'("# t, S, E, Ih, Ic, Rh, Rc, Dh, Dc, Ih+0.2Ic, R, D, Cr")')
    write(30,'("# t, S, E, Ih, Ic, Rh, Rc, Dh, Dc, Ih+0.2Ic, R, D, Cr")')
    do i = 0, nt
      Cr    = U(5,i)+Fa*U(6,i)+U(7,i)+U(8,i)              ! Cases (reported)
      dCrdt = dUdt(5,i)+Fa*dUdt(6,i)+dUdt(7,i)+dUdt(8,i)  ! Cases (reported)
      write(10,10) t(i), U(:,i), U(3,i)+Fa*U(4,i), U(5,i)+U(6,i), &
                   U(7,i)+U(8,i), Cr 
      write(20,10) t(i), U(:,i)/P, (U(3,i)+Fa*U(4,i))/P, (U(5,i)+U(6,i))/P, &
                   (U(7,i)+U(8,i))/P, Cr/P
      write(30,10) t(i), dUdt(:,i), dUdt(3,i)+Fa*dUdt(4,i), &
                   dUdt(5,i)+dUdt(6,i), dUdt(7,i)+dUdt(8,i), dCrdt
    end do
    close(10)
    close(20)
    close(30)
  end do

10 format(12(1pe16.8E3,', '),1pe16.8e3)
  stop
end program sanseir

!=============================================================================!
!> \brief Erlang probability distribution 
!> \param[in] x value for which probability is sought
!> \param[in] k Erlang shape parameter 
!> \param[in] mean Mean value of the distribution 
!> \return probability from Erlang_k distribution 
function Erlang(x, k, mean)
!=============================================================================!
  implicit none
  integer k
  real Erlang, x, mean, lambda, invfact
!=============================================================================!
  lambda = k/mean
  invfact = 1.0/Gamma(real(k))
  Erlang = lambda**k*x**(k-1)*exp(-lambda*x)*invfact
end function Erlang

!=============================================================================!
!> \brief Sample from an Erlang distribution 
!> \param[in] k Erlang shape parameter 
!> \param[in] mean Mean value of the distribution 
!> \return sample from Erlang_k distribution 
function Erlang_sample(k, mean)
!=============================================================================!
  implicit none
  integer k, i
  real Erlang_sample, r(k), prod, mean, lambda
!=============================================================================!
  if (k.eq.0) then
    erlang_sample = mean 
    return
  end if
  lambda = k/mean
  call random_number(r)
  prod = r(1)
  do i = 2, k
    prod = prod*r(i)
  end do
  Erlang_sample = -1.0/lambda*log(prod) 
end function Erlang_sample

!=============================================================================!
!> \brief Compute the time-derivative of the SEIR model
!> \param[in] neq number of equations
!> \param[in] U state fector at time t
!> \param[in] t current time
!> \param[out] dUdt time-derivative at time t 
subroutine seir(neq, U, t, dUdt)
!=============================================================================!
  use constants
  use seir_model 
  implicit none
  real, external :: erlang_sample
  integer neq, i, k
  real U(neq), t, dUdt(neq), b, r
  real N, Ni, S, E, Ih, Ic, Rh, Rc, It, ct
!=============================================================================!

! compute the effective total population

  N = zero
  do i = 1, 6 
    N = N + U(i)
  end do 
  Ni = one/N  ! inverse of total effective population

  S  = U(1)   ! susceptable population
  E  = U(2)   ! exposed population
  Ih = U(3)   ! Infected population in a hospital
  Ic = U(4)   ! Infected population in community
  Rh = U(5)   ! Recovered population in hospital
  Rc = U(6)   ! Recovered population in community

! use time varying contact rate and testing fractions

  b = beta(1) 
  do k = 1, nk
    if (t.gt.tk(k)) then
      b = beta(k)         ! current contact rate
      r = rho(k)          ! current testing fraction
    endif
  end do

  if (.not.use_cfr) then 
    dUdt(1) = -b*Ic*S*Ni
    dUdt(2) =  b*Ic*S*Ni - alphat*E
    dUdt(3) =  alphat*r*E - gammat*Ih 
    dUdt(4) =  alphat*(one-r)*E - gammat*Ic 
    dUdt(5) =  gammat*(one-c)*Ih 
    dUdt(6) =  gammat*(one-c)*Ic
    dUdt(7) =  gammat*c*Ih
    dUdt(8) =  gammat*c*Ic
  else
    It = (Ih + Fa*Ic)        ! confirmed cases
    ct = cfr * It/(Ih+Ic)    ! convert CFR to c
    dUdt(1) = -b*Ic*S*Ni
    dUdt(2) =  b*Ic*S*Ni - alphat*E
    dUdt(3) =  alphat*r*E - gammat*Ih 
    dUdt(4) =  alphat*(one-r)*E - gammat*Ic 
    dUdt(5) =  gammat*(one-ct)*Ih 
    dUdt(6) =  gammat*(one-ct)*Ic
    dUdt(7) =  gammat*ct*Ih
    dUdt(8) =  gammat*ct*Ic
  endif

  return
end subroutine seir

!=============================================================================!
!> \brief Advance one time step using Runge-Kutta Cash-Karp 4/5 method
!> \param[in] neq number of equations
!> \param[in] yo initial value
!> \param[out] yf final value
!> \param[in] to intial time
!> \param[in] h time step
!> \param[in] FUNC function to integrate
subroutine RKCK45(neq, yo, yf, to, h, FUNC)
!=============================================================================!
!
! Advance one time step Runge-Kutta Cash-Karp method
!
!=============================================================================!
  implicit none
  external FUNC
  integer  neq, i, j, k, m, n
  real     to, h, t
  real     yo(neq), yf(neq), yt(neq)
  real     yk(neq,6), ye(neq)
!=============================================================================!
  real b(6,5)
  real a(6), c(6), d(6)
  data a / 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 /
  data b / 0.0, 0.2, 0.075, 0.3, -0.2037037037037037, &
 &         0.029495804398148147, &
 &         0.0, 0.0, 0.225, -0.9, 2.5, 0.341796875, &
 &         0.0, 0.0, 0.0, 1.2, -2.5925925925925926, &
 &         0.041594328703703706, &
 &         0.0, 0.0, 0.0, 0.0, 1.2962962962962963, &
 &         0.40034541377314814, &
 &         0.0, 0.0, 0.0, 0.0, 0.0, 0.061767578125 /
  data c / 0.09788359788359788, 0.0, 0.4025764895330113, &
 &         0.21043771043771045, 0.0, 0.2891022021456804 /
  data d / -0.004293774801587311, 0.0, 0.018668586093857853, &
 &         -0.034155026830808066, -0.019321986607142856, &
 &         0.03910220214568039 /

! Test data

#ifdef DEBUG
  do i = 1, 6
    do j = 1, 5
      write(*,*) i, j, b(i,j)
    end do
  end do
  stop
#endif

! Stage 1 - 6

  do m = 1, 6
    t = to + a(m)*h
    do n = 1, neq
      yt(n) = yo(n)
    end do
    do k = 1, m-1
      do n = 1, neq
        yt(n) = yt(n) + b(m,k)*yk(n,k)
      end do
    end do
    call FUNC(neq, yt, t, yk(1,m))
    do n = 1, neq
      yk(n,m) = h * yk(n,m)
    end do
  end do

! Final solution and error

  do n = 1, neq
    yf(n) = yo(n)
    ye(n) = 0.0
  end do
  do k = 1, 6
    do n = 1, neq
      yf(n) = yf(n) + c(k)*yk(n,k)
      ye(n) = ye(n) + d(k)*yk(n,k)
    end do
  end do

  return
end subroutine RKCK45

!=============================================================================!
!> \brief Advance one time step using fourth order (real) Runge-Kutta
!> \param[in] neq number of equations
!> \param[in] yo initial value
!> \param[out] yf final value
!> \param[in] to intial time
!> \param[in] h time step
!> \param[in] FUNC function to integrate
subroutine RK4(neq, yo, yf, to, h, FUNC)
!=============================================================================!
!
! Advance one time step using fourth order (real) Runge-Kutta
!
!=============================================================================!
  external FUNC
  integer  neq
  real     to, h
  real     yo(neq), yf(neq)
  real     f(neq), k1(neq), k2(neq), k3(neq), k4(neq), q(neq)

  call FUNC(neq, yo, to, f)
  do j = 1 , neq
    k1(j) = h*f(j)
    q(j) = yo(j) + 0.5*k1(j)
  end do
  call FUNC(neq, q, to+0.5*h, f)
  do j = 1 , neq
    k2(j) = h*f(j)
    q(j) = yo(j) + 0.5*k2(j)
  end do
  call FUNC(neq, q, to+0.5*h, f)
  do j = 1 , neq
    k3(j) = h*f(j)
    q(j) = yo(j) + k3(j)
  end do
  call FUNC(neq, q, to+h, f)
  do j = 1 , neq
    k4(j) = h*f(j)
    yf(j) = yo(j)+k1(j)/6.+(k2(j)+k3(j))/3.+k4(j)/6.
  end do

  return
end subroutine RK4
