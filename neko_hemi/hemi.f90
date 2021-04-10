module user
  use neko
  implicit none

  type(field_t) :: up
  
contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_usr_ic => user_ic
    u%fluid_usr_if => user_inflow_eval
!    u%usr_msh_setup => user_mesh_scale
!    u%usr_chk => usr_calc_quantities
  end subroutine user_setup


  subroutine set_prof
    integer icalld, i, ntot
    save    icalld
    data    icalld  /0/
    real(kind=dp) :: u0, rad, visc, delta, hk, uk, rk, re
    if (icalld.ne.0) return
    icalld = icalld + 1
    !    Compute the Blasius profile
    u0   = 1d0
    visc = 1d0/1400d0 !param(2)  !!!<M<<<<<<

    rad   = 0.5d0
    delta = 1.2d0*rad  ! Thorsten

    re = 1d0/visc
    hk = 0.55d0
    uk = blasius_sin(hk,delta,u0)
    rk = uk*hk/visc
    do i=1, up%dof%n_dofs
       up%x(i,1,1,1)=blasius_sin(up%dof%z(i,1,1,1),delta,u0)
    enddo
    return
  end subroutine set_prof

  function blasius(y,delta,u)
    !    Return the velocity at a given y value for specified delta and U.
    integer icalld, i
    save    icalld
    data    icalld /0/
    real(kind=dp) :: blasius, y, delta, u, eta, vel
    real(kind=dp), dimension(45) :: u0, y0, work
    save      u0    ,y0    ,work
    data      u0 / 0.00000d0 , 0.06641d0 , 0.13277d0 , 0.19894d0 , 0.26471d0 &
                 , 0.32979d0 , 0.39378d0 , 0.45627d0 , 0.51676d0 , 0.57477d0 &
                 , 0.62977d0 , 0.68132d0 , 0.72899d0 , 0.77246d0 , 0.81152d0 &
                 , 0.84605d0 , 0.87609d0 , 0.90177d0 , 0.92333d0 , 0.94112d0 &
                 , 0.95552d0 , 0.96696d0 , 0.97587d0 , 0.98269d0 , 0.98779d0 &
                 , 0.99155d0 , 0.99425d0 , 0.99616d0 , 0.99748d0 , 0.99838d0 &
                 , 0.99898d0 , 0.99937d0 , 0.99961d0 , 0.99977d0 , 0.99987d0 &
                 , 0.99992d0 , 0.99996d0 , 0.99998d0 , 0.99999d0 , 1.00000d0 &
                 , 1.00000d0 , 1.00000d0 , 1.00000d0 , 1.00000d0 , 1.00000d0 /
    if (icalld.eq.0) then
       !        Initialize Blasius profile and spline fitting routine.
       icalld=1
       do i=1,45
          y0(i)=dble(i-1)/5.0d0
       end do
       call spline(y0,u0,45,work)
    endif
    eta=5.0d0*y/delta
    if (eta.gt.8.5d0) then
       blasius=u
    else
       call splint(y0,u0,work,45,eta,vel)
       blasius=vel*u
    endif
  end function blasius

  function blasius_sin(y,delta,U)

    !     Return the velocity at a given y value for specified delta and U.
    real(kind=dp) :: blasius_sin, y, delta, U
    real(kind=dp) :: pi,one, arg
    save pi,one
    data pi,one /0d0,1d0/

    if (pi.eq.0) pi=4d0*atan(one)

    arg = .5d0*pi*y/delta
    blasius_sin = sin(arg)
    
    if (arg.gt.0.5d0*pi) blasius_sin = 1d0
    
    blasius_sin = U*blasius_sin
    
    end function blasius_sin

    subroutine spline(x,y,n,y2)
      integer, parameter :: nmax=100
      integer :: n, i, ir, il, k
      real(kind=dp) , dimension(n) :: x,y,y2
      real(kind=dp) :: u(nmax), sig, p, qn,un
      y2(1)=0.0d0
      u(1) =0.0d0
      do  i=2,n-1
         ir=i+1
         il=i-1
         sig=(x(i)-x(il))/(x(ir)-x(il))
         p=sig*y2(il)+2d0
         y2(i)=(sig-1.)/p
         u(i)= ( 6d0* &
              ( (y(ir)-y(i))/(x(ir)-x(i))-(y(i)-y(il))/ (x(i)-x(il) ) ) &
              / (x(ir)-x(il)) &
              - sig*u(il) )/p
      end do
      qn=0d0
      un=0d0
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
      end do
    end subroutine spline

    subroutine splint(xa,ya,y2a,n,x,y)
!     p. 88-89, numerical recipes
      integer :: n, klo, khi, k
      real(kind=dp) xa(n),ya(n),y2a(n), x, y, h, a, b
      klo=1
      khi=n
1     if ((khi-klo).gt.1) then
         k=(khi+klo)/2
         if (xa(k).gt.x) then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+ &
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
    end subroutine splint
  
  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    integer :: i
    real(kind=dp) :: uvw(3)
    type(file_t) :: ff

    ff = file_t('up.fld')
    call field_init(up, u%dof, 'up')

    call set_prof

    call ff%write(up)

    do i = 1, u%dof%size()
       uvw = 0d0 !tgv_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = up%x(i, 1, 1, 1)
       v%x(i,1,1,1) = 0d0
       w%x(i,1,1,1) = 0d0
    end do
    p = 0d0
  end subroutine user_ic

  subroutine user_inflow_eval(u, v, w, x, y, z, nx, ny, nz, ix, iy, iz, ie)
    real(kind=dp), intent(inout) :: u
    real(kind=dp), intent(inout) :: v
    real(kind=dp), intent(inout) :: w
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: y
    real(kind=dp), intent(in) :: z
    real(kind=dp), intent(in) :: nx
    real(kind=dp), intent(in) :: ny
    real(kind=dp), intent(in) :: nz
    integer, intent(in) :: ix
    integer, intent(in) :: iy
    integer, intent(in) :: iz
    integer, intent(in) :: ie

    u = up%x(ix, iy, iz, ie)
    v = 0d0
    w = 0d0
    
  end subroutine user_inflow_eval
  


end module user
