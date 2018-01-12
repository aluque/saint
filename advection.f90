module cmp
contains
  
  subroutine deriv_u_upwind(dx, c, u, du)
    double precision, intent(in) :: c, dx, u(:)
    double precision, intent(out) :: du(:)
    
    double precision :: f(size(du) + 1)
  
    f(:) = c * u(1: size(u) - 1)
    du(:) = -(f(2: size(f)) - f(1: size(f) - 1)) / dx
  end subroutine deriv_u_upwind
  
  subroutine forward_euler_serial(dx, dt, c, u0, u, n, m)
    ! Note the inverted order of dimensions
    double precision, intent(in) :: dx, dt, c, u0(n)
    double precision, intent(inout) :: u(n + 2, m)
    integer, intent(in) :: n, m

    !f2py integer intent(hide), depend(u) :: n=size(u,1)-2, m=size(u,2)

    double precision :: du(n)
  
    u(2:size(u, 1) - 1, 1) = u0
    do j = 2, m
       ! Periodic b.c.
       u(1, j - 1) = u(size(u, 1) - 1, j - 1)
       u(size(u, 1), j - 1) = u(1, j - 1)

       call deriv_u_upwind(dx, c, u(:, j - 1), du)
       do i=2, size(u, 1) - 1
          u(i, j) = u(i, j - 1) + dt * du(j - 1)
       end do
    end do
  end subroutine forward_euler_serial

  subroutine forward_euler_omp(dx, dt, c, u0, u, n, m)
    ! Note the inverted order of dimensions
    double precision, intent(in) :: dx, dt, c, u0(n)
    double precision, intent(inout) :: u(n + 2, m)
    integer, intent(in) :: n, m

    !f2py integer intent(hide), depend(u) :: n=size(u,1)-2, m=size(u,2)

    double precision :: du(n)
  
    u(2:size(u, 1) - 1, 1) = u0
    do j = 2, m
       ! Periodic b.c.
       u(1, j - 1) = u(size(u, 1) - 1, j - 1)
       u(size(u, 1), j - 1) = u(1, j - 1)

       call deriv_u_upwind(dx, c, u(:, j - 1), du)
       !$OMP PARALLEL DO 
       do i=2, size(u, 1) - 1
          u(i, j) = u(i, j - 1) + dt * du(j - 1)
       end do
       !$OMP END PARALLEL DO
    end do
  end subroutine forward_euler_omp

  subroutine forward_euler_vector(dx, dt, c, u0, u, n, m)
    ! Note the inverted order of dimensions
    double precision, intent(in) :: dx, dt, c, u0(n)
    double precision, intent(inout) :: u(n + 2, m)
    integer, intent(in) :: n, m

    !f2py integer intent(hide), depend(u) :: n=size(u,1)-2, m=size(u,2)

    double precision :: du(n)
  
    u(2:size(u, 1) - 1, 1) = u0
    do j = 2, m
       ! Periodic b.c.
       u(1, j - 1) = u(size(u, 1) - 1, j - 1)
       u(size(u, 1), j - 1) = u(1, j - 1)

       call deriv_u_upwind(dx, c, u(:, j - 1), du)
       u(2:size(u, 1) - 1, j) = u(2:size(u, 1) - 1, j - 1) + dt * du
     
    end do
  end subroutine forward_euler_vector
end module cmp