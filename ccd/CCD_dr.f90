    Module UEG_CCD_dr
    Use Precision
    Use Constants
    Implicit None

    Contains

    Subroutine CCD_Res(f, u, t, r2, no, nao)
    Implicit None
    Integer,            Intent(in)  :: no, nao
    Real (kind=pr),     Intent(in)  :: f(nao,nao), u(nao,nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: t(no+1:nao,no,no+1:nao,no)
    Real (kind=pr),     Intent(out) :: r2(no+1:nao,no,no+1:nao,no)
    Integer                         :: a,b,c,d,i,j,k,l,nv
    nv = nao - no
    !$omp parallel default(shared)

    !$omp single
    r2 = 0.0
    !$omp end single

    ! ladder
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                t(c, i, d, j) * u(a, d, b, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ladder
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                2 * t(c, i, d, j) * u(a, c, b, d)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ladder
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do l=1, no
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                t(a, k, b, l) * u(l, i, k, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ladder
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do l=1, no
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                2 * t(a, k, b, l) * u(k, i, l, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ladder
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        t(a, k, b, l) * t(c, i, d, j) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ladder
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, k, b, l) * t(c, i, d, j) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(a, i, c, k) * u(b, c, k, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(a, k, c, i) * u(b, j, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(b, j, c, k) * u(a, c, k, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(b, k, c, j) * u(a, i, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                t(b, k, c, j) * u(a, c, k, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                4 * t(a, i, c, k) * u(b, j, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                4 * t(b, j, c, k) * u(a, i, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                t(a, k, c, i) * u(b, c, k, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        4 * t(a, i, c, k) * t(b, j, d, l) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        4 * t(a, i, c, k) * t(b, l, d, j) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        4 * t(a, k, c, i) * t(b, j, d, l) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        t(a, k, c, i) * t(b, l, d, j) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, i, c, k) * t(b, l, d, j) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, k, c, i) * t(b, j, d, l) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, k, c, i) * t(b, l, d, j) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        8 * t(a, i, c, k) * t(b, j, d, l) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(a, j, c, k) * u(b, i, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(a, k, c, j) * u(b, c, k, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(b, i, c, k) * u(a, j, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                2 * t(b, k, c, i) * u(a, c, k, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                t(a, j, c, k) * u(b, c, k, i)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                t(a, k, c, j) * u(b, i, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                t(b, i, c, k) * u(a, c, k, j)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        do k=1, no
                            r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                t(b, k, c, i) * u(a, j, k, c)&
                            )
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        4 * t(a, j, c, k) * t(b, i, d, l) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        t(a, j, c, k) * t(b, l, d, i) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        t(a, k, c, j) * t(b, i, d, l) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                                        t(a, k, c, j) * t(b, l, d, i) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, j, c, k) * t(b, i, d, l) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, j, c, k) * t(b, l, d, i) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, k, c, j) * t(b, i, d, l) * u(k, c, l, d)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    ! x ring
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do d=no + 1, no + nv
                        do c=no + 1, no + nv
                            do l=1, no
                                do k=1, no
                                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                                        2 * t(a, k, c, j) * t(b, l, d, i) * u(k, d, l, c)&
                                    )
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    4 * t(a, i, b, k) * t(c, l, d, j) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    4 * t(a, i, c, j) * t(b, k, d, l) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    4 * t(a, k, b, j) * t(c, i, d, l) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    4 * t(a, k, c, l) * t(b, j, d, i) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    t(a, j, b, k) * t(c, i, d, l) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    t(a, j, c, i) * t(b, k, d, l) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    t(a, k, b, i) * t(c, l, d, j) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) - ( &
    !                                    t(a, k, c, l) * t(b, i, d, j) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, i, b, k) * t(c, l, d, j) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, i, c, j) * t(b, k, d, l) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, j, b, k) * t(c, i, d, l) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, j, c, i) * t(b, k, d, l) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, k, b, i) * t(c, l, d, j) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, k, b, j) * t(c, i, d, l) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, k, c, l) * t(b, i, d, j) * u(k, d, l, c)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do
    !! mosaic
    !!$omp do schedule(static)
    !do j=1, no
    !    do b=no + 1, no + nv
    !        do i=1, no
    !            do a=no + 1, no + nv
    !                do d=no + 1, no + nv
    !                    do c=no + 1, no + nv
    !                        do l=1, no
    !                            do k=1, no
    !                                r2(a, i, b, j) = r2(a, i, b, j) + ( &
    !                                    2 * t(a, k, c, l) * t(b, j, d, i) * u(k, c, l, d)&
    !                                )
    !                            end do
    !                        end do
    !                    end do
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
    !!$omp end do

    ! driving
    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
    
                    r2(a, i, b, j) = r2(a, i, b, j) - ( &
                        u(a, j, b, i)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
    
                    r2(a, i, b, j) = r2(a, i, b, j) + ( &
                        2 * u(a, i, b, j)&
                    )
    
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        r2(a, i, b, j) = r2(a, i, b, j) - ( &
                            f(a, c) * t(b, i, c, j)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        r2(a, i, b, j) = r2(a, i, b, j) - ( &
                            f(b, c) * t(a, j, c, i)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        r2(a, i, b, j) = r2(a, i, b, j) + ( &
                            2 * f(a, c) * t(b, j, c, i)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do c=no + 1, no + nv
                        r2(a, i, b, j) = r2(a, i, b, j) + ( &
                            2 * f(b, c) * t(a, i, c, j)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do k=1, no
                        r2(a, i, b, j) = r2(a, i, b, j) - ( &
                            2 * f(k, i) * t(a, k, b, j)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do k=1, no
                        r2(a, i, b, j) = r2(a, i, b, j) - ( &
                            2 * f(k, j) * t(a, i, b, k)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do k=1, no
                        r2(a, i, b, j) = r2(a, i, b, j) + ( &
                            f(k, i) * t(a, j, b, k)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp do schedule(static)
    do j=1, no
        do b=no + 1, no + nv
            do i=1, no
                do a=no + 1, no + nv
                    do k=1, no
                        r2(a, i, b, j) = r2(a, i, b, j) + ( &
                            f(k, j) * t(a, k, b, i)&
                        )
                    end do
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp end parallel
    End Subroutine

    End Module
