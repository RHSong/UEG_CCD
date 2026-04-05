    Module qconserv
    Use Precision
    Use Constants
    Implicit None

    Integer, allocatable :: g_map(:,:,:)
    Integer :: max_g
    Real (kind=pr), allocatable :: gvecs(:,:)

    Contains

    Subroutine init_g_map(rgvecs, nao)
    Implicit None
    Integer,            Intent(in)  :: nao
    Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
    Real (kind=pr)                  :: max_val
    Integer                         :: idx, nx, ny, nz
    allocate(gvecs(nao,3))
    gvecs = rgvecs

    max_val = maxval(abs(rgvecs))
    max_g = nint(max_val) * 5 + 1
    allocate(g_map(-max_g:max_g, -max_g:max_g, -max_g:max_g))
    g_map = 0
                                                                                   
    Do idx = 1, nao
        nx = nint(rgvecs(idx, 1))
        ny = nint(rgvecs(idx, 2))
        nz = nint(rgvecs(idx, 3))
        g_map(nx, ny, nz) = idx
    End Do
    End Subroutine
                                                                                   
    Subroutine qconserv2(i, a, j, b)
    Implicit None
    Integer,            Intent(in)  :: i, j, a
    Integer,            Intent(out) :: b
    integer :: bx, by, bz
    bx = nint(gvecs(i,1) + gvecs(j,1) - gvecs(a,1))
    by = nint(gvecs(i,2) + gvecs(j,2) - gvecs(a,2))
    bz = nint(gvecs(i,3) + gvecs(j,3) - gvecs(a,3))
                                                                                   
    if (abs(bx) <= max_g .and. abs(by) <= max_g .and. abs(bz) <= max_g) then
        b = g_map(bx, by, bz)
    else
        b = 0
    end if
    End Subroutine
                                                                                   
    Subroutine qconserv_c(i, j, k, a, b, c)
    Implicit None
    Integer,            Intent(in)  :: i, j, k, a, b
    Integer,            Intent(out) :: c
    integer :: cx, cy, cz
    cx = nint(gvecs(i,1) + gvecs(j,1) + gvecs(k,1) - gvecs(a,1) - gvecs(b,1))
    cy = nint(gvecs(i,2) + gvecs(j,2) + gvecs(k,2) - gvecs(a,2) - gvecs(b,2))
    cz = nint(gvecs(i,3) + gvecs(j,3) + gvecs(k,3) - gvecs(a,3) - gvecs(b,3))
                                                                                   
    if (abs(cx) <= max_g .and. abs(cy) <= max_g .and. abs(cz) <= max_g) then
        c = g_map(cx, cy, cz)
    else
        c = 0
    end if
    End Subroutine

    End Module
