    Module qconserv
    Use Precision
    Use Constants
    Implicit None

    Integer, allocatable        :: g_map(:,:,:)
    Integer, allocatable        :: igvecs(:,:)
    Integer                     :: max_g
    Real (kind=pr), allocatable :: gvecs(:,:)

    Contains

    Subroutine init_g_map(rgvecs, nao)
    Implicit None
    Integer,            Intent(in)  :: nao
    Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
    Real (kind=pr)                  :: max_val
    Integer                         :: idx, nx, ny, nz
    
    ! Allocate both real and integer arrays
    allocate(gvecs(nao,3))
    allocate(igvecs(nao,3))
    gvecs = rgvecs

    max_val = maxval(abs(rgvecs))
    max_g = nint(max_val) * 5 + 1
    allocate(g_map(-max_g:max_g, -max_g:max_g, -max_g:max_g))
    g_map = 0

    Do idx = 1, nao
        nx = nint(rgvecs(idx, 1))
        ny = nint(rgvecs(idx, 2))
        nz = nint(rgvecs(idx, 3))
        
        igvecs(idx, 1) = nx
        igvecs(idx, 2) = ny
        igvecs(idx, 3) = nz
        
        g_map(nx, ny, nz) = idx
    End Do
    End Subroutine


    Pure Subroutine qconserv2(i, a, j, b)
    Implicit None
    Integer,            Intent(in)  :: i, j, a
    Integer,            Intent(out) :: b
    Integer                         :: bx, by, bz
    
    bx = igvecs(i,1) + igvecs(j,1) - igvecs(a,1)
    by = igvecs(i,2) + igvecs(j,2) - igvecs(a,2)
    bz = igvecs(i,3) + igvecs(j,3) - igvecs(a,3)

    if (abs(bx) <= max_g .and. abs(by) <= max_g .and. abs(bz) <= max_g) then
        b = g_map(bx, by, bz)
    else
        b = 0
    end if
    End Subroutine


    Pure Subroutine qconserv_c(i, j, k, a, b, c)
    Implicit None
    Integer,            Intent(in)  :: i, j, k, a, b
    Integer,            Intent(out) :: c
    Integer                         :: cx, cy, cz
    
    cx = igvecs(i,1) + igvecs(j,1) + igvecs(k,1) - igvecs(a,1) - igvecs(b,1)
    cy = igvecs(i,2) + igvecs(j,2) + igvecs(k,2) - igvecs(a,2) - igvecs(b,2)
    cz = igvecs(i,3) + igvecs(j,3) + igvecs(k,3) - igvecs(a,3) - igvecs(b,3)

    if (abs(cx) <= max_g .and. abs(cy) <= max_g .and. abs(cz) <= max_g) then
        c = g_map(cx, cy, cz)
    else
        c = 0
    end if
    End Subroutine

    End Module
