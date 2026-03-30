    Module UEG_CCD_T
    Use Precision
    Use Constants
    Implicit None

    Integer, allocatable :: g_map(:,:,:)
    Integer :: max_g

    Contains

    Subroutine init_g_map(rgvecs, nao)
    Implicit None
    Integer,            Intent(in)  :: nao
    Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
    Real (kind=pr)                  :: max_val
    Integer                         :: idx, nx, ny, nz
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

    Subroutine qconserv_c(rgvecs, i, j, k, a, b, c, nao)
    Implicit None
    Integer,            Intent(in)  :: nao, i, j, k, a, b
    Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
    Integer,            Intent(out) :: c
    integer :: cx, cy, cz
    cx = nint(rgvecs(i,1) + rgvecs(j,1) + rgvecs(k,1) - rgvecs(a,1) - rgvecs(b,1))
    cy = nint(rgvecs(i,2) + rgvecs(j,2) + rgvecs(k,2) - rgvecs(a,2) - rgvecs(b,2))
    cz = nint(rgvecs(i,3) + rgvecs(j,3) + rgvecs(k,3) - rgvecs(a,3) - rgvecs(b,3))

    if (abs(cx) <= max_g .and. abs(cy) <= max_g .and. abs(cz) <= max_g) then
        c = g_map(cx, cy, cz)
    else
        c = 0
    end if
    End Subroutine

    Subroutine build_imd(qconserv, eri, t2, imd, a, i, b, j, c, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: imd
    Integer                         :: f, m
    imd = Zero
    f = qconserv(a,i,b) + 1
    If (f > nocc) then
        imd = imd + t2(c,k,j) * eri(a,i)
    End If

    m = qconserv(i,a,j) + 1
    If (m <= nocc .and. m > 0) then
        imd = imd - t2(b,m,k) * eri(a,i)
    End If
    End Subroutine

    Subroutine buildW(qconserv, eri, t2, W, a, b, c, i, j, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: W
    Real (kind=pr)                  :: imd
    W = Zero
    Call build_imd(qconserv, eri, t2, imd, a, i, b, j, c, k, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, b, j, c, k, a, i, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, c, k, a, i, b, j, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, a, i ,c, k, b, j, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, b, j, a, i, c, k, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, c, k, b ,j, a, i, nocc, nao)
    W = W + imd
    End Subroutine

    Subroutine T_ene(rgvecs, qconserv, eri, t2, moe, ene, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
    Real (kind=pr),     Intent(in)  :: moe(nao), eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: ene
    Integer                         :: a,b,c,i,j,k
    Real (kind=pr)                  :: Waibjck, W, denom, fac
    Real (kind=pr)                  :: symm

    Call init_g_map(rgvecs, nao)
    ene = Zero

    !$omp parallel default(shared)
    !$omp do schedule(dynamic) private(Waibjck, W, denom, c, fac, symm) reduction(+:ene)
    Do a = nocc+1, nao
    Do b = nocc+1, a       
        Do i = 1, nocc     
        Do j = 1, nocc     
        Do k = 1, nocc     
            Call qconserv_c(rgvecs, i, j, k, a, b, c, nao)
            If (c <= nocc) cycle
            If (c > b) cycle
            denom = moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c)

            If (a == c) then
                denom = denom * 6.0_pr
            Else If (a == b .or. b == c) then
                denom = denom * 2.0_pr
            End If

            Call buildW(qconserv, eri, t2, W, a, b, c, i, j, k, nocc, nao)
            Waibjck = W / denom
            ene = ene + 4.0_pr * Waibjck * W

            Call buildW(qconserv, eri, t2, W, a, b, c, j, k, i, nocc, nao)
            ene = ene + Waibjck * W

            Call buildW(qconserv, eri, t2, W, a, b, c, k, i, j, nocc, nao)
            ene = ene + Waibjck * W

            Call buildW(qconserv, eri, t2, W, a, b, c, j, i, k, nocc, nao)
            ene = ene - 2.0_pr * Waibjck * W

            Call buildW(qconserv, eri, t2, W, a, b, c, i, k, j, nocc, nao)
            ene = ene - 2.0_pr * Waibjck * W

            Call buildW(qconserv, eri, t2, W, a, b, c, k, j, i, nocc, nao)
            ene = ene - 2.0_pr * Waibjck * W
        End Do
        End Do
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    ene = ene * 2.0_pr
    Deallocate(g_map)
    End Subroutine

    Subroutine buildDW(qconserv, t2, DW, a, i, b, j, c, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: DW
    Integer                         :: e, m
    DW = Zero
    e = qconserv(i,a,j) + 1
    If (e > nocc) then
        DW = DW + t2(a,i,j)
    End If

    m = qconserv(a,i,b) + 1
    If (m <= nocc .and. m > 0) then
        DW = DW - t2(a,i,m)
    End If
    End Subroutine

    Subroutine T_StructFac(rgvecs, qconserv, eri, t2, moe, Sq, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
    Real (kind=pr),     Intent(in)  :: moe(nao), eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: Sq(nocc+1:nao, nocc)
    Integer                         :: a,b,c,i,j,k
    Real (kind=pr)                  :: Waibjck, W, DW, denom

    Call init_g_map(rgvecs, nao)
    Sq = Zero
    Do a = nocc+1, nao
    Do b = nocc+1, nao
        Do i = 1, nocc
        Do j = 1, nocc
        Do k = 1, nocc
            Call qconserv_c(rgvecs, i, j, k, a, b, c, nao)
            If (c <= nocc) cycle
            denom = moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c)

            Call buildW(qconserv, eri, t2, W, a, b, c, i, j, k, nocc, nao)
            Waibjck = 4.0_pr * W
            Call buildW(qconserv, eri, t2, W, a, b, c, j, k, i, nocc, nao)
            Waibjck = Waibjck + W
            Call buildW(qconserv, eri, t2, W, a, b, c, k, i, j, nocc, nao)
            Waibjck = Waibjck + W
            Call buildW(qconserv, eri, t2, W, a, b, c, j, i, k, nocc, nao)
            Waibjck = Waibjck - 2.0_pr * W
            Call buildW(qconserv, eri, t2, W, a, b, c, i, k, j, nocc, nao)
            Waibjck = Waibjck - 2.0_pr * W
            Call buildW(qconserv, eri, t2, W, a, b, c, k, j, i, nocc, nao)
            Waibjck = Waibjck - 2.0_pr * W
            Waibjck = Waibjck / denom

            Call buildDW(qconserv, t2, DW, a, i, b, j, c, k, nocc, nao)
            Sq(c,k) = Sq(c,k) + DW * Waibjck

            Call buildDW(qconserv, t2, DW, b, j, a, i, c, k, nocc, nao)
            Sq(c,k) = Sq(c,k) + DW * Waibjck
            
            Call buildDW(qconserv, t2, DW, a, i, c, k, b, j, nocc, nao)
            Sq(b,j) = Sq(b,j) + DW * Waibjck

            Call buildDW(qconserv, t2, DW, c, k, a, i, b, j, nocc, nao)
            Sq(b,j) = Sq(b,j) + DW * Waibjck

            Call buildDW(qconserv, t2, DW, b, j, c, k, a, i, nocc, nao)
            Sq(a,i) = Sq(a,i) + DW * Waibjck

            Call buildDW(qconserv, t2, DW, c, k, b, j, a, i, nocc, nao)
            Sq(a,i) = Sq(a,i) + DW * Waibjck

        End Do
        End Do
        End Do
    End Do
    End Do

    Deallocate(g_map)
    End Subroutine

    End Module
