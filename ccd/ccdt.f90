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

    Subroutine build_imd(qconserv, eri, t2, imd, a, i, b, j, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, i, j
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: imd
    Integer                         :: c,d,l
    imd = Zero
    d = qconserv(i,a,j) + 1
    If (d > nocc) then
        imd = imd + t2(a,i,j) * eri(b,d)
    End If

    l = qconserv(a,i,b) + 1
    If (l <= nocc .and. l > 0) then
        imd = imd - t2(a,i,l) * eri(l,j)
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
    Call build_imd(qconserv, eri, t2, imd, a, i, b, j, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, b, j, c, k, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, c, k, a, i, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, a, i ,c, k, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, b, j, a, i, nocc, nao)
    W = W + imd
    Call build_imd(qconserv, eri, t2, imd, c, k, b ,j, nocc, nao)
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
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: Waibjck, W, denom
    Call init_g_map(rgvecs, nao)
    ene = Zero
    Do a = nocc+1, nao
    Do b = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
    Do k = 1, nocc
        Call qconserv_c(rgvecs, i, j, k, a, b, c, nao)
        If (c <= nocc) cycle
        denom = moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c)
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
    ene = ene / 3.0_pr
    End Subroutine

    End Module
