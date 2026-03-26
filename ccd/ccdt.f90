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

    !Subroutine build_imd(qconserv, eri, t2, imd, a, i, b, j, c, k, nocc, nao)
    !Implicit None
    !Integer,            Intent(in)  :: nocc, nao
    !Integer,            Intent(in)  :: a, b, c, i, j, k
    !Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    !Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    !Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    !Real (kind=pr),     Intent(out) :: imd
    !Integer                         :: f, m, tmp

    !imd = Zero

    !! --- Term 1: \sum_f t2(a,i,j) * (bc|fk) ---
    !! T2 implied virtual f = k_i + k_j - k_a
    !f = qconserv(i,a,j) + 1
    !! ERI required virtual f = k_b - k_c + k_k
    !tmp = qconserv(b,c,k) + 1
    !If (f > nocc .and. f <= nao .and. f == tmp) then
    !    imd = imd + t2(a,i,j) * eri(b,c)
    !End If

    !! --- Term 2: \sum_m t2(a,i,m) * (mc|jk) ---
    !! T2 implied occupied m = k_a + k_b - k_i
    !m = qconserv(a,i,b) + 1
    !! ERI required occupied m = k_c + k_j - k_k
    !tmp = qconserv(c,k,j) + 1
    !If (m <= nocc .and. m > 0 .and. m == tmp) then
    !    imd = imd - t2(a,i,m) * eri(j,k)
    !End If
    !End Subroutine

    Subroutine build_imd(qconserv, eri, t2, imd, a, i, b, j, c, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: imd
    Integer                         :: f, m, tmp
    imd = Zero
    f = qconserv(a,i,b) + 1
    tmp = qconserv(k,c,j) + 1
    If (f > nocc .and. f == tmp) then
        imd = imd + t2(c,k,j) * eri(a,i)
    End If

    m = qconserv(a,i,j) + 1
    tmp = qconserv(b,k,c) + 1
    If (m <= nocc .and. m > 0 .and. m == tmp) then
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
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: Waibjck, W, denom
    !Call init_g_map(rgvecs, nao)
    ene = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(Waibjck, W, denom) reduction(+:ene)
    Do a = nocc+1, nao
    Do b = nocc+1, nao
    Do c = nocc+1, nao
        Do i = 1, nocc
        Do j = 1, nocc
        Do k = 1, nocc
            !Call qconserv_c(rgvecs, i, j, k, a, b, c, nao)
            !If (c <= nocc) cycle
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
    End Do
    !$omp end do
    !$omp end parallel
    ene = ene / 3.0_pr
    End Subroutine

    End Module
