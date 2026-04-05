    Module UEG_CCD
    Use Precision
    Use Constants
    Use qconserv
    Implicit None

    Contains

    Subroutine build_denom(moe, denom, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: moe(nao)
    Real (kind=pr),     Intent(out) :: denom(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    denom = 1e8_pr
    !$omp parallel default(shared)
    !$omp do schedule(static) private(b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        Do j = 1, nocc
            Call qconserv2(i, a, j, b)
            If (b <= nocc) cycle
            denom(a,i,j) = moe(a) + moe(b) - moe(i) - moe(j)
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    Subroutine CCD_ene(eri, t2, ed, ex, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: ed, ex
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: iajb, ibja
    ed = Zero
    ex = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) reduction(+:ed, ex), &
    !$omp& private(iajb, ibja, b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        iajb = eri(i,a)
        Do j = 1, nocc
            Call qconserv2(i, a, j, b)
            ibja = eri(i,b)
            If (b <= nocc) cycle
            
            ed = ed + 2.0_pr * iajb * t2(a,i,j)
            ex = ex - ibja * t2(a,i,j)
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    Subroutine CCD_Res(eri, denom, t2, r2, nocc, nao, &
                       DoMosaic, DoRing, DoLadder, DoxRing)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao), denom(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Logical,            Intent(in)  :: DoMosaic, DoRing
    Logical,            Intent(in)  :: DoLadder, DoxRing
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: aibj, ajbi, r2tmp(nocc+1:nao,nocc,nocc)
    r2 = Zero
    ! always include driving terms
    !$omp parallel default(shared)
    !$omp do schedule(static) private(aibj, ajbi, b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        aibj = eri(a,i)
        Do j = 1, nocc
            ajbi = eri(a,j)
            Call qconserv2(i, a, j, b)
            If (b <= nocc) cycle
            r2(a,i,j) = r2(a,i,j) + 2.0_pr * aibj - ajbi
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * t2(a,i,j) - t2(a,j,i)) * denom(a,i,j)
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    ! different channels
    If (DoMosaic) then
        Call Mosaic(eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If

    If (DoRing) then
        Call Ring(eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If

    If (DoLadder) then
        Call Ladder(eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If

    If (DoxRing) then
        Call xRing(eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If
    End Subroutine

    Subroutine CCFock(eri, t2, Foo, Fvv, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: Foo(nocc), Fvv(nocc+1:nao)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: kcld, kdlc
    Foo = Zero
    Fvv = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(kcld, kdlc, d)
    Do k = 1, nocc
        Do c = nocc+1, nao
        Do l = 1, nocc
            Call qconserv2(k, c, l, d)
            If (d <= nocc) cycle
            kcld = eri(k,c)
            kdlc = eri(k,d)
            Foo(k) = Foo(k) + t2(c,k,l) * (2.0_pr * kcld - kdlc)
        End Do
        End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(kcld, kdlc, d)
    Do c = nocc+1, nao
        Do k = 1, nocc
        Do l = 1, nocc
            Call qconserv2(k, c, l, d)
            If (d <= nocc) cycle
            kcld = eri(k,c)
            kdlc = eri(k,d)
            Fvv(c) = Fvv(c) - t2(c,k,l) * (2.0_pr * kcld - kdlc)
        End Do
        End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    Subroutine Mosaic(eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: kcld, kdlc
    Real (kind=pr)                  :: Foo(nocc), Fvv(nocc+1:nao)
    r2 = Zero
    Call CCFock(eri, t2, Foo, Fvv, nocc, nao)
    !$omp parallel default(shared)
    !$omp do schedule(static) private(b)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        Call qconserv2(i, a, j, b)
        If (b <= nocc) cycle
        r2(a,i,j) = r2(a,i,j) + (Fvv(a) + Fvv(b) - Foo(i) - Foo(i)) * (2.0_pr * t2(a,i,j) - t2(a,j,i))
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine Ring(eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: bjkc, bckj, aikc, acki
    Real (kind=pr)                  :: kcld, kdlc
    Real (kind=pr)                  :: imd(nocc+1:nao,nocc,nocc)
    r2 = Zero
    imd = Zero

    !$omp parallel default(shared)
    !$omp do schedule(static) private(kcld, kdlc, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do k = 1, nocc
        Call qconserv2(i, a, k, c)
        If (c <= nocc) cycle
        Do l = 1, nocc
            Call qconserv2(k, c, l, d)
            If (d <= nocc) cycle
            kcld = eri(k,c)
            kdlc = eri(k,d)
            imd(a,i,l) = imd(a,i,l) + (2.0_pr * kcld - kdlc) * (2.0_pr * t2(a,i,k) - t2(a,k,i))
        End Do
    End Do
    End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(bjkc, bckj, aikc, acki), &
    !$omp& private(b, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        Call qconserv2(i, a, j, b)
        If (b <= nocc) cycle

        Do k = 1, nocc
            Call qconserv2(i, a, k, c)
            If (c <= nocc) cycle
            bjkc = eri(b,j)
            bckj = eri(b,c)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * bjkc - bckj) * (2.0_pr * t2(a,i,k) - t2(a,k,i))
        End Do

        Do k = 1, nocc
            Call qconserv2(a, i, k, c)
            If (c <= nocc) cycle
            aikc = eri(a,i)
            acki = eri(a,c)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * aikc - acki + imd(a,i,k)) * (2.0_pr * t2(b,j,k) - t2(b,k,j))
        End Do

    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine Ladder(eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: acbd, adbc
    Real (kind=pr)                  :: kilj, kjli
    Real (kind=pr)                  :: kcld, kdlc
    Real (kind=pr)                  :: imd(nocc,nocc,nocc)
    r2 = Zero
    imd = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(kcld, kdlc, d, l)
    Do i = 1, nocc
    Do j = 1, nocc
    Do c = nocc+1, nao
        Call qconserv2(i, c, j, d)
        If (d <= nocc) cycle
        Do k = 1, nocc
            Call qconserv2(c, k, d, l)
            If (l > nocc .or. l == 0) cycle
            kcld = eri(k,c)
            kdlc = eri(k,d)
            imd(k,i,j) = imd(k,i,j) + (2.0_pr * kcld - kdlc) * t2(c,i,j)
        End Do
    End Do
    End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(acbd, adbc, kilj, kjli), &
    !$omp& private(b, d, l)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        Call qconserv2(i, a, j, b)
        If (b <= nocc) cycle

        Do c = nocc+1, nao
            Call qconserv2(a, c, b, d)
            If (d <= nocc) cycle
            acbd = eri(a,c)
            adbc = eri(a,d)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * acbd - adbc) * t2(c,i,j)
        End Do

        Do k = 1, nocc
            Call qconserv2(i, k, j, l)
            If (l > nocc .or. l == 0)  cycle
            kilj = eri(k,i)
            kjli = eri(k,j)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * kilj - kjli + imd(k,i,j)) * t2(a,k,l)
        End Do

    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine xRing(eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: bikc, bcki, ajkc, ackj
    Real (kind=pr)                  :: kcld, kdlc
    Real (kind=pr)                  :: imd1(nocc+1:nao,nocc,nocc)
    Real (kind=pr)                  :: imd2(nocc+1:nao,nocc,nocc)
    r2 = Zero
    imd1 = Zero
    imd2 = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(kcld, kdlc, c, d)
    Do b = nocc+1, nao
    Do i = 1, nocc
    Do l = 1, nocc
        Call qconserv2(i, b, l, d)
        If (d <= nocc) cycle
        Do k = 1, nocc
            Call qconserv2(k, d, l, c)
            If (c <= nocc) cycle
            kcld = eri(k,c)
            kdlc = eri(k,d)
            imd1(b,i,k) = imd1(b,i,k) + (2.0_pr * kcld - kdlc) * t2(b,l,i)
            imd1(b,i,k) = imd1(b,i,k) - 2.0_pr * (2.0_pr * kcld - kdlc) * t2(b,i,l)
            imd2(b,k,i) = imd2(b,k,i) + (2.0_pr * kdlc - kcld) * t2(b,l,i)
            imd2(b,k,i) = imd2(b,k,i) + (2.0_pr * kcld - kdlc) * t2(b,i,l)
        End Do
    End Do
    End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(bikc, bcki, ajkc, ackj), &
    !$omp& private(b, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        Call qconserv2(i, a, j, b)
        If (b <= nocc) cycle

        Do k = 1, nocc
            Call qconserv2(b, i, k, c)
            If (c <= nocc) cycle
            bikc = eri(b,i)
            bcki = eri(b,c)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * bikc - bcki - imd1(b,i,k)) * t2(a,j,k)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * bcki - bikc - imd2(b,k,i)) * t2(a,k,j)
        End Do

        Do k = 1, nocc
            Call qconserv2(a, j, k, c)
            If (c <= nocc) cycle
            ajkc = eri(a,j)
            ackj = eri(a,c)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * ajkc - ackj) * t2(b,i,k)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * ackj - ajkc) * t2(b,k,i)
        End Do

    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine drCCD_res(eri, denom, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: eri(nao,nao), denom(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: aibj, bjkc, aikc, kcld
    Real (kind=pr)                  :: imd(nocc+1:nao,nocc,nocc)
    r2 = Zero
    imd = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(aibj, b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        aibj = eri(a,i)
        Do j = 1, nocc
            Call qconserv2(i, a, j, b)
            If (b <= nocc) cycle
            r2(a,i,j) = r2(a,i,j) + aibj
            r2(a,i,j) = r2(a,i,j) + t2(a,i,j) * denom(a,i,j)
        End Do
    End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(kcld, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do k = 1, nocc
        Call qconserv2(i, a, k, c)
        If (c <= nocc) cycle
        Do l = 1, nocc
            Call qconserv2(k, c, l, d)
            If (d <= nocc) cycle
            kcld = eri(k,c)
            imd(a,i,l) = imd(a,i,l) + 2.0_pr * kcld * t2(a,i,k)
        End Do
    End Do
    End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(bjkc, aikc, b, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        Call qconserv2(i, a, j, b)
        If (b <= nocc) cycle

        Do k = 1, nocc
            Call qconserv2(i, a, k, c)
            If (c <= nocc) cycle
            bjkc = eri(b,j)
            r2(a,i,j) = r2(a,i,j) + 2.0_pr * bjkc * t2(a,i,k)
        End Do

        Do k = 1, nocc
            Call qconserv2(a, i, k, c)
            If (c <= nocc) cycle
            aikc = eri(a,i)
            r2(a,i,j) = r2(a,i,j) + 2.0_pr * (aikc + imd(a,i,k)) * t2(b,j,k)
        End Do
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    Subroutine strucfac_t2(t2, Sqd, Sqx, Lq, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: Sqd(nocc+1:nao,nocc)
    Real (kind=pr),     Intent(out) :: Sqx(nocc+1:nao,nocc)
    Real (kind=pr),     Intent(out) :: Lq(nocc+1:nao,nocc,3)
    Integer                         :: a,b,c,d,i,j,k,l
    Sqd = Zero
    Sqx = Zero
    Lq = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        Lq(a,i,:) = gvecs(a,:) - gvecs(i,:)
        Do j = 1, nocc
            Call qconserv2(i, a, j, b)
            If (b <= nocc) cycle
            Sqd(a,i) = Sqd(a,i) + 2 * t2(a,i,j)
            Sqx(a,i) = Sqx(a,i) - t2(a,j,i)
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    End Module

