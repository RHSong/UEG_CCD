    Module UEG_CCD
    Use Precision
    Use Constants
    Implicit None

    Contains

    Subroutine build_denom(qconserv, moe, denom, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: moe(nao)
    Real (kind=pr),     Intent(out) :: denom(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    denom = 1e8_pr
    !$omp parallel default(shared)
    !$omp do schedule(static) private(b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        Do j = 1, nocc
            b = qconserv(i,a,j) + 1
            If (b <= nocc) cycle
            denom(a,i,j) = moe(a) + moe(b) - moe(i) - moe(j)
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    Subroutine CCD_ene(qconserv, eri, t2, ed, ex, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
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
            b = qconserv(i,a,j) + 1
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

    Subroutine CCD_Res(qconserv, eri, denom, t2, r2, nocc, nao, &
                       DoMosaic, DoRing, DoLadder, DoxRing)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
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
            b = qconserv(i,a,j) + 1
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
        Call Mosaic(qconserv, eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If

    If (DoRing) then
        Call Ring(qconserv, eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If

    If (DoLadder) then
        Call Ladder(qconserv, eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If

    If (DoxRing) then
        Call xRing(qconserv, eri, t2, r2tmp, nocc, nao)
        r2 = r2 + r2tmp
    End If
    End Subroutine

    Subroutine CCFock(qconserv, eri, t2, Foo, Fvv, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
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
            d = qconserv(k,c,l) + 1
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
            d = qconserv(k,c,l) + 1
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

    Subroutine Mosaic(qconserv, eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: kcld, kdlc
    Real (kind=pr)                  :: Foo(nocc), Fvv(nocc+1:nao)
    r2 = Zero
    Call CCFock(qconserv, eri, t2, Foo, Fvv, nocc, nao)
    !$omp parallel default(shared)
    !$omp do schedule(static) private(b)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        b = qconserv(i,a,j) + 1
        If (b <= nocc) cycle
        r2(a,i,j) = r2(a,i,j) + (Fvv(a) + Fvv(b) - Foo(i) - Foo(i)) * (2.0_pr * t2(a,i,j) - t2(a,j,i))
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine Ring(qconserv, eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: bjkc, bckj, aikc, acki
    Real (kind=pr)                  :: kcld, kdlc
    r2 = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(bjkc, bckj, aikc, acki), &
    !$omp& private(kcld, kdlc, b, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        b = qconserv(i,a,j) + 1
        If (b <= nocc) cycle

        Do k = 1, nocc
            c = qconserv(i,a,k) + 1
            If (c <= nocc) cycle
            bjkc = eri(b,j)
            bckj = eri(b,c)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * bjkc - bckj) * (2.0_pr * t2(a,i,k) - t2(a,k,i))
        End Do

        Do k = 1, nocc
            c = qconserv(a,i,k) + 1
            If (c <= nocc) cycle
            aikc = eri(a,i)
            acki = eri(a,c)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * aikc - acki) * (2.0_pr * t2(b,j,k) - t2(b,k,j))
        End Do

        Do k = 1, nocc
            c = qconserv(i,a,k) + 1
            If (c <= nocc) cycle
            Do l = 1, nocc
                d = qconserv(j,b,l) + 1
                If (d <= nocc) cycle
                kcld = eri(k,c)
                kdlc = eri(k,d)
                r2(a,i,j) = r2(a,i,j) + (2.0_pr * kcld - kdlc) * (2.0_pr * t2(a,i,k) - t2(a,k,i)) * (2.0_pr * t2(b,j,l) - t2(b,l,j))
            End Do
        End Do
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine Ladder(qconserv, eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: acbd, adbc
    Real (kind=pr)                  :: kilj, kjli
    Real (kind=pr)                  :: kcld, kdlc
    r2 = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(acbd, adbc, kilj, kjli), &
    !$omp& private(kcld, kdlc, b, d, l)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        b = qconserv(i,a,j) + 1
        If (b <= nocc) cycle

        Do c = nocc+1, nao
            d = qconserv(a,c,b) + 1
            If (d <= nocc) cycle
            acbd = eri(a,c)
            adbc = eri(a,d)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * acbd - adbc) * t2(c,i,j)
        End Do

        Do k = 1, nocc
            l = qconserv(i,k,j) + 1
            If (l > nocc .or. l == 0)  cycle
            kilj = eri(k,i)
            kjli = eri(k,j)
            r2(a,i,j) = r2(a,i,j) + (2.0_pr * kilj - kjli) * t2(a,k,l)
        End Do

        Do k = 1, nocc
            l = qconserv(a,k,b) + 1
            If (l > nocc .or. l == 0) cycle
            Do c = nocc+1, nao
                d = qconserv(k,c,l) + 1
                If (d <= nocc) cycle
                kcld = eri(k,c)
                kdlc = eri(k,d)
                r2(a,i,j) = r2(a,i,j) + (2.0_pr * kcld - kdlc) * t2(a,k,l) * t2(c,i,j)
            End Do
        End Do
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine xRing(qconserv, eri, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: bikc, bcki, ajkc, ackj
    Real (kind=pr)                  :: kcld, kdlc
    r2 = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(bikc, bcki, ajkc, ackj), &
    !$omp& private(kcld, kdlc, b, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        b = qconserv(i,a,j) + 1
        If (b <= nocc) cycle

        Do k = 1, nocc
            c = qconserv(b,i,k) + 1
            If (c <= nocc) cycle
            bikc = eri(b,i)
            bcki = eri(b,c)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * bikc - bcki) * t2(a,j,k)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * bcki - bikc) * t2(a,k,j)
        End Do

        Do k = 1, nocc
            c = qconserv(a,j,k) + 1
            If (c <= nocc) cycle
            ajkc = eri(a,j)
            ackj = eri(a,c)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * ajkc - ackj) * t2(b,i,k)
            r2(a,i,j) = r2(a,i,j) - (2.0_pr * ackj - ajkc) * t2(b,k,i)
        End Do

        Do k = 1, nocc
            c = qconserv(j,a,k) + 1
            If (c <= nocc) cycle
            Do l = 1, nocc
                d = qconserv(i,b,l) + 1
                If (d <= nocc) cycle
                kcld = eri(k,c)
                kdlc = eri(k,d)
                r2(a,i,j) = r2(a,i,j) + (2.0_pr * kcld - kdlc) * (t2(a,j,k) * t2(b,l,i) + t2(a,k,j) * t2(b,i,l))
                r2(a,i,j) = r2(a,i,j) + (2.0_pr * kdlc - kcld) * t2(a,k,j) * t2(b,l,i)
                r2(a,i,j) = r2(a,i,j) - 2.0_pr * (2.0_pr * kcld - kdlc) * t2(a,j,k) * t2(b,i,l)
            End Do
        End Do
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    End Subroutine

    Subroutine drCCD_res(qconserv, eri, denom, t2, r2, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: qconserv(nao,nao,nao)
    Real (kind=pr),     Intent(in)  :: eri(nao,nao), denom(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: r2(nocc+1:nao,nocc,nocc)
    Integer                         :: a,b,c,d,i,j,k,l
    Real (kind=pr)                  :: aibj, bjkc, aikc, kcld
    r2 = Zero
    !$omp parallel default(shared)
    !$omp do schedule(static) private(aibj, b)
    Do a = nocc+1, nao
    Do i = 1, nocc
        aibj = eri(a,i)
        Do j = 1, nocc
            b = qconserv(i,a,j) + 1
            If (b <= nocc) cycle
            r2(a,i,j) = r2(a,i,j) + aibj
            r2(a,i,j) = r2(a,i,j) + t2(a,i,j) * denom(a,i,j)
        End Do
    End Do
    End Do
    !$omp end do

    !$omp do schedule(static) private(bjkc, aikc, kcld, b, c, d)
    Do a = nocc+1, nao
    Do i = 1, nocc
    Do j = 1, nocc
        b = qconserv(i,a,j) + 1
        If (b <= nocc) cycle

        Do k = 1, nocc
            c = qconserv(i,a,k) + 1
            If (c <= nocc) cycle
            bjkc = eri(b,j)
            r2(a,i,j) = r2(a,i,j) + 2.0_pr * bjkc * t2(a,i,k)
        End Do

        Do k = 1, nocc
            c = qconserv(a,i,k) + 1
            If (c <= nocc) cycle
            aikc = eri(a,i)
            r2(a,i,j) = r2(a,i,j) + 2.0_pr * aikc * t2(b,j,k)
        End Do

        Do k = 1, nocc
            c = qconserv(i,a,k) + 1
            If (c <= nocc) cycle
            Do l = 1, nocc
                d = qconserv(j,b,l) + 1
                If (d <= nocc) cycle
                kcld = eri(k,c)
                r2(a,i,j) = r2(a,i,j) + 4.0_pr * kcld * t2(a,i,k) * t2(b,j,l)
            End Do
        End Do
    End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    End Module
