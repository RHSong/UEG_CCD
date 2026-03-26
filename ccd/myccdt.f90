Module UEG_CCD_T
Use Precision
Use Constants
Implicit None

Contains

Subroutine build_w_tensor(qconserv, eri, t2, wt, a, b, c, nocc, nao)
Implicit None
Integer,            Intent(in)  :: nocc, nao
Integer,            Intent(in)  :: a, b, c
Integer,            Intent(in)  :: qconserv(nao,nao,nao)
Real (kind=pr),     Intent(in)  :: eri(nao,nao)
Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
Real (kind=pr),     Intent(out) :: wt(nocc,nocc,nocc)
Integer                         :: i, j, k, f, m
Real (kind=pr)                  :: eai

wt = Zero
Do i = 1, nocc
    Do j = 1, nocc
        Do k = 1, nocc
            eai = eri(a,i)

            f = qconserv(a,i,b) + 1
            If (f > nocc .and. f <= nao) then
                If (qconserv(k,c,j) + 1 == f) then
                    wt(i,j,k) = wt(i,j,k) + t2(c,k,j) * eai
                End If
            End If

            m = qconserv(a,i,j) + 1
            If (m >= 1 .and. m <= nocc) then
                If (qconserv(m,b,k) + 1 == c) then
                    wt(i,j,k) = wt(i,j,k) - t2(b,m,k) * eai
                End If
            End If
        End Do
    End Do
End Do
End Subroutine

Subroutine build_z_tensor(wt, zt, moe, a, b, c, nocc, nao)
Implicit None
Integer,            Intent(in)  :: nocc, nao
Integer,            Intent(in)  :: a, b, c
Real (kind=pr),     Intent(in)  :: wt(nocc,nocc,nocc)
Real (kind=pr),     Intent(in)  :: moe(nao)
Real (kind=pr),     Intent(out) :: zt(nocc,nocc,nocc)
Integer                         :: i, j, k
Real (kind=pr)                  :: d3, symm

symm = One
If (a == c) then
    symm = 6.0_pr
Else If (a == b .or. b == c) then
    symm = 2.0_pr
End If

Do i = 1, nocc
    Do j = 1, nocc
        Do k = 1, nocc
            d3 = symm * (moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c))
            zt(i,j,k) = (4.0_pr * wt(i,j,k) + wt(j,k,i) + wt(k,i,j) &
                       - 2.0_pr * wt(k,j,i) - 2.0_pr * wt(i,k,j) - 2.0_pr * wt(j,i,k)) / d3
        End Do
    End Do
End Do
End Subroutine

Subroutine T_ene(rgvecs, qconserv, eri, t2, moe, ene, nocc, nao)
Implicit None
Integer,            Intent(in)  :: nocc, nao
Integer,            Intent(in)  :: qconserv(nao,nao,nao)
Real (kind=pr),     Intent(in)  :: rgvecs(nao, 3)
Real (kind=pr),     Intent(in)  :: moe(nao), eri(nao,nao)
Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
Real (kind=pr),     Intent(out) :: ene
Integer                         :: a, b, c, i, j, k
Real (kind=pr), allocatable     :: wabc(:,:,:), wacb(:,:,:), wbac(:,:,:), wbca(:,:,:), wcab(:,:,:), wcba(:,:,:)
Real (kind=pr), allocatable     :: zabc(:,:,:), zacb(:,:,:), zbac(:,:,:), zbca(:,:,:), zcab(:,:,:), zcba(:,:,:)

allocate(wabc(nocc,nocc,nocc), wacb(nocc,nocc,nocc), wbac(nocc,nocc,nocc))
allocate(wbca(nocc,nocc,nocc), wcab(nocc,nocc,nocc), wcba(nocc,nocc,nocc))
allocate(zabc(nocc,nocc,nocc), zacb(nocc,nocc,nocc), zbac(nocc,nocc,nocc))
allocate(zbca(nocc,nocc,nocc), zcab(nocc,nocc,nocc), zcba(nocc,nocc,nocc))

ene = Zero
Do a = nocc + 1, nao
    Do b = nocc + 1, a
        Do c = nocc + 1, b
            Call build_w_tensor(qconserv, eri, t2, wabc, a, b, c, nocc, nao)
            Call build_w_tensor(qconserv, eri, t2, wacb, a, c, b, nocc, nao)
            Call build_w_tensor(qconserv, eri, t2, wbac, b, a, c, nocc, nao)
            Call build_w_tensor(qconserv, eri, t2, wbca, b, c, a, nocc, nao)
            Call build_w_tensor(qconserv, eri, t2, wcab, c, a, b, nocc, nao)
            Call build_w_tensor(qconserv, eri, t2, wcba, c, b, a, nocc, nao)

            Call build_z_tensor(wabc, zabc, moe, a, b, c, nocc, nao)
            Call build_z_tensor(wacb, zacb, moe, a, b, c, nocc, nao)
            Call build_z_tensor(wbac, zbac, moe, a, b, c, nocc, nao)
            Call build_z_tensor(wbca, zbca, moe, a, b, c, nocc, nao)
            Call build_z_tensor(wcab, zcab, moe, a, b, c, nocc, nao)
            Call build_z_tensor(wcba, zcba, moe, a, b, c, nocc, nao)

            Do i = 1, nocc
                Do j = 1, nocc
                    Do k = 1, nocc
                        ene = ene + wabc(i,j,k) * zabc(i,j,k)
                        ene = ene + wacb(i,k,j) * zabc(i,j,k)
                        ene = ene + wbac(j,i,k) * zabc(i,j,k)
                        ene = ene + wbca(j,k,i) * zabc(i,j,k)
                        ene = ene + wcab(k,i,j) * zabc(i,j,k)
                        ene = ene + wcba(k,j,i) * zabc(i,j,k)

                        ene = ene + wacb(i,j,k) * zacb(i,j,k)
                        ene = ene + wabc(i,k,j) * zacb(i,j,k)
                        ene = ene + wcab(j,i,k) * zacb(i,j,k)
                        ene = ene + wcba(j,k,i) * zacb(i,j,k)
                        ene = ene + wbac(k,i,j) * zacb(i,j,k)
                        ene = ene + wbca(k,j,i) * zacb(i,j,k)

                        ene = ene + wbac(i,j,k) * zbac(i,j,k)
                        ene = ene + wbca(i,k,j) * zbac(i,j,k)
                        ene = ene + wabc(j,i,k) * zbac(i,j,k)
                        ene = ene + wacb(j,k,i) * zbac(i,j,k)
                        ene = ene + wcba(k,i,j) * zbac(i,j,k)
                        ene = ene + wcab(k,j,i) * zbac(i,j,k)

                        ene = ene + wbca(i,j,k) * zbca(i,j,k)
                        ene = ene + wbac(i,k,j) * zbca(i,j,k)
                        ene = ene + wcba(j,i,k) * zbca(i,j,k)
                        ene = ene + wcab(j,k,i) * zbca(i,j,k)
                        ene = ene + wabc(k,i,j) * zbca(i,j,k)
                        ene = ene + wacb(k,j,i) * zbca(i,j,k)

                        ene = ene + wcab(i,j,k) * zcab(i,j,k)
                        ene = ene + wcba(i,k,j) * zcab(i,j,k)
                        ene = ene + wacb(j,i,k) * zcab(i,j,k)
                        ene = ene + wabc(j,k,i) * zcab(i,j,k)
                        ene = ene + wbca(k,i,j) * zcab(i,j,k)
                        ene = ene + wbac(k,j,i) * zcab(i,j,k)

                        ene = ene + wcba(i,j,k) * zcba(i,j,k)
                        ene = ene + wcab(i,k,j) * zcba(i,j,k)
                        ene = ene + wbca(j,i,k) * zcba(i,j,k)
                        ene = ene + wbac(j,k,i) * zcba(i,j,k)
                        ene = ene + wacb(k,i,j) * zcba(i,j,k)
                        ene = ene + wabc(k,j,i) * zcba(i,j,k)
                    End Do
                End Do
            End Do
        End Do
    End Do
End Do

ene = 2.0_pr * ene

deallocate(wabc, wacb, wbac, wbca, wcab, wcba)
deallocate(zabc, zacb, zbac, zbca, zcab, zcba)
End Subroutine

End Module
