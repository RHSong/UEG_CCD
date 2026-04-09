    Module UEG_CCD_T
    Use Precision
    Use Constants
    Use qconserv
    Implicit None

    Contains

    Subroutine build_imd(eri, t2, imd, a, i, b, j, c, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: imd
    Integer                         :: f, m
    imd = Zero
    Call qconserv2(a, i, b, f)
    If (f > nocc) then
        imd = imd + t2(c,k,j) * eri(a,i)
    End If

    Call qconserv2(i, a, j, m)
    If (m <= nocc .and. m > 0) then
        imd = imd - t2(b,m,k) * eri(a,i)
    End If
    End Subroutine

    Subroutine buildW(eri, t2, W, a, b, c, i, j, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Real (kind=pr),     Intent(in)  :: eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: W
    Real (kind=pr)                  :: imd
    W = Zero
    Call build_imd(eri, t2, imd, a, i, b, j, c, k, nocc, nao)
    W = W + imd
    Call build_imd(eri, t2, imd, b, j, c, k, a, i, nocc, nao)
    W = W + imd
    Call build_imd(eri, t2, imd, c, k, a, i, b, j, nocc, nao)
    W = W + imd
    Call build_imd(eri, t2, imd, a, i ,c, k, b, j, nocc, nao)
    W = W + imd
    Call build_imd(eri, t2, imd, b, j, a, i, c, k, nocc, nao)
    W = W + imd
    Call build_imd(eri, t2, imd, c, k, b ,j, a, i, nocc, nao)
    W = W + imd
    End Subroutine

    Subroutine T_ene(eri, t2, moe, ene, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: moe(nao), eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: ene
    Integer                         :: a, b, c, i, j, k
    Real (kind=pr)                  :: denom
    Real (kind=pr)                  :: W1, W2, W3, W4, W5, W6
    Real (kind=pr)                  :: Wsq_sum, Seven, Sodd
    Real (kind=pr)                  :: Etot

    ene = Zero

    !$omp parallel default(shared) &
    !$omp private(b, c, i, j, k, denom, &
    !$omp         W1, W2, W3, W4, W5, W6, &
    !$omp         Seven, Sodd, Wsq_sum, Etot)
    !$omp do schedule(dynamic) reduction(+:ene)
    Do a = nocc+1, nao
    Do b = nocc+1, a
        Do i = 1, nocc
        Do j = 1, i 
        Do k = 1, j 
            Call qconserv_c(i, j, k, a, b, c)
            If (c <= nocc) cycle
            If (c > b) cycle
            denom = moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c)

            If (a == c) then
                denom = denom * 6.0_pr
            Else If (a == b .or. b == c) then
                denom = denom * 2.0_pr
            End If

            If (i == k) then
                denom = denom * 6.0_pr
            Else If (i == j .or. j == k) then
                denom = denom * 2.0_pr
            End If

            Call buildW(eri, t2, W1, a, b, c, i, j, k, nocc, nao)
            Call buildW(eri, t2, W2, a, b, c, j, k, i, nocc, nao)
            Call buildW(eri, t2, W3, a, b, c, k, i, j, nocc, nao)
            Call buildW(eri, t2, W4, a, b, c, j, i, k, nocc, nao)
            Call buildW(eri, t2, W5, a, b, c, i, k, j, nocc, nao)
            Call buildW(eri, t2, W6, a, b, c, k, j, i, nocc, nao)

            Seven = W1 + W2 + W3
            Sodd  = W4 + W5 + W6
            Wsq_sum = W1*W1 + W2*W2 + W3*W3 + W4*W4 + W5*W5 + W6*W6

            ene = ene + (3.0_pr * Wsq_sum + Seven*Seven + Sodd*Sodd - 4.0_pr * Seven * Sodd) / denom
        End Do
        End Do
        End Do
    End Do
    End Do
    !$omp end do
    !$omp end parallel

    ene = ene * 2.0_pr
    End Subroutine

    Subroutine buildDW(t2, DW, a, i, b, j, c, k, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Integer,            Intent(in)  :: a, b, c, i, j, k
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: DW
    Integer                         :: f, m
    DW = Zero
    Call qconserv2(i, a, j, m)
    If (m <= nocc .and. m > 0) then
        DW = DW - t2(b,m,k)
    End If

    Call qconserv2(a, i, b, f)
    If (f > nocc) then
        DW = DW + t2(c,k,j)
    End If
    End Subroutine

    !Subroutine StructFac_T3(eri, t2, moe, Sq, nocc, nao)
    !Implicit None
    !Integer,            Intent(in)  :: nocc, nao
    !Real (kind=pr),     Intent(in)  :: moe(nao), eri(nao,nao)
    !Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    !Real (kind=pr),     Intent(out) :: Sq(nocc+1:nao, nocc)
    !Integer                         :: a,b,c,i,j,k
    !Real (kind=pr)                  :: Waibjck, W, DW, denom
    !                                                                    
    !Sq = Zero
    !                                                                    
    !!$omp parallel default(shared)
    !!$omp do schedule(dynamic) reduction(+:Sq) &
    !!$omp& private(c, denom, W, Waibjck, DW)
    !Do a = nocc+1, nao
    !Do b = nocc+1, a
    !    Do i = 1, nocc
    !    Do j = 1, nocc
    !    Do k = 1, nocc
    !        Call qconserv_c(i, j, k, a, b, c)
    !        If (c <= nocc) cycle
    !        If (c > b) cycle
    !        denom = moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c)
    !                                                                    
    !        If (a == c) then
    !            denom = denom
    !        Else If (a == b .or. b == c) then
    !            denom = denom / 3.0_pr
    !        Else
    !            denom = denom / 6.0_pr
    !        End If
    !                                                                    
    !        Call buildW(eri, t2, W, a, b, c, i, j, k, nocc, nao)
    !        Waibjck = 4.0_pr * W
    !        Call buildW(eri, t2, W, a, b, c, j, k, i, nocc, nao)
    !        Waibjck = Waibjck + W
    !        Call buildW(eri, t2, W, a, b, c, k, i, j, nocc, nao)
    !        Waibjck = Waibjck + W
    !        Call buildW(eri, t2, W, a, b, c, j, i, k, nocc, nao)
    !        Waibjck = Waibjck - 2.0_pr * W
    !        Call buildW(eri, t2, W, a, b, c, i, k, j, nocc, nao)
    !        Waibjck = Waibjck - 2.0_pr * W
    !        Call buildW(eri, t2, W, a, b, c, k, j, i, nocc, nao)
    !        Waibjck = Waibjck - 2.0_pr * W
    !        Waibjck = Waibjck / denom
    !                                                                    
    !        Call buildDW(t2, DW, a, i, b, j, c, k, nocc, nao)
    !        Sq(a,i) = Sq(a,i) + DW * Waibjck
    !                                                          
    !        Call buildDW(t2, DW, b, j, c, k, a, i, nocc, nao)
    !        Sq(b,j) = Sq(b,j) + DW * Waibjck
    !                                                          
    !        Call buildDW(t2, DW, c, k, a, i, b, j, nocc, nao)
    !        Sq(c,k) = Sq(c,k) + DW * Waibjck
    !                                                          
    !        Call buildDW(t2, DW, a, i, c, k, b, j, nocc, nao)
    !        Sq(a,i) = Sq(a,i) + DW * Waibjck
    !                                                          
    !        Call buildDW(t2, DW, b, j, a, i, c, k, nocc, nao)
    !        Sq(b,j) = Sq(b,j) + DW * Waibjck
    !                                                          
    !        Call buildDW(t2, DW, c, k, b, j, a, i, nocc, nao)
    !        Sq(c,k) = Sq(c,k) + DW * Waibjck

    !    End Do
    !    End Do
    !    End Do
    !End Do
    !End Do
    !!$omp end do
    !!$omp end parallel
    !End Subroutine

    Subroutine StructFac_T3(eri, t2, moe, Sq, nocc, nao)
    Implicit None
    Integer,            Intent(in)  :: nocc, nao
    Real (kind=pr),     Intent(in)  :: moe(nao), eri(nao,nao)
    Real (kind=pr),     Intent(in)  :: t2(nocc+1:nao,nocc,nocc)
    Real (kind=pr),     Intent(out) :: Sq(nocc+1:nao, nocc)
    
    Integer                         :: a, b, c, i, j, k
    Integer                         :: m, f, p, I_idx, J_idx, K_idx
    Integer                         :: p_i(6), p_j(6), p_k(6)
    Real (kind=pr)                  :: denom, perm_weight
    Real (kind=pr)                  :: W1, W2, W3, W4, W5, W6, W_cur, DW_total
    Real (kind=pr)                  :: Waib(6)

    Sq = Zero

    !$omp parallel default(shared)
    !$omp do schedule(dynamic) reduction(+:Sq) &
    !$omp& private(b, c, i, j, k, m, f, denom, perm_weight, &
    !$omp& W1, W2, W3, W4, W5, W6, Waib, p_i, p_j, p_k, p, I_idx, J_idx, K_idx, W_cur, DW_total)
    Do a = nocc+1, nao
        Do b = nocc+1, a
            
            ! RESTRICTED LOOPS: 6x fewer iterations
            Do i = 1, nocc
            Do j = 1, i
            Do k = 1, j
                
                Call qconserv_c(i, j, k, a, b, c)
                If (c <= nocc) cycle
                If (c > b) cycle
                
                denom = moe(i) + moe(j) + moe(k) - moe(a) - moe(b) - moe(c)

                ! a, b, c symmetry handling (from your original Sq code)
                If (a == c) then
                    ! denom = denom
                Else If (a == b .or. b == c) then
                    denom = denom / 3.0_pr
                Else
                    denom = denom / 6.0_pr
                End If

                ! i, j, k permutation weight (prevents overcounting boundaries)
                If (i == k) then
                    perm_weight = 1.0_pr / 6.0_pr
                Else If (i == j .or. j == k) then
                    perm_weight = 0.5_pr
                Else
                    perm_weight = 1.0_pr
                End If

                ! Evaluate buildW exactly ONCE for the unique triplet
                Call buildW(eri, t2, W1, a, b, c, i, j, k, nocc, nao)
                Call buildW(eri, t2, W2, a, b, c, j, k, i, nocc, nao)
                Call buildW(eri, t2, W3, a, b, c, k, i, j, nocc, nao)
                Call buildW(eri, t2, W4, a, b, c, j, i, k, nocc, nao)
                Call buildW(eri, t2, W5, a, b, c, i, k, j, nocc, nao)
                Call buildW(eri, t2, W6, a, b, c, k, j, i, nocc, nao)
                
                ! Algebraically construct Waibjck for all 6 permutations
                Waib(1) = (4.0_pr * W1 + W2 + W3 - 2.0_pr * (W4 + W5 + W6)) / denom * perm_weight
                Waib(2) = (4.0_pr * W2 + W3 + W1 - 2.0_pr * (W6 + W4 + W5)) / denom * perm_weight
                Waib(3) = (4.0_pr * W3 + W1 + W2 - 2.0_pr * (W5 + W6 + W4)) / denom * perm_weight
                Waib(4) = (4.0_pr * W4 + W5 + W6 - 2.0_pr * (W1 + W2 + W3)) / denom * perm_weight
                Waib(5) = (4.0_pr * W5 + W6 + W4 - 2.0_pr * (W3 + W1 + W2)) / denom * perm_weight
                Waib(6) = (4.0_pr * W6 + W4 + W5 - 2.0_pr * (W2 + W3 + W1)) / denom * perm_weight

                p_i = (/ i, j, k, j, i, k /)
                p_j = (/ j, k, i, i, k, j /)
                p_k = (/ k, i, j, k, j, i /)

                ! Execute the memory updates for the 6 permutations mapped correctly
                Do p = 1, 6
                    I_idx = p_i(p)
                    J_idx = p_j(p)
                    K_idx = p_k(p)
                    W_cur = Waib(p)

                    ! --- Updates for Sq(a, I_idx) ---
                    DW_total = 0.0_pr
                    
                    Call qconserv2(I_idx, a, J_idx, m)
                    If (m <= nocc .and. m > 0) DW_total = DW_total - t2(b,m,K_idx)
                    Call qconserv2(a, I_idx, b, f)
                    If (f > nocc) DW_total = DW_total + t2(c,K_idx,J_idx)
                    
                    Call qconserv2(I_idx, a, K_idx, m)
                    If (m <= nocc .and. m > 0) DW_total = DW_total - t2(c,m,J_idx)
                    Call qconserv2(a, I_idx, c, f)
                    If (f > nocc) DW_total = DW_total + t2(b,J_idx,K_idx)
                    
                    Sq(a,I_idx) = Sq(a,I_idx) + DW_total * W_cur

                    ! --- Updates for Sq(b, J_idx) ---
                    DW_total = 0.0_pr
                    
                    Call qconserv2(J_idx, b, K_idx, m)
                    If (m <= nocc .and. m > 0) DW_total = DW_total - t2(c,m,I_idx)
                    Call qconserv2(b, J_idx, c, f)
                    If (f > nocc) DW_total = DW_total + t2(a,I_idx,K_idx)
                    
                    Call qconserv2(J_idx, b, I_idx, m)
                    If (m <= nocc .and. m > 0) DW_total = DW_total - t2(a,m,K_idx)
                    Call qconserv2(b, J_idx, a, f)
                    If (f > nocc) DW_total = DW_total + t2(c,K_idx,I_idx)
                    
                    Sq(b,J_idx) = Sq(b,J_idx) + DW_total * W_cur

                    ! --- Updates for Sq(c, K_idx) ---
                    DW_total = 0.0_pr
                    
                    Call qconserv2(K_idx, c, I_idx, m)
                    If (m <= nocc .and. m > 0) DW_total = DW_total - t2(a,m,J_idx)
                    Call qconserv2(c, K_idx, a, f)
                    If (f > nocc) DW_total = DW_total + t2(b,J_idx,I_idx)
                    
                    Call qconserv2(K_idx, c, J_idx, m)
                    If (m <= nocc .and. m > 0) DW_total = DW_total - t2(b,m,I_idx)
                    Call qconserv2(c, K_idx, b, f)
                    If (f > nocc) DW_total = DW_total + t2(a,I_idx,J_idx)
                    
                    Sq(c,K_idx) = Sq(c,K_idx) + DW_total * W_cur

                End Do
                
            End Do
            End Do
            End Do
            
        End Do
    End Do
    !$omp end do
    !$omp end parallel
    End Subroutine

    End Module
