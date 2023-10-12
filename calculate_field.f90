

SUBROUTINE Get_Field(nconfig, nmoltypes, nmols, natoms, which_is_wat, rmax, L, &
                    & rO, r1, r2, rmol, charges, dot1, dot2, eOH1, eOH2)
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: nconfig, nmoltypes, nmols(:), natoms(:), which_is_wat
    REAL, INTENT(IN) :: rmax, L(3)
    REAL, INTENT(IN) :: rO(:,:,:), r1(:,:,:), r2(:,:,:), rmol(:,:,:,:,:)
    REAL, INTENT(IN) :: charges(:,:)

    ! Output variables
    REAL, DIMENSION(:,:) :: dot1, dot2
    REAL, DIMENSION(:,:,:) :: eOH1, eOH2
    REAL, DIMENSION(3) :: ef1_tmp, ef2_tmp

    !REAL, DIMENSION(nmols(which_is_wat),3,nconfig) :: efield1, efield2

    INTEGER imol, k, z, p, type, jatom

    REAL :: norm1, norm2
    REAL :: dist1o, dist2o, dist1, dist2
    REAL, DIMENSION(3) :: rtmp1o, rtmp2o, rtmp1, rtmp2, ratom

    REAL :: angperau
    


    ! Convert Units
    angperau = 0.52917721092d0

    !efield1 = 0.0; efield2 = 0.0
    eOH1 = 0.0; eOH2 = 0.0
    !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE),  SCHEDULE(STATIC), &
    !$OMP SHARED( rO, r1, r2, rmol, charges), & 
    !$OMP SHARED(eOH1, eOH2), &
    !$OMP SHARED(dot1, dot2)
    DO z=1, nconfig
        ! Loop over the water molecules to get the electric field
        DO imol=1, nmols(which_is_wat)

            ! Get the OH vectors for the water molecule
            CALL OH_Vector(r1(imol,:,z),rO(imol,:,z),eOH1(imol,:,z))
            CALL OH_Vector(r2(imol,:,z),rO(imol,:,z),eOH2(imol,:,z))
            
            ef1_tmp = 0.0; ef2_tmp = 0.0
            ! Calculate the field contribution...
            DO type=1, nmoltypes
            ! ... for water
            IF (type == which_is_wat) THEN
                DO p=1, nmols(type)
                    IF (p .eq. imol) CYCLE ! Don't calculate if same mol
                    ! Check oxygen distances
                    CALL PBC_DIST(rO(p,:,z), r1(imol,:,z), L, dist1o, rtmp1o(:))
                    CALL PBC_DIST(rO(p,:,z), r2(imol,:,z), L, dist2o, rtmp2o(:))

                    ! Chekc if the O -> h(imol) distance is less than rmax
                    IF (dist1o .le. rmax) THEN
                        ! Add field contribution from O
                        DO k=1,3
                            ef1_tmp(k) = ef1_tmp(k) &
                            & + charges(type,1) * (r1(imol,k,z) - rtmp1o(k))/(dist1o**3)
                            !& + E_Cont(charges(type,1), r1(imol,k,z), rtmp1o(k), dist1o)
                        ENDDO

                        ! Add field contribution from H1
                        dist1 = 0.0
                        DO k=1, 3
                            rtmp1(k) = r1(p,k,z) - L(k)*anint((r1(p,k,z)-r1(imol,k,z))/L(k))
                            dist1 = dist1 + (r1(imol,k,z) - rtmp1(k))**2
                        ENDDO
                        dist1 = Sqrt(dist1)
                        DO k=1,3
                            ef1_tmp(k) = ef1_tmp(k) &
                            & + charges(type,2) * (r1(imol,k,z) - rtmp1(k))/(dist1**3)
                            !& + E_Cont(charges(type,2), r1(imol,k,z), rtmp1(k), dist1)
                        ENDDO
                        
                        ! Add field contribution from H2
                        dist1 = 0.0
                        DO k=1, 3
                            rtmp1(k) = r2(p,k,z) - L(k)*anint((r2(p,k,z)-r1(imol,k,z))/L(k))
                            dist1 = dist1 + (r1(imol,k,z) - rtmp1(k))**2
                        ENDDO
                        dist1 = Sqrt(dist1)
                        DO k=1,3
                            ef1_tmp(k) = ef1_tmp(k) &
                            & + charges(type,3) * (r1(imol,k,z) - rtmp1(k))/(dist1**3)
                            !& + E_Cont(charges(type,3), r1(imol,k,z), rtmp1(k), dist1)
                        ENDDO
                    ENDIF ! (dist1o .le. rmax)

                    IF (dist2o .le. rmax) THEN
                        ! Add field contribution from O
                        DO k=1,3
                            ef2_tmp(k) = ef2_tmp(k) &
                            & + charges(type,1) * (r2(imol,k,z) - rtmp2o(k))/(dist2o**3)
                            !& + E_Cont(charges(type,1), r2(imol,k,z), rtmp2o(k), dist2o)
                        ENDDO

                        ! Add field contribution from H1
                        dist2 = 0.0
                        DO k=1, 3
                            rtmp2(k) = r1(p,k,z) - L(k)*anint((r1(p,k,z)-r2(imol,k,z))/L(k))
                            dist2 = dist2 + (r2(imol,k,z) - rtmp2(k))**2
                        ENDDO
                        dist2 = Sqrt(dist2)
                        DO k=1,3
                            ef2_tmp(k) = ef2_tmp(k) &
                            & + charges(type,2) * (r2(imol,k,z) - rtmp2(k))/(dist2**3)
                            !& + E_Cont(charges(type,2), r2(imol,k,z), rtmp2(k), dist2)
                        ENDDO
                        
                        ! Add field contribution from H2
                        dist2 = 0.0
                        DO k=1, 3
                            rtmp2(k) = r2(p,k,z) - L(k)*anint((r2(p,k,z)-r2(imol,k,z))/L(k))
                            dist2 = dist2 + (r2(imol,k,z) - rtmp2(k))**2
                        ENDDO
                        dist2 = Sqrt(dist2)
                        DO k=1,3
                            ef2_tmp(k) = ef2_tmp(k) &
                            & + charges(type,3) * (r2(imol,k,z) - rtmp2(k))/(dist2**3)
                            !& + E_Cont(charges(type,3), r2(imol,k,z), rtmp2(k), dist2)
                        ENDDO
                    ENDIF ! (dist1o .le. rmax)
                ENDDO ! p
            ! ... for the rest
            ELSE
                DO p=1, nmols(type)
                DO jatom=1, natoms(type)
                    ratom = rmol(type, p, jatom, :, z)
                    ! Contribution of atom on H1 of imol
                    dist1 = 0.0
                    do k=1,3 
                        rtmp1(k) = ratom(k) - L(k)*anint((ratom(k)-r1(imol,k,z))/L(k))
                        dist1 = dist1 + (r1(imol,k,z) - rtmp1(k))**2
                    ENDDO ! k
                    dist1 = sqrt(dist1)
                    IF (dist1 .le. rmax) THEN
                        DO k=1,3
                            ef1_tmp(k) = ef1_tmp(k) &
                            & + charges(type,jatom) * (r1(imol,k,z) - rtmp1(k))/(dist1**3)
                            !& + E_Cont(charges(type,jatom), r1(imol,k,z), rtmp1(k), dist1)
                        ENDDO
                    ENDIF ! (dist1 .le. rmax)

                    ! Contribution of atom on H2 of imol
                    dist2 = 0.0
                    do k=1,3 
                        rtmp2(k) = ratom(k) - L(k)*anint((ratom(k)-r2(imol,k,z))/L(k))
                        dist2 = dist2 + (r2(imol,k,z) - rtmp2(k))**2
                    ENDDO ! k
                    dist2 = sqrt(dist2)
                    IF (dist2 .le. rmax) THEN
                        DO k=1,3
                            ef2_tmp(k) = ef2_tmp(k) &
                            & + charges(type,jatom) * (r2(imol,k,z) - rtmp2(k))/(dist2**3)
                            !& + E_Cont(charges(type,jatom), r2(imol,k,z), rtmp2(k), dist2)
                            
                        ENDDO
                    ENDIF ! (dist1 .le. rmax)
                ENDDO ! jatom
                ENDDO ! p
            ENDIF ! (type == which_is_wat)
            ENDDO ! type
            
            ! Convert units and project the 
            DO k=1,3
                ef1_tmp(k) = angperau**2*ef1_tmp(k)
                ef2_tmp(k) = angperau**2*ef2_tmp(k) 
            ENDDO
            IF (z == 1 .and. imol == 1) THEN
                WRITE(*,*) eOH1(1,1,1), ef1_tmp(1), "test"
            END IF
            dot1(imol,z) = Dot_Product(eOH1(imol,:,z), ef1_tmp(:))
            dot2(imol,z) = Dot_Product(eOH2(imol,:,z), ef2_tmp(:))
        ENDDO !imol
    ENDDO ! z
    !$OMP END PARALLEL DO
END SUBROUTINE

SUBROUTINE OH_Vector(ra, rb, eOH)

    IMPLICIT NONE

    REAL, DIMENSION(3), INTENT(IN) :: ra, rb
    REAL, DIMENSION(3), INTENT(OUT) :: eOH

    REAL :: norm
    INTEGER :: k

    eOH = 0.0
    DO k=1,3
        eOH(k) = ra(k) - rb(k)
        norm = norm + eOH(k)**2
    ENDDO ! k
    norm = SQRT(norm)
    DO k=1,3
        eOH(k) = eOH(k) / norm
    ENDDO ! k

END SUBROUTINE OH_Vector

SUBROUTINE PBC_Dist(ra, rb, L, dist, vector)

    IMPLICIT NONE

    REAL, DIMENSION(3), INTENT(IN) :: ra, rb, L
    REAL, INTENT(OUT) :: dist
    REAL, DIMENSION(3), INTENT(OUT) :: vector

    INTEGER :: k

    dist =0.0

    DO k=1,3
        vector(k) = ra(k) - L(k)*ANINT((ra(k)-rb(k))/L(k))
        dist = dist + (rb(k) - vector(k))**2
    ENDDO
    dist = SQRT(dist)

END SUBROUTINE PBC_DIST
