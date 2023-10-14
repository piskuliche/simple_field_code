

SUBROUTINE Get_Field(nconfig, nmoltypes, nmols, natoms, which_is_wat, rmax, L, &
                    & rO, r1, r2, rmol, charges, dot1, dot2, eOH1, eOH2)
    USE ieee_arithmetic, ONLY: ieee_is_finite
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: maxconfig = 100
    INTEGER :: num_chunks
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

    INTEGER :: imol, k, z, p, type, jatom, chunk
    INTEGER :: iconfig

    REAL :: norm1, norm2
    REAL :: dist1o, dist2o, dist1, dist2
    REAL, DIMENSION(3) :: rtmp1o, rtmp2o, rtmp1, rtmp2, ratom

    REAL :: angperau
    
    WRITE(*,*) "Get field without samples"

    ! Convert Units
    angperau = 0.52917721092d0

    num_chunks = ceiling(real(nconfig)/ real(maxconfig))

    DO chunk=1, num_chunks
        DO z=1, maxconfig
            IF (chunk*maxconfig + z > nconfig) then
                EXIT
            ENDIF
            CALL Read_XYZ_Frame(12, nmoltypes, nmols, natoms, which_is_wat, L, rO(:,:,z), r1(:,:,z), r2(:,:,z), rmol(:,:,:,:,z))
        ENDDO ! z 

        eOH1 = 0.0; eOH2 = 0.0
        !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE),  SCHEDULE(STATIC), &
        !$OMP SHARED( rO, r1, r2, rmol, charges), & 
        !$OMP SHARED(eOH1, eOH2), &
        !$OMP SHARED(dot1, dot2)
        DO z=1, maxconfig
            iconfig = chunk*maxconfig + z
            ! Loop over the water molecules to get the electric field
            DO imol=1, nmols(which_is_wat)
                ! Zero the eoh
                ef1_tmp = 0.0; ef2_tmp = 0.0

                ! Get the OH vectors for the water molecule
                CALL OH_Vector(r1(imol,:,z),rO(imol,:,z),eOH1(imol,:,iconfig))
                CALL OH_Vector(r2(imol,:,z),rO(imol,:,z),eOH2(imol,:,iconfig))

                ! Calculate the field contribution...
                DO type=1, nmoltypes
                ! ... for water
                IF (type == which_is_wat) THEN
                    DO p=1, nmols(type)
                        IF (p .eq. imol) CYCLE ! Don't calculate if same mol
                        ! Check oxygen distances
                        CALL PBC_Dist(rO(p,:,z), r1(imol,:,z), L, dist1o, rtmp1o(:))
                        CALL PBC_Dist(rO(p,:,z), r2(imol,:,z), L, dist2o, rtmp2o(:))

                        ! Contributions to H1
                        IF (dist1o .le. rmax) THEN

                            ! Add field contribution from O
                            CALL Field_Contribution(charges(type,1), r1(imol,:,z), rtmp1o(:), dist1o, ef1_tmp(:))

                            ! Add field contribution from H1
                            CALL PBC_Dist(r1(p,:,z), r1(imol,:,z), L, dist1, rtmp1(:))
                            CALL Field_Contribution(charges(type,2), r1(imol,:,z), rtmp1(:), dist1, ef1_tmp(:))

                            ! Add field contribution from H2
                            CALL PBC_Dist(r2(p,:,z), r1(imol,:,z), L, dist1, rtmp1(:))
                            CALL Field_Contribution(charges(type,3), r1(imol,:,z), rtmp1(:), dist1, ef1_tmp(:))
                            
                        ENDIF ! (dist1o .le. rmax)

                        ! Contributions to H2
                        IF (dist2o .le. rmax) THEN

                            ! Add field contribution from O
                            CALL Field_Contribution(charges(type,1), r2(imol,:,z), rtmp2o(:), dist2o, ef2_tmp(:))
                            ! Add field contribution from H1
                            CALL PBC_Dist(r1(p,:,z), r2(imol,:,z), L, dist2, rtmp2(:))
                            CALL Field_Contribution(charges(type,2), r2(imol,:,z), rtmp2(:), dist2, ef2_tmp(:))

                            ! Add field contribution from H2
                            CALL PBC_Dist(r2(p,:,z), r2(imol,:,z), L, dist2, rtmp2(:))
                            CALL Field_Contribution(charges(type,3), r2(imol,:,z), rtmp2(:), dist2, ef2_tmp(:))

                            IF (.NOT. ieee_is_finite(ef2_tmp(1))) THEN
                                WRITE(*,*) "ef2_tmp not finite"
                                WRITE(*,*) dist2o, "r2r2", ef2_tmp(1)
                                WRITE(*,*) L, dist2
                                WRITE(*,*) ef2_tmp(:)
                                ERROR STOP "Error: Exiting"
                            ENDIF

                        ENDIF ! (dist1o .le. rmax)
                    ENDDO ! p
                ! ... for the rest
                ELSE
                    DO p=1, nmols(type)
                        DO jatom=1, natoms(type)
                            ratom = rmol(type, p, jatom, :, z)
                            ! Contribution from jatom to H1
                            CALL PBC_Dist(ratom(:), r1(imol,:,z), L, dist1, rtmp1(:))
                            IF (dist1 .le. rmax) THEN
                                CALL Field_Contribution(charges(type,jatom), r1(imol,:,z), rtmp1(:), dist1, ef1_tmp(:))
                            ENDIF

                            ! Contribtuion from jatom to H2
                            CALL PBC_Dist(ratom(:), r2(imol,:,z), L, dist2, rtmp2(:))
                            IF (dist2 .le. rmax) THEN
                                CALL Field_Contribution(charges(type,jatom), r2(imol,:,z), rtmp2(:), dist2, ef2_tmp(:))
                            ENDIF
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
                    WRITE(*,*) eOH2(1,1,1), ef2_tmp(1), "test"
                END IF
                dot1(imol,iconfig) = Dot_Product(eOH1(imol,:,iconfig), ef1_tmp(:))
                dot2(imol,iconfig) = Dot_Product(eOH2(imol,:,iconfig), ef2_tmp(:))
                IF (z == 1 .and. imol == 1) THEN
                    WRITE(*,*) dot1(1,1), dot2(1,1), "testdot"
                END IF
            ENDDO !imol
        ENDDO ! z
        !$OMP END PARALLEL DO
    ENDDO ! Chunks
END SUBROUTINE Get_Field

SUBROUTINE Get_Field_Samples(nconfig, nmoltypes, nmols, natoms, which_is_wat, rmax, L, samples, &
                    & rO, r1, r2, rmol, charges, dot1, dot2, eOH1, eOH2)
    IMPLICIT NONE

    INTEGER, PARAMETER :: maxconfig = 100
    INTEGER :: num_chunks 

    ! Input variables
    INTEGER, INTENT(IN) :: nconfig, nmoltypes, nmols(:), natoms(:), which_is_wat
    REAL, INTENT(IN) :: rmax, L(3)
    REAL, INTENT(IN) :: rO(:,:,:), r1(:,:,:), r2(:,:,:), rmol(:,:,:,:,:)
    REAL, INTENT(IN) :: charges(:,:)
    INTEGER, INTENT(IN) :: samples(:)
    

    ! Output variables

    REAL, DIMENSION(:,:) :: dot1, dot2
    REAL, DIMENSION(:,:,:) :: eOH1, eOH2
    REAL, DIMENSION(3) :: ef1_tmp, ef2_tmp

    !REAL, DIMENSION(nmols(which_is_wat),3,nconfig) :: efield1, efield2

    INTEGER :: imol, ival, k, z, p, type, jatom, chunk
    INTEGER :: iconfig

    REAL :: norm1, norm2
    REAL :: dist1o, dist2o, dist1, dist2
    REAL, DIMENSION(3) :: rtmp1o, rtmp2o, rtmp1, rtmp2, ratom

    REAL :: angperau
    INTEGER :: nsamples
    
    ! Get the number of elements in samples
    nsamples = SIZE(samples)

    ! Convert Units
    angperau = 0.52917721092d0

    num_chunks = ceiling(real(nconfig)/ real(maxconfig))

    DO chunk=1, num_chunks
        DO z=1, maxconfig
            IF (chunk*maxconfig + z > nconfig) then
                EXIT
            ENDIF
            CALL Read_XYZ_Frame(12, nmoltypes, nmols, natoms, which_is_wat, L, rO(:,:,z), r1(:,:,z), r2(:,:,z), rmol(:,:,:,:,z))
        ENDDO ! z 

        !efield1 = 0.0; efield2 = 0.0
        eOH1 = 0.0; eOH2 = 0.0
        !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE),  SCHEDULE(STATIC), &
        !$OMP SHARED( rO, r1, r2, rmol, charges), & 
        !$OMP SHARED(eOH1, eOH2), &
        !$OMP SHARED(dot1, dot2)
        DO z=1, nconfig
            iconfig = chunk*maxconfig + z
            ! Loop over the water molecules to get the electric field
            DO ival=1, nsamples
                imol = samples(ival)

                ! Zero the eoh
                ef1_tmp = 0.0; ef2_tmp = 0.0

                ! Get the OH vectors for the water molecule
                CALL OH_Vector(r1(imol,:,z), rO(imol,:,z), eOH1(imol,:,iconfig))
                CALL OH_Vector(r2(imol,:,z), rO(imol,:,z), eOH2(imol,:,iconfig))

                IF (z == 1 .and. ival == 1) THEN
                    WRITE(*,*) ival, imol
                    WRITE(*,*) r1(imol,:,z), rO(imol,:,z)
                    WRITE(*,*) eOH1(1,1,1), eOH2(1,1,1), "test"
                END IF

                ! Calculate the field contribution...
                DO type=1, nmoltypes
                ! ... for water
                IF (type == which_is_wat) THEN
                    DO p=1, nmols(type)
                        IF (p .eq. imol) CYCLE ! Don't calculate if same mol
                        ! Check oxygen distances
                        CALL PBC_Dist(rO(p,:,z), r1(imol,:,z), L, dist1o, rtmp1o(:))
                        CALL PBC_Dist(rO(p,:,z), r2(imol,:,z), L, dist2o, rtmp2o(:))

                        ! Contributions to H1
                        IF (dist1o .le. rmax) THEN

                            ! Add field contribution from O
                            CALL Field_Contribution(charges(type,1), r1(imol,:,z), rtmp1o(:), dist1o, ef1_tmp(:))

                            ! Add field contribution from H1
                            CALL PBC_Dist(r1(p,:,z), r1(imol,:,z), L, dist1, rtmp1(:))
                            CALL Field_Contribution(charges(type,2), r1(imol,:,z), rtmp1(:), dist1, ef1_tmp(:))

                            ! Add field contribution from H2
                            CALL PBC_Dist(r2(p,:,z), r1(imol,:,z), L, dist1, rtmp1(:))
                            CALL Field_Contribution(charges(type,3), r1(imol,:,z), rtmp1(:), dist1, ef1_tmp(:))
                            
                        ENDIF ! (dist1o .le. rmax)

                        ! Contributions to H2
                        IF (dist2o .le. rmax) THEN

                            ! Add field contribution from O
                            CALL Field_Contribution(charges(type,1), r2(imol,:,z), rtmp2o(:), dist2o, ef2_tmp(:))

                            ! Add field contribution from H1
                            CALL PBC_Dist(r1(p,:,z), r2(imol,:,z), L, dist2, rtmp2(:))
                            CALL Field_Contribution(charges(type,2), r2(imol,:,z), rtmp2(:), dist2, ef2_tmp(:))

                            ! Add field contribution from H2
                            CALL PBC_Dist(r2(p,:,z), r2(imol,:,z), L, dist2, rtmp2(:))
                            CALL Field_Contribution(charges(type,3), r2(imol,:,z), rtmp2(:), dist2, ef2_tmp(:))

                        ENDIF ! (dist1o .le. rmax)
                    ENDDO ! p
                ! ... for the rest
                ELSE
                    DO p=1, nmols(type)
                        DO jatom=1, natoms(type)
                            ratom = rmol(type, p, jatom, :, z)
                            ! Contribution from jatom to H1
                            CALL PBC_Dist(ratom(:), r1(imol,:,z), L, dist1, rtmp1(:))
                            IF (dist1 .le. rmax) THEN
                                CALL Field_Contribution(charges(type,jatom), r1(imol,:,z), rtmp1(:), dist1, ef1_tmp(:))
                            ENDIF

                            ! Contribtuion from jatom to H2
                            CALL PBC_Dist(ratom(:), r2(imol,:,z), L, dist2, rtmp2(:))
                            IF (dist2 .le. rmax) THEN
                                CALL Field_Contribution(charges(type,jatom), r2(imol,:,z), rtmp2(:), dist2, ef2_tmp(:))
                            ENDIF
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
                    WRITE(*,*) eOH2(1,1,1), ef2_tmp(1), "test"
                END IF
                dot1(imol,iconfig) = Dot_Product(eOH1(imol,:,iconfig), ef1_tmp(:))
                dot2(imol,iconfig) = Dot_Product(eOH2(imol,:,iconfig), ef2_tmp(:))
                IF (z == 1 .and. imol == 1) THEN
                    WRITE(*,*) dot1(1,1), dot2(1,1), "testdot"
                END IF
            ENDDO !imol
        ENDDO ! z
        !$OMP END PARALLEL DO
    ENDDO ! Chunks
END SUBROUTINE Get_Field_Samples



SUBROUTINE OH_Vector(ra, rb, eOH)

    IMPLICIT NONE

    REAL, DIMENSION(3), INTENT(IN) :: ra, rb
    REAL, DIMENSION(3), INTENT(OUT) :: eOH

    REAL :: norm
    INTEGER :: k

    eOH = 0.0
    norm = 0.0
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

END SUBROUTINE PBC_Dist

SUBROUTINE Field_Contribution(q, ra, vector, dist, efield)
    
    IMPLICIT NONE

    REAL, INTENT(IN) :: q, dist
    REAL, DIMENSION(3), INTENT(IN) :: ra, vector
    REAL, DIMENSION(3), INTENT(INOUT) :: efield

    INTEGER :: k
    DO k=1,3
        efield(k) = efield(k) + q * (ra(k) - vector(k))/(dist**3)
    ENDDO

END SUBROUTINE

SUBROUTINE All_Distances(na, nb, ra, rb, L, dist, dr_vec)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: na, nb
    REAL, DIMENSION(na, 3), INTENT(IN) :: ra
    REAL, DIMENSION(nb, 3), INTENT(IN) :: rb
    REAL, DIMENSION(3), INTENT(IN) :: L
    REAL, DIMENSION(na, nb) :: dist
    REAL, DIMENSION(na, nb, 3), INTENT(OUT) :: dr_vec

    INTEGER :: i, j, k

    DO i=1, na
        DO j=1,nb
            DO k=1,3
                dr_vec(i,j,k) = ra(i,k) - rb(j,k) - L(k)*ANINT((ra(i,k)-rb(j,k))/L(k))
            ENDDO
            dist(i,j) = SQRT(DOT_PRODUCT(dr_vec(i,j,:), dr_vec(i,j,:)))
        ENDDO
    ENDDO
    
END SUBROUTINE All_Distances
