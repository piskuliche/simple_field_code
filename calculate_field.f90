

SUBROUTINE Get_Field(nconfig, nmoltypes, nmols, natoms, which_is_wat, rmax, L, &
                    & rO, r1, r2, rmol, charges, dot1, dot2, eOH1, eOH2)
    use field_module
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: nconfig, nmoltypes, nmols(:), natoms(:), which_is_wat
    REAL, INTENT(IN) :: rmax, L(3)
    REAL, INTENT(IN) :: rO(:,:,:), r1(:,:,:), r2(:,:,:), rmol(:,:,:,:,:)
    REAL, INTENT(IN) :: charges(:,:)

    ! Output variables
    REAL, DIMENSION(:,:) :: dot1, dot2
    REAL, DIMENSION(:,:,:) :: eOH1, eOH2

    REAL, DIMENSION(nmols(which_is_wat),3,nconfig) :: efield1, efield2

    INTEGER imol, k, z, p, type, jatom

    REAL :: norm1, norm2
    REAL :: dist1o, dist2o, dist1, dist2
    REAL, DIMENSION(3) :: rtmp1o, rtmp2o, rtmp1, rtmp2, ratom

    REAL :: angperau


    angperau = 0.52917721092d0

    efield1 = 0.0; efield2 = 0.0
    DO z=1, nconfig
        ! Loop over the water molecules to get the electric field
        DO imol=1, nmols(which_is_wat)
            ! Get OH vector
            norm1 = 0.0; norm2 = 0.0
            DO k=1, 3
                eOH2(imol,k,z) = r2(imol,k,z) - rO(imol,k,z)
                eOH1(imol,k,z) = r1(imol,k,z) - rO(imol,k,z)
                norm1 = norm1 + eOH1(imol,k,z)**2
                norm2 = norm2 + eOH2(imol,k,z)**2
            ENDDO ! k
            norm1 = SQRT(norm1); norm2 = SQRT(norm2)

            ! Normalize the vectors
            DO k=1,3
                eOH1(imol,k,z) = eOH1(imol,k,z) / norm1
                eOH2(imol,k,z) = eOH2(imol,k,z) / norm2
            ENDDO ! k

            ! Calculate the field contribution...
            DO type=1, nmoltypes
            ! ... for water
            IF (type == which_is_wat) THEN
                DO p=1, nmols(type)
                    IF (p .eq. imol) CYCLE ! Don't calculate if same mol
                    ! Check oxygen distances
                    dist1o = 0.0; dist2o = 0.0
                    DO k=1, 3
                        rtmp1o(k) = rO(p,k,z) - L(k)*anint((rO(p,k,z)-r1(imol,k,z))/L(k))
                        rtmp2o(k) = rO(p,k,z) - L(k)*anint((rO(p,k,z)-r2(imol,k,z))/L(k))
                        dist1o = dist1o + (r1(imol,k,z) - rtmp1o(k))**2
                        dist2o = dist2o + (r2(imol,k,z) - rtmp2o(k))**2
                    ENDDO
                    dist1o = Sqrt(dist1o); dist2o = Sqrt(dist2o)

                    ! Chekc if the O -> h(imol) distance is less than rmax
                    IF (dist1o .le. rmax) THEN
                        ! Add field contribution from O
                        DO k=1,3
                            efield1(imol,k,z) = efield1(imol,k,z) &
                            & + E_Cont(charges(type,1), r1(imol,k,z), rtmp1o(k), dist1o)
                        ENDDO

                        ! Add field contribution from H1
                        dist1 = 0.0
                        DO k=1, 3
                            rtmp1(k) = r1(p,k,z) - L(k)*anint((r1(p,k,z)-r1(imol,k,z))/L(k))
                            dist1 = dist1 + (r1(imol,k,z) - rtmp1(k))**2
                        ENDDO
                        dist1 = Sqrt(dist1)
                        DO k=1,3
                            efield1(imol,k,z) = efield1(imol,k,z) &
                            & + E_Cont(charges(type,2), r1(imol,k,z), rtmp1(k), dist1)
                        ENDDO
                        
                        ! Add field contribution from H2
                        dist1 = 0.0
                        DO k=1, 3
                            rtmp1(k) = r2(p,k,z) - L(k)*anint((r2(p,k,z)-r1(imol,k,z))/L(k))
                            dist1 = dist1 + (r2(imol,k,z) - rtmp1(k))**2
                        ENDDO
                        dist1 = Sqrt(dist1)
                        DO k=1,3
                            efield1(imol,k,z) = efield1(imol,k,z) &
                            & + E_Cont(charges(type,3), r1(imol,k,z), rtmp1(k), dist1)
                        ENDDO
                    ENDIF ! (dist1o .le. rmax)

                    IF (dist2o .le. rmax) THEN
                        ! Add field contribution from O
                        DO k=1,3
                            efield2(imol,k,z) = efield2(imol,k,z) &
                            & + E_Cont(charges(type,1), r2(imol,k,z), rtmp2o(k), dist2o)
                        ENDDO

                        ! Add field contribution from H1
                        dist2 = 0.0
                        DO k=1, 3
                            rtmp2(k) = r1(p,k,z) - L(k)*anint((r1(p,k,z)-r2(imol,k,z))/L(k))
                            dist2 = dist2 + (r2(imol,k,z) - rtmp2(k))**2
                        ENDDO
                        dist2 = Sqrt(dist2)
                        DO k=1,3
                            efield2(imol,k,z) = efield2(imol,k,z) &
                            & + E_Cont(charges(type,2), r2(imol,k,z), rtmp2(k), dist2)
                        ENDDO
                        
                        ! Add field contribution from H2
                        dist2 = 0.0
                        DO k=1, 3
                            rtmp2(k) = r2(p,k,z) - L(k)*anint((r2(p,k,z)-r2(imol,k,z))/L(k))
                            dist2 = dist2 + (r2(imol,k,z) - rtmp2(k))**2
                        ENDDO
                        dist2 = Sqrt(dist2)
                        DO k=1,3
                            efield2(imol,k,z) = efield2(imol,k,z) &
                            & + E_Cont(charges(type,3), r2(imol,k,z), rtmp2(k), dist2)
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
                            efield1(imol,k,z) = efield1(imol,k,z) &
                            & + E_Cont(charges(type,jatom), r1(imol,k,z), rtmp1(k), dist1)
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
                            efield2(imol,k,z) = efield2(imol,k,z) &
                            & + E_Cont(charges(type,jatom), r2(imol,k,z), rtmp2(k), dist2)
                        ENDDO
                    ENDIF ! (dist1 .le. rmax)
                ENDDO ! jatom
                ENDDO ! p
            ENDIF ! (type == which_is_wat)
            ENDDO ! type
            
            ! Convert units and project the 
            DO k=1,3
                efield1(imol,k,z) = angperau**2*efield1(imol,k,z)
                efield2(imol,k,z) = angperau**2*efield2(imol,k,z)
            ENDDO

            dot1(imol,z) = Dot_Product(eOH1(imol,:,z), efield1(imol,:,z))
            dot2(imol,z) = Dot_Product(eOH2(imol,:,z), efield2(imol,:,z))
        ENDDO !imol
    ENDDO ! z
END SUBROUTINE
