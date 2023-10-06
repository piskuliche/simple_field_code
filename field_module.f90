MODULE field_module
    IMPLICIT NONE
    CONTAINS

    REAL FUNCTION E_Cont(q, r1, r2, dist) RESULT(contribution)
        IMPLICIT NONE
        REAL :: q, r1, r2, dist
        contribution = q * (r1 - r2) / dist**3
    END FUNCTION E_Cont

    include read_field_input.f90
    include calculate_field.f90

END MODULE