!
!    Calculate field lines of structure
!

PROGRAM FieldLines
    IMPLICIT NONE
    INTEGER, PARAMETER:: r=SELECTED_REAL_KIND(P=15)
    REAL, PARAMETER:: pi=4*atan(1.)

    !Cell Parameters 
    REAL (KIND=r)::  cell_a, cell_b, cell_c
    REAL (KIND=r)::  cell_alpha, cell_beta, cell_gamma
    REAL (KIND=r)::  sin_alpha, sin_beta, sin_gamma
    REAL (KIND=r)::  cos_alpha, cos_beta, cos_gamma
    REAL (KIND=r)::  SIGMA
    
    
    INTEGER:: i, j, k, l, n_atoms, n_points, flag, id
    REAL (KIND=r), ALLOCATABLE:: u_atom(:), v_atom(:), w_atom(:), q_stru(:)
    REAL (KIND=r), ALLOCATABLE:: x_atom(:), y_atom(:), z_atom(:), q_atom(:)
    REAL (KIND=r):: E_x, E_y, E_z
    REAL (KIND=r):: E_total(3), E_aux(3), pos(3), dir(3), dr, distance
    REAL (KIND=r):: aux_x, aux_y, aux_z       
    
    CHARACTER (LEN=9):: Coordinates

    CHARACTER (LEN=6):: atom
    CHARACTER (LEN=23):: endlines, lines
    REAL (KIND=r):: cx, cy, cz

    REAL (KIND=r):: Cordinated_Radius, tol_E, m_E_total

!STRUCTURE PARAMETERS
    READ*, cell_a
    READ*, cell_b
    READ*, cell_c
    READ*, cell_alpha
    READ*, cell_beta
    READ*, cell_gamma

!INPUT PARAMETERS
    READ*, Coordinates    
    READ*, pos(1)
    READ*, pos(2)
    READ*, pos(3)
    READ*, n_atoms
    READ*, dr 
    READ*, n_points
  
    sin_alpha=sin(cell_alpha*pi/180)
    sin_beta=sin(cell_beta*pi/180)
    sin_gamma=sin(cell_gamma*pi/180)

    cos_alpha=cos(cell_alpha*pi/180)
    cos_beta=cos(cell_beta*pi/180)
    cos_gamma=cos(cell_gamma*pi/180)

    SIGMA=cell_a*cell_b*cell_c*SQRT(1-cos_alpha*cos_alpha-cos_beta*cos_beta-cos_gamma*cos_gamma+2*cos_alpha*cos_beta*cos_gamma)

    IF (Coordinates.eq."Crystallo") THEN
            pos(1)=cell_a*(pos(1)+2)+cell_b*cos_gamma*(pos(2)+2)+cell_c*cos_beta*(pos(3)+2)                             !lo convierto cartersianas de la 3 3 3
            pos(2)=cell_b*sin_gamma*(pos(2)+2)+cell_c*((cos_alpha-cos_beta*cos_gamma)/(sin_gamma))*(pos(3)+2)
            pos(3)=SIGMA/(cell_a*cell_b*sin_gamma)*(pos(3)+2)
    ELSE
     cx=1/cell_a*pos(1)-(cos_gamma)/(cell_a*sin_gamma)*pos(2)+cell_b*cell_c*(cos_alpha*cos_gamma-cos_beta)/(SIGMA*sin_gamma)*pos(3)           !lo paso a cristalografico
            cy=1/(cell_b*sin_gamma)*pos(2)+cell_a*cell_c*(cos_beta*cos_gamma-cos_alpha)/(SIGMA*sin_alpha)*pos(3)
            cz=cell_a*cell_b*sin_gamma/SIGMA*pos(3)

            pos(1)=cell_a*(cx+2)+cell_b*cos_gamma*(cy+2)+cell_c*cos_beta*(cz+2)                 !lo convierto a cartesianas en la 3 3 3
            pos(2)=cell_b*sin_gamma*(cy+2)+cell_c*((cos_alpha-cos_beta*cos_gamma)/(sin_gamma))*(cz+2)
            pos(3)=SIGMA/(cell_a*cell_b*sin_gamma)*(cz+2)
    END IF

    aux_x=pos(1)       ! en coordenadas cartesianas en la 3 3 3
    aux_y=pos(2)
    aux_z=pos(3)
    
    tol_E=0.1 
    ALLOCATE (u_atom(1:n_atoms), v_atom(1:n_atoms), w_atom(1:n_atoms), q_stru(1:n_atoms))       !cif
    ALLOCATE (x_atom(1:125*n_atoms), y_atom(1:125*n_atoms), z_atom(1:125*n_atoms), q_atom(1:125*n_atoms))       !replica 5x5x5
    
    OPEN(2,file="field_line.pdb",STATUS='NEW', ACTION='WRITE')


    atom="ATOM"
    id=416
    lines=" F   MOL       "
    endlines=" 1.00  0.00           F"

    
    CALL readatoms
    CALL replicaCelda

    DO l=1, n_points    
        flag=0
        cx=pos(1)
        cy=pos(2)
        cz=pos(3)
        
        CALL traslada
        
        WRITE(2,'(A6,I5,A18,F9.3,X,F7.3,F7.2,A21)')atom,id,lines,cx,cy,cz,endlines
        Cordinated_Radius=50.
        E_total=0.
        E_aux=0.
        DO WHILE (flag.eq.0)
            DO i=1, 125*n_atoms
        
                distance=(sqrt((pos(1)-x_atom(i))**2+(pos(2)-y_atom(i))**2+(pos(3)-z_atom(i))**2))
!                print*, distance
                IF (distance.le.Cordinated_Radius) THEN       
                        E_x=q_atom(i)/(distance**3)*(pos(1)-x_atom(i))
                        E_y=q_atom(i)/(distance**3)*(pos(2)-y_atom(i))
                        E_z=q_atom(i)/(distance**3)*(pos(3)-z_atom(i))
        
                        E_total(1)=E_total(1)+E_x
                        E_total(2)=E_total(2)+E_y    
                        E_total(3)=E_total(3)+E_z
!                        print*, E_x, E_y, E_z, q_atom(i), i, Cordinated_Radius
                END IF
            END DO
            m_E_total=SQRT((E_total(1))**2+(E_total(2))**2+(E_total(3))**2)

!            print*,m_E_total,  SQRT((E_total(1)-E_aux(1))**2+(E_total(2)-E_aux(2))**2+(E_total(3)-E_aux(3))**2)/m_E_total
            IF (SQRT((E_total(1)-E_aux(1))**2+(E_total(2)-E_aux(2))**2+(E_total(3)-E_aux(3))**2)/m_E_total.lt.tol_E) THEN
                    
                    dir(1)=E_total(1)/sqrt((E_total(1)**2+E_total(2)**2+E_total(3)**2))
                    dir(2)=E_total(2)/sqrt((E_total(1)**2+E_total(2)**2+E_total(3)**2))
                    dir(3)=E_total(3)/sqrt((E_total(1)**2+E_total(2)**2+E_total(3)**2))
        
                    pos(1)=pos(1)+dr*dir(1)
                    pos(2)=pos(2)+dr*dir(2)
                    pos(3)=pos(3)+dr*dir(3)
                    flag=1
            ELSE
                    Cordinated_Radius=Cordinated_radius+5
                    E_aux=E_total
            END IF

        END DO          !do while             
    END DO              !do l

    pos(1)=aux_x
    pos(2)=aux_y
    pos(3)=aux_z
    dr=-dr

        DO l=1, n_points    
        flag=0
        cx=pos(1)
        cy=pos(2)
        cz=pos(3)
        
        CALL traslada

        WRITE(2,'(A6,I5,A18,F9.3,X,F7.3,F7.2,A21)')atom,id,lines,cx,cy,cz,endlines
        
        Cordinated_Radius=50.
        E_total=0.
        E_aux=0.
        DO WHILE (flag.eq.0)
            DO i=1, 125*n_atoms
        
                distance=(sqrt((pos(1)-x_atom(i))**2+(pos(2)-y_atom(i))**2+(pos(3)-z_atom(i))**2))
                IF (distance.le.Cordinated_Radius) THEN       
                        E_x=q_atom(i)/(distance**3)*(pos(1)-x_atom(i))
                        E_y=q_atom(i)/(distance**3)*(pos(2)-y_atom(i))
                        E_z=q_atom(i)/(distance**3)*(pos(3)-z_atom(i))
        
                        E_total(1)=E_total(1)+E_x
                        E_total(2)=E_total(2)+E_y    
                        E_total(3)=E_total(3)+E_z
              
               END IF
            END DO
            m_E_total=SQRT((E_total(1))**2+(E_total(2))**2+(E_total(3))**2)

            IF (SQRT((E_total(1)-E_aux(1))**2+(E_total(2)-E_aux(2))**2+(E_total(3)-E_aux(3))**2)/m_E_total.lt.tol_E) THEN
                    
                    dir(1)=E_total(1)/sqrt((E_total(1)**2+E_total(2)**2+E_total(3)**2))
                    dir(2)=E_total(2)/sqrt((E_total(1)**2+E_total(2)**2+E_total(3)**2))
                    dir(3)=E_total(3)/sqrt((E_total(1)**2+E_total(2)**2+E_total(3)**2))
        
                    pos(1)=pos(1)+dr*dir(1)
                    pos(2)=pos(2)+dr*dir(2)
                    pos(3)=pos(3)+dr*dir(3)
                    flag=1
            ELSE
                    Cordinated_Radius=Cordinated_radius+5
                    E_aux=E_total
            END IF

        END DO          !do while             
    END DO              !do l

    CONTAINS
      
    SUBROUTINE readatoms
        
        OPEN(1,FILE="atom-positions.dat",STATUS='OLD', ACTION='READ')
        DO i=1, n_atoms    !recorro los atomos
            READ(1,*) u_atom(i), v_atom(i), w_atom(i), q_stru(i)
        END DO
        CLOSE (UNIT=1) 
    END SUBROUTINE readatoms

    SUBROUTINE replicaCelda
    IMPLICIT NONE
        INTEGER:: ll
        OPEN(25,file="UnitCell.pdb",STATUS='NEW', ACTION='WRITE')
        ll=0
        DO i=0, 4 
                DO j=0, 4 
                        DO k=0, 4 
                                DO l=1, n_atoms
x_atom(ll*n_atoms+l)=cell_a*(u_atom(l)+i)+cell_b*cos_gamma*(v_atom(l)+j)+cell_c*cos_beta*(w_atom(l)+k)
y_atom(ll*n_atoms+l)=cell_b*sin_gamma*(j+v_atom(l))+cell_c*((cos_alpha-cos_beta*cos_gamma)/(sin_gamma))*(w_atom(l)+k)
z_atom(ll*n_atoms+l)=SIGMA/(cell_a*cell_b*sin_gamma)*(k+w_atom(l))
q_atom(ll*n_atoms+l)=q_stru(l)
IF (ll.eq.0) THEN
WRITE(25,'(A6,I5,A18,F9.3,X,F7.3,F7.2,A21)')atom,id,lines,x_atom(l),y_atom(l),z_atom(l),endlines               !la replica en cartesianas
END IF
!print*, x_atom(ll+l),y_atom(ll+l),z_atom(ll+l),q_atom(ll+l), ll*n_atoms+l
                                END DO  !l
                        ll=ll+1
!                        print*, ll, ll*n_atoms
                        END DO  !k
                END DO  !j
        END DO  !i
    CLOSE(UNIT=25)    

    END SUBROUTINE replicaCelda


    SUBROUTINE traslada
        REAL (KIND=r)::  aux_u, aux_v, aux_w
        !lo paso a cristalografico
        aux_u=1/cell_a*cx-(cos_gamma)/(cell_a*sin_gamma)*cy+cell_b*cell_c*(cos_alpha*cos_gamma-cos_beta)/(SIGMA*sin_gamma)*cz
        aux_v=1/(cell_b*sin_gamma)*cy+cell_a*cell_c*(cos_beta*cos_gamma-cos_alpha)/(SIGMA*sin_alpha)*cz
        aux_w=cell_a*cell_b*sin_gamma/SIGMA*cz
       
        aux_u=aux_u-2
        aux_v=aux_v-2
        aux_w=aux_w-2

        cx=cell_a*aux_u+cell_b*cos_gamma*aux_v+cell_c*cos_beta*aux_w
        cy=cell_b*sin_gamma*aux_v+cell_c*((cos_alpha-cos_beta*cos_gamma)/(sin_gamma))*aux_w
        cz=SIGMA/(cell_a*cell_b*sin_gamma)*aux_w

    END SUBROUTINE traslada


        
END PROGRAM FieldLines
