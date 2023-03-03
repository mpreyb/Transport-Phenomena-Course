program estanque
     implicit none
     integer,parameter :: M=21
     real : : dx , dy , Lx , Ly , dt , mu, nu , rho
     integer : : nit , nt , i , j , k , t , nx , ny
     real, dimension(M,M) : : ux , v , U, P, Pn , un , vn

     ! --------------------------------------------PARAMETROS DEL SISTEMA-----------------------------------------------------------
     Lx=2
     Ly=2                                        ! Longitud del estanque
     nx=M
     ny = M                                      ! Numero de nodos
     dx=Lx/(nx-1)
     dy= Ly/(ny-1)
     nt = 500                                    ! Numero de iteraciones para calcular las velocidades
     nit = 50                                    ! Numero de iteraciones para calcular la presion
     dt = 0.01                                   ! Paso temporal
     mu = 0.1                                    ! Viscosidad
     rho = 1.0                                   ! Densidad
     nu = mu/ rho                                ! Viscosidad cinematica

     do i = 1,M
        do j = 1,M
           ux(i, j)=0.0
           v(i, j)=0.0
           U(i, j)=0.0
           P(i, j)=0.0
        end do
     end do
     ! -----------------------------------------------------------------------------------------------------------------------------

     DO t = 1 , nt
        un = ux
        vn = v
   
     ! Condiciones de frontera para la presi¢n
        DO i =1,ny
          P(i,nx) = P(i,nx -1)           ! dP/ dx = 0 en x = 2
          P(1,i) = P(2,i)                 ! dP/dy = 0 en y = 0
          P(i,1) = P(i,2)                 ! dP/ dx = 0 en x = 0
          P(ny,i) = 0                       ! P = 0 en y = 2
        END DO
   
     ! Se calcula la presi¢n en cada nodo de la malla
        DO k = 1 , n i t
           DO i =2,ny-1
              DO j =2,nx-1
                 Pn = P
                 P(i,j) = (((Pn(i, j+1) + Pn(i, j-1))*dy**2 + &
                 (Pn( i +1, j ) + Pn( i-1,j))*dx**2) / &
                 ( 2*( dx**2 + dy**2)) - &
                 dx **2 * dy **2 / (2*( dx**2 + dy**2)) * &
                 (rho*(1/dt * &
                 (( ux (i, j+1) - ux(i, j-1)) / &
                 (2*dx) + (v(i +1,j) - v(i-1, j)) /(2*dy))&
                 -((ux (i,j+1)- ux (i,j-1)) / (2*dx)) * *2 - &
                 2*((ux (i+1, j) - ux (i-1,j)) / (2 * dy ) * &
                 ( v(i , j+1) - v (i ,j-1)) / (2*dx ))- &
                 (( v(i+1, j) - v (i-1, j )) / (2*dy )) * * 2 ) ) )
              END DO
           END DO
        END DO

        ! Se calcula la velocidad ux y uy en cada nodo de la malla .
        DO i =2,ny-1
           DO j =2,nx-1

           ! Calculo de velocidad en x
              ux (i,j) = ( un(i,j)- &
              un (i,j) * dt/dx * &
              (un (i,j) - un (i,j-1)) - &
              vn (i,j) * dt/dy * &
              (un (i,j) - un ( i -1, j ) ) - &
              dt / (2 *rho*dx ) * (P(i,j+1) - P(i,j-1))+&
              nu * ( dt / dx **2 * &
              (un (i,j+1) - 2 * un(i,j) + un(i,j-1)) + &
              dt / dy **2 * &
              (un (i+1,j) - 2 * un(i,j) + un(i-1,j))))

            ! Calculo de velocidad en y
              v(i,j) = ( vn(i,j) - &
              un(i,j) * dt/dx * &
              (vn (i,j) - vn ( i , j -1)) - &
              vn (i,j) * dt / dy * &
              (vn (i,j) - vn (i-1,j)) - &
              dt / (2*rho*dy ) * (P(i+1,j) - P(i-1,j)) +&
              nu * (dt/dx**2 * &
              (vn (i,j+1) - 2 * vn(i,j) + vn(i,j-1)) + &
              dt / dy **2 * &(vn (i+1, j) - 2 * vn (i,j) + vn(i-1,j))))
           END DO
        END DO

        ! Condiciones de frontera para la velocidad
        DO i =1,ny
           v (1,i) = 0.0 ; v(ny,i) = 0.0 ; v(i,1) = 0.0 ; v(i,nx) = 0
           ux (1,i) = 0.0 : ux (i,1) = 0.0 ; ux(i,nx ) = 0
           ux (ny,i) = 1 ! velocidad de la placa en movimiento
        END DO
     END DO

     ! Velocidad total
     DO i =1,ny
        DO j =1,nx
           U(i,j) = sqrt(ux(i,j)**2 + v(i,j)**2)
        END DO
     END DO

     ! Guardamos datos para graficaci¢n
     OPEN(2 , file = 'resultados_U.txt', status='unknown')
     DO i = 1 , ny
        DO j = 1 , nx
           write(2,*)(j-1)*dx , (i-1)*dy , U(i,j) , ux(i ,j) , v(i,j)
        END DO
     END DO
     CLOSE(2)

     OPEN(3 , 'resultados_ux.txt', status='unknown')
     DO i = 1 , ny
        DO j = 1 , nx
           write(3,*)(j-1)*dx , (i-1)*dy , ux(i,j) , ux(i,j) , v(i,j)
        END DO
     END DO
     CLOSE(3)

     OPEN(4 , file = 'resultados_v.txt', status='unknown')
     DO i = 1 , ny
        DO j = 1 , nx
           write(4,*)(j-1)*dx , (i-1)*dy , v(i,j) , ux(i,j) , v(i,j)
        END DO
     END DO
     CLOSE(4)

     OPEN(5 , 'resultados_P.txt', status='unknown')
     DO i = 1 , ny
        DO j = 1 , nx
           write(5,*)(j-1)*dx ,(i-1)*dy , P(i,j)
        END DO
     END DO
     CLOSE(5)

end program estanque