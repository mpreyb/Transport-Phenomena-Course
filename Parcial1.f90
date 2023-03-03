program tareauno

    ! Nombre de archivo: kirchhoff_diffraction.py
    ! Autor: Maria Paula Rey
    ! Email : mpreyb@eafit.edu.co
    ! Ultima modificacion: 05/03/2022
    ! Version Force: 2.0

    !----------------------------------- Variables y constantes ---------------------------------------------------
    ! El subindice 1 corresponde al agua, 2 corresponde al aceite

    integer,parameter :: nmi=100000, nd=40
    ! mni = numero maximo de iteraciones.
    ! nd = numero de divisiones en cada fase

    integer :: i,icont,iter
    real :: temp, S, Dy1, Dy2, mu1, mu2, yint,yO, gpre,ytemp
    real,dimension(2*nd+1)::u,dely

    naprox=9     !Orden de aproximacion
    mu1=10       !Viscosidad del agua [Pa*s]
    mu2=30       !Viscosidad del aceite [Pa*s]
    yO=0.5       !Distancia entre placas [m]
    yint=0.1     !Altura de la interfase [m]
    gpre=500     !Gradiente de presion al que se somete el fluido [Pa*s/m]
    !----------------------------------------------------------------------------------------------------------------

    aprox=10.0**(- naprox)
    ! Se calcula el diferencial de altura en cada una de las fases dependiendo del numero de divisiones (nd).
    Dy1=yint/real(nd)
    Dy2=(yO-yint)/real(nd)

    !----------------------------------- Condiciones de frontera ---------------------------------------------------
    u(1)=0.0        !Velocidad de la placa inferior
    u(2*nd+1)=0.0   !Velocidad de la placa superior
    dely(1)=0.0     !Altura inicial (de la placa inferior)
    dely(2*nd+1)=yO !Altura de la placa superior
    !----------------------------------------------------------------------------------------------------------------


    !Inicializacion del vector perfil de velocidad
    do i=2,2*nd
        u(i)=0.0
    end do

    !Iniciacion de las iteraciones
    do iter=1,nmi-1
        icont=0                               !bandera de deteccion para interrumpir la operacion cuando haya convergencia
        do i=2,2*nd
            !En los dos condicionales siguientes se calcula la altura de cada nodo, y el termino S independiente de cada fase.
            if (i .gt. nd) then               !Constantes del modelo para el aceite
                ytemp=yint+(i-nd-1)*Dy2
                S=(Dy2**2)*gpre/mu2
            end if
            if (i.le.nd) then                 !Constantes del modelo para el agua
                ytemp=(i-1)*Dy1
                S=(Dy1**2)*gpre/mu1
            end if

            dely(i)=ytemp                     !Se llena el vector de alturas
            temp=( u(i-1) + u(i+1)+S)*0.5     !Se calcula el valor temporal de velocidad
            if(  (abs(u(i)-temp)/temp) .le. aprox) go to 2
                icont=icont+1
            2 u(i)=temp
        end do
    end do


    200 format(E10.3,2X,F6.5)                   !Formato para escribir los datos en el archivo de soluciones
    open(57,file='perfil.dat',status='unknown') !Se abre el archivo con las soluciones.

    !Se escriben los resultados calculados en un archivo de nombre 'perfil.dat'
    do i=1,2*nd+1
        write(*,*) u(i), dely(i)
        write(57,200) u(i), dely(i)
    end do
    close(57)
    write(*,*) 'A continuacion use el archivo "graficando.m" para visualizar los datos'
    read(*,*)
end program tareauno
