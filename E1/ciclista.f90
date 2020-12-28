	!Programa que calcula a velocidade de um ciclista pelo método de euler e compará-a com a velocidade real

	program ciclista
	IMPLICIT NONE

	!Variáveis
	real*8 :: v    ! Velocidade que procuramos
	!real*8 :: ve  ! velocidade exata
	real*8 :: t    ! tempo total decorrido
	real*8 :: dt   ! variação do tempo
	real*8 :: m    ! massa do ciclista
	real*8 :: r    ! densidade
	real*8 :: A    ! área
	real*8 :: P    ! Potência
	real*8 :: K    ! energia cinética
	real*8 :: n    ! número de iterações
	real*8 :: S    ! Espaço (distância)
	real*8 :: vi   ! Velocidade anterior 
	integer*8 :: i ! contador


	! abrindo um arquivo para imprimir resultados
	!open(80,file = "ciclista_resistente.dat")
	!open(20, file = "ciclista_exato.dat")
	!open(50, file = "ciclista_Areas3.dat")
	! Arquivo para salvar a distância percorrida
	open(97,file = "distancia_res.dat")

	
	! Definindo constantes:
	m = 70.d0 !dado (kg)
	P = 400.d0 !dado (W)
	t = 300.d0 ! 5 mins = 300 s
	dt = 0.01d0 !dado (s)
	r = 1.3d0   !dado
	A = 0.3d0   !dado
	vi = 0.d0 ! Consideramos a velocidade anterior como zero para corrigirmos a integração do espaço
	S = 0.d0

	! Decidindo o número de iterações:
	! O número de iterações será o tempo total dividido pelos pequenos intervalos de tempo dt
	n = t/dt

	! Escolherndo a velocidade inicial:
	! A velocidade inical não pode ser 0 pela natureza do método, portanto escolhemos números próximos de 0
	!v = 1.d0*(10.d0**(-4.d0))
	! Porém se estivermos fazendo por energia cinética, podemos considerar v = 0
	! como K = 1/2(mv²) => K = 1/2(m)
	K = 0.5d0*m


	! Descobrindo a velocidade:
	Do i=1,int(n,8)
	
		! Realizar o cálculo do método de euler
		
		! Por velocidade
		!v = v + ((P/(m*v))*dt) - (((r*A*(v**2.d0))/2*m)*dt)

		! Por energia cinética
		K = K + dt*(P - (sqrt(2.d0/(m**3.d0))*r*A*(K**(1.5d0))))

		v = sqrt((2.d0*K)/m)

		! Valor exato
		!ve = sqrt((2.d0*(P*(dt*i)))/m)

		! A distância percorrida precisa ser calculada por integração
		! por motivos de simplicidade utilizaremos a integral de trapézio, utilizando o próprio valor de dt como h
		S = S + 0.5d0*(dt*(v+vi))

		vi = v

		!write(80,*)(dt*i),v
		!write(20,*)(dt*i),ve
		!write(50,*)(dt*i),v


	enddo

		write(97,*)S
	
	
	end program ciclista
