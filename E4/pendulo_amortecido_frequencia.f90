! Programa que calcula o menor teta em função da frequência de um pêndulo amortecido e influênciado por uma forç externa na vertival pelo método de Verlet

	program pendulo_amortecido_frequencia
	IMPLICIT NONE

	! Variáveis:
	real*8 :: alpha ! amplitude do movimento
	real*8 :: beta  ! coef. de Stokes
	real*8 :: t ! tempo adimensional
	real*8 :: dt ! variação do tempo
	integer*8 :: freq ! frequencia do movimento vertical
	real*8 :: w ! velocidade angular
	real*8 :: teta ! angulo
	real*8 :: a ! aceleração
	real*8 :: teta_min
	real*8 :: tetai_1
	real*8 :: tetai
	real*8 :: aux
	real, parameter :: g = 9.8d0
	real, parameter :: pi = 3.1415927

	! Inicializando variáveis:
	alpha = 0.1d0
	beta = 0.05d0
	dt = 0.0001d0

	! Abrindo arquivos para imprimir resultados
	open(87,file = "pendulo_amortecido_frequencia.dat")


	!Para o ângulo e velocidade angular teta1 e w1 precisamos usar o método de Euler:
	tetai = teta

	DO freq = 0,20

		teta = 0.9999d0*pi
		w = 0.d0


		!Para o ângulo e velocidade angular teta1 e w1 precisamos usar o método de Euler:
		tetai = teta

		! Devemos calcular o teta em um certo intervalo de tempo
		DO WHILE (t <= 30.d0*2.d0*pi)

			! Vamos utilizar o método de Verlet para calcular o movimento
			! Começamos com uma aproximação da aceleração
			tetai_1 = tetai
			tetai = teta
			a = -(beta * w) - ((1.d0 - (alpha * (freq * freq) * sin(freq * t))) * sin(tetai))
			teta = ((2.d0*tetai)-tetai_1)+(a * (dt * dt))
			w = (tetai - tetai_1)/dt

			IF (tetai > teta) THEN
				teta_min = teta
			endif

			t = t+dt

		enddo

		write(87,*)freq,teta_min

	enddo

	end program pendulo_amortecido_frequencia