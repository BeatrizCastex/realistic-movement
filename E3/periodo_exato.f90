! Código que calcula o período exato de um pêndulo para qualquer teta

	program periodo_exato
	IMPLICIT NONE

	! Variáveis
	real*8 :: m ! massa
	real*8 :: l ! comprimento do fio
	real*8 :: teta_max ! maior ângulo que o pendulo pode ter
	real*8 :: teta ! angulo do pendulo
	real*8 :: t ! tempo
	real*8 :: dt ! variação do tempo
	real*8 :: w ! velocidade angulas
	real*8 :: Em ! Energia mecânica
	real*8 :: Po ! Período inicial
	real*8 :: P ! Período
	real*8 :: I ! integral
	real*8 :: tetai
	real, parameter :: g = 9.8d0
	real, parameter :: pi = 3.1415927
	real*8 :: h
	real*8, external :: f !função



	! Inicializando variáveis
	m = 1.d0 !kg, dado
	l = 1.d0 !m, dado
	teta_max = pi/2.d0 ! dado
	Po = 2.00709 ! s, dado
	dt = 0.005d0 !s, dado
	w = 0.d0
	t = 0.d0
	teta = teta_max
	I = Po
	tetai = 0.d0
	
	! Abrindo arquivos para imprimir resultados
	! Arquivo para energia mecânica 
	open(95,file = "periodo_exato.dat")




	! Devemos calcular o teta em um intervalo de tempo 0 <= t <= 40 s
	DO WHILE (t <= 40)


		! Usamos o método de Euler-Cromer para calcular o ângulo e a velocidade angular

		w = w - (sin(teta)*dt*(g/l))
		teta = teta + (w * dt)


		! Usamos a integral do trapézio para calcular o período exato
		I = I + 0.5*(f(teta+tetai)*dt)
		P = 4*sqrt(l/g)*I
		tetai = teta
		t = t + dt

		write(95,*)teta,P


	enddo

	end program periodo_exato

	real*8 function f(p)
	IMPLICIT NONE
	
	!Variáveis:
	real, parameter :: pi = 3.1415927
	real*8, intent(in)  :: p !Variável de entrada

	f = 1/sqrt(1 - ((sin((pi/2)/2)**2)*((sin(p))**2)))
	
	RETURN 
	end