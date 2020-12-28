	!Programa que calcula a o ângulo e energia mecânica de um pêndulo pelo método de Euler

	program pendulo_a
	IMPLICIT NONE

	! Variáveis
	real*8 :: m ! massa
	real*8 :: l ! comprimento do fio
	real*8 :: teta_max ! maior ângulo que o pendulo pode ter
	real*8 :: teta ! angulo do pendulo
	real*8 :: t ! tempo
	real*8 :: dt ! variação do tempo
	real*8 :: w ! velocidade angular
	real*8 :: Em ! Energia mecânica
	real*8 :: Po ! Período inicial
	real*8 :: aux
	real, parameter :: g = 9.8d0
	real, parameter :: pi = 3.1415927


	! Inicializando variáveis
	m = 1.d0 !kg, dado
	l = 1.d0 !m, dado
	teta_max = pi/2.d0 ! dado
	Po = 2.00709 ! s, dado
	dt = 0.005d0 !s, dado
	w = 0.d0
	t = 0.d0
	teta = teta_max
	
	! Abrindo arquivos para imprimir resultados
	! Arquivo para energia mecânica 
	open(42,file = "pendulo_Em_a.dat")
	! Arquivo para ângulo
	open(46,file = "pendulo_teta_a.dat")




	! Devemos calcular o teta e a energia mecânica em um intervalo de tempo 0 <= t <= 40 s
	DO WHILE (t <= 40)


		! Usamos o método de Euler para calcular o ângulo e a velocidade angular
		aux = teta
		teta = teta + (w * dt)
		w = w - (sin(aux)*dt*(g/l))


		! A energia mecanica é a energia potencial mais a energia cinética 
		Em = (m*g*(L*(1 - cos(teta)))) + (0.5d0*m*((w*l)**2.d0))

		t = t + dt

		write(42,*)t,Em
		write(46,*)t,teta


	enddo




	end program pendulo_a