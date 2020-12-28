	! Programa que calcula o angulo e velocidade de um pêndulo pelo método de verlet

	program pendulo_graph_verlet
	IMPLICIT NONE

	! Variáveis
	real*8 :: m ! massa
	real*8 :: l ! comprimento do fio
	real*8 :: teta_max ! maior ângulo que o pendulo pode ter
	real*8 :: teta ! angulo do pendulo
	real*8 :: t ! tempo
	real*8 :: dt ! variação do tempo
	real*8 :: w ! velocidade angulas
	real*8 :: Po ! Período inicial
	real*8 :: Em !Energia mecanica
	real*8 :: Ep ! Energia potencial
	real*8 :: Ek ! Energia cinética
	real*8 :: w_exato
	real*8 :: tetai_1
	real*8 :: tetai
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
	teta = pi/2
	
	! Abrindo arquivos para imprimir resultados
	! Arquivo para valores aproximados
	open(56,file = "pendulo_graph_verlet.dat")
	! Arquivo para valores exatos
	open(82,file = "pendulo_graph_verlet_exato.dat")

	!Para o ângulo e velocidade angular teta1 e w1 precisamos usar o método de Euler:
	tetai = teta

	
	DO WHILE (t <= 40)


		! Usamos o método de Verlet para calcular o ângulo e a velocidade angular
		tetai_1 = tetai
		tetai = teta
		teta = ((2.d0*tetai)-tetai_1)-((g/l)*(sin(tetai)*(dt*dt)))
		w = w - (sin(tetai)*dt*(g/l))


		! Valor exato
		! Para o valor exato da velocidade primeiro calculamos a energia mecânica, que se conserva:
		Em = m*g*l

		! Calculamos então a energia potêncial no tempo t
		Ep = m*g*(l*(1-cos(teta)))

		! Calculamos a energia cinética
		Ek = Em - Ep

		!Pela energia cinética calculamos a velocidade angular
		w_exato = sqrt((2*Ek)/m)/l

		t = t + dt

		write(56,*)teta,abs(w)
		write(82,*)teta,w_exato




	enddo




	end program pendulo_graph_verlet 
