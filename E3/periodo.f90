	! Programa que calcula o período de um pendulo em função de teta-max pelo metodo de euler-cromer

	program periodo
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
	open(92,file = "periodo.dat")




	! Devemos calcular o teta em um intervalo de tempo 0 <= t <= 40 s
	DO WHILE (t <= 40)


		! Usamos o método de Euler-Cromer para calcular o ângulo e a velocidade angular

		w = w - (sin(teta)*dt*(g/l))
		teta = teta + (w * dt)


		! o período pode ser dado através de uma aproximação encontrada em:
		! http://users.df.uba.ar/sgil/physics_paper_doc/papers_phys/mechan/Pendulo2.pdf
		P = 2.d0*pi*sqrt((l/g)*cos(teta/2.d0))

		t = t + dt

		write(92,*)teta,P


	enddo

	end program periodo