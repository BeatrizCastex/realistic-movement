	! Programa que calcula o período de um pendulo em função de teta-max pelo metodo de verlet

	program periodo_verlet
	IMPLICIT NONE

	! Variáveis
	real*8 :: m ! massa
	real*8 :: l ! comprimento do fio
	real*8 :: teta_max ! maior ângulo que o pendulo pode ter
	real*8 :: teta ! angulo do pendulo
	real*8 :: t ! tempo
	real*8 :: dt ! variação do tempo
	real*8 :: Em ! Energia mecânica
	real*8 :: Po ! Período inicial
	real*8 :: P ! Período
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
	t = 0.d0
	teta = teta_max
	tetai = teta
	
	! Abrindo arquivos para imprimir resultados
	! Arquivo para energia mecânica 
	open(99,file = "periodo_verlet.dat")




	! Devemos calcular o teta em um intervalo de tempo 0 <= t <= 40 s
	DO WHILE (t <= 40)


		! Usamos o método de verlet para calcular o ângulo
		tetai_1 = tetai
		tetai = teta
		teta = ((2.d0*tetai)-tetai_1)-((g/l)*(sin(tetai)*(dt*dt)))


		! o período pode ser dado através de uma aproximação encontrada em:
		! http://users.df.uba.ar/sgil/physics_paper_doc/papers_phys/mechan/Pendulo2.pdf
		P = 2.d0*pi*sqrt((l/g)*cos(teta/2.d0))

		t = t + dt

		write(99,*)teta,P


	enddo

	end program periodo_verlet
