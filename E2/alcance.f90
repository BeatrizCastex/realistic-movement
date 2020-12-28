	!Programa que calcula o alcance de um projétil
	
	program alcance
	IMPLICIT NONE

	!Variáveis
	real*8 :: vo !velocidade inicial
	real*8 :: vx !velocidade em x
	real*8 :: vy !velocidade em y
	real*8 :: y !posição y
	real*8 :: x !posição x
	real*8 :: xi !posição x no passo anterior
	real*8 :: teta !angulo da velocidade com o eixo x
	real*8 :: teta_max !ângulo que maximiza o alcance
	real*8 :: m !massa
	real*8 :: g !gravidade
	real*8 :: r !densidade
	real*8 :: A !área
	real*8 :: CD !drag coefficient
	real*8 :: dt !Vaiação de tempo
	real*8 :: cal !calibre
	real*8 :: aux ! auxiliar
	integer*8 :: c !contador
	real, parameter :: pi = 3.1415927

	!Inicalizando variáveis:
	vo = 377.d0 !m/s, dado
	m = 42.d0 !kg, dado
	cal = 149.1d0 !mm, dado
	g = 9.8d0 !m/s², dado
	r = 1.3d0 !kg/m³, dado
	CD = 0.295d0 !dado
	dt = 0.01d0 !dado
	xi = 0.d0

	
	!Calculando a área:
	A = pi*(((cal*0.001d0)/2.d0)**2.d0) !pi*r²

	!Abrindo arquivo para impressão de resultados
	open(20,file = "alcance.dat")

	!Precisamos calcular o alcance para vários tetas diferentes, para fazer um gráfico.
	Do c = 0, 10000

		!mantemos teta de 0 à pi/2
		teta = real(c,8)*(pi/2.d0)*1d-4

		!precisamos saber das velocidades para encontrarmos as posições, temos o vetor inicial vo, que separamos em componentes x e y
		vx = vo*(cos(teta))
		vy = vo*(sin(teta))
		x = 0.d0
		y = 0.d0

		!Seguimos o método de Euler para calcular o alcance e a altura para cada teta

		DO WHILE(y >= 0)
			
			aux = vx
		
			vx = vx - (((r*A*CD)/(2.d0*m))*sqrt((aux**2.d0)+(vy**2.d0))*aux*dt)
		
			vy = vy - ((g+(((r*A*CD)/(2.d0*m))*sqrt((aux**2.d0)+(vy**2.d0))*vy))*dt)

			x = x + vx*dt
			y = y + vy*dt

			!Vendo se cada alcance é maior do que o anterior, nós descobrimos o angulo de alcance máximo.

			IF (x .GE. xi) THEN
				xi = x
				teta_max = teta
			endif

		enddo

		write(20,*)teta,x
	enddo
		write(20,*)teta_max

	end program alcance
