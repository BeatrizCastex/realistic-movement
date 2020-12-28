	!Programa que calcula a trajetória de um projétil
	
	program projetil
	IMPLICIT NONE

	!Variáveis
	real*8 :: vo !velocidade inicial
	real*8 :: vx !velocidade em x
	real*8 :: vy !velocidade em y
	real*8 :: y !posição y
	real*8 :: x !posição x
	real*8 :: teta !angulo da velocidade com o eixo x
	real*8 :: m !massa
	real*8 :: g !gravidade
	real*8 :: r !densidade
	real*8 :: A !área
	real*8 :: CD !drag coefficient
	real*8 :: dt !Variação de tempo
	real*8 :: aux ! auxiliar
	real*8 :: cal !calibre
	real, parameter :: pi = 3.1415927

	!Inicalizando variáveis:
	vo = 377.d0 !m/s, dado
	m = 42.d0 !kg, dado
	cal = 149.1d0 !mm, dado
	g = 9.8d0 !m/s², dado
	r = 0.d0 !kg/m³, dado
	CD = 0.295d0 !dado
	dt = 0.01d0 !dado
	teta = pi/3.d0
	vx = vo*(cos(teta))
	vy = vo*(sin(teta))
	x = 0.d0
	y = 0.d0
	
	!Calculando a área:
	A = pi*(((cal*0.001d0)/2.d0)**2.d0) !pi*r²

	!Abrindo arquivo para impressão de sesultados
	open(19,file = "projetil_6.dat")
	write(19,*)"teta: ", teta," densidade: ", r

	!Calculamos o movimento de lançamento assumindo que o objeto começa na mesma altura que termina, ou seja, zero.

	DO WHILE(y >= 0)

		!Utilizamos o método de Euler para calcular a haltura e alcance do objeto a cada pedaço de tempo, e imprimimos a trajetória.

		aux = vx
	
		vx = vx - (((r*A*CD)/(2*m))*sqrt((vx**2)+(vy**2))*vx*dt)
		
		vy = vy - ((g+(((r*A*CD)/(2*m))*sqrt((aux**2)+(vy**2))*vy))*dt)

		x = x + vx*dt
		y = y + vy*dt

		write(19,*)x,y
	enddo
	
	end program projetil

